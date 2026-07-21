# -*- coding: utf-8 -*-
"""Distribution diagnostics for a generated dimer ensemble.

Pure, framework-free analysis intended to be shared by two consumers:

  * the SEAMM ``dimer_builder_step`` -- called from ``analyze()`` to write the
    charts automatically whenever the step runs; and
  * a standalone notebook script -- to re-analyze an exported ensemble offline.

Keeping the math in one place means the in-step charts and the offline analysis
never drift.  The module depends only on numpy for the metrics; **plotly** is
imported lazily inside ``make_dashboard``, so ``compute_metrics`` works in a
minimal / headless environment.  ``make_dashboard`` returns the combined panel
as one plotly ``go.Figure``; ``make_panels`` returns the *same* panels as
individual ``go.Figure`` objects (a dict keyed by panel name) for callers that
want each chart on its own -- both share one set of trace builders, so the
combined and separated views can never drift.  SEAMM uses plotly throughout, so
the step writes any of these with the same machinery it uses for every other
graph (and the offline notebook can ``write_image`` / ``write_html`` them).

Design notes
------------
* **General, not motif-specific.**  Every default metric is defined from
  geometry + (optional) interaction energy alone -- separation, approach
  direction, contact distances, relative orientation, and the energy
  distribution.  No hydrogen-bond / donor-acceptor knowledge is assumed, in
  keeping with the energy-stratified sampling design (water-electrolyte plan,
  decision 5).  A motif overlay can be layered on by a caller that knows the
  contact atoms; it is deliberately not baked in here.
* **Degrades gracefully.**  If the ensemble carries no energies (e.g. the
  van-der-Waals contact method was used), the energy panels are simply omitted.
* **Cheap.**  All metrics are O(N * nA * nB); fine for the hundreds-to-thousands
  of dimers a build produces.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = [
    "Dimer",
    "DimerMetrics",
    "compute_metrics",
    "energy_histogram",
    "approach_concentration",
    "summarize",
    "make_dashboard",
    "make_panels",
    "panel_separation",
    "panel_contact",
    "panel_approach",
    "panel_orientation",
    "panel_energy",
    "panel_energy_vs_R",
]

# --------------------------------------------------------------------------- #
# Minimal element data (extend as needed).  Masses in u; Bondi vdW radii in Å.
# --------------------------------------------------------------------------- #
_MASS = {
    "H": 1.008, "Li": 6.941, "B": 10.811, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Na": 22.990, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "K": 39.098,
}
_VDW = {
    "H": 1.20, "Li": 1.82, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "Na": 2.27, "P": 1.80, "S": 1.80, "Cl": 1.75, "K": 2.75,
}


def _masses(symbols) -> np.ndarray:
    return np.array([_MASS.get(s, 12.0) for s in symbols], dtype=float)


# --------------------------------------------------------------------------- #
# Ensemble record
# --------------------------------------------------------------------------- #
@dataclass
class Dimer:
    """One generated dimer.

    Coordinates are (n, 3) arrays in Å; ``energy`` is the (surrogate)
    interaction energy in kJ/mol, or ``None`` if not available.
    """

    symbols_A: list
    xyz_A: np.ndarray
    symbols_B: list
    xyz_B: np.ndarray
    energy: Optional[float] = None
    label: Optional[str] = None
    # scan variables the builder may already have saved as properties:
    separation: Optional[float] = None   # nominal COM separation, Å
    orientation: Optional[int] = None    # orientation-seed index

    def __post_init__(self):
        self.xyz_A = np.asarray(self.xyz_A, dtype=float).reshape(-1, 3)
        self.xyz_B = np.asarray(self.xyz_B, dtype=float).reshape(-1, 3)


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def _com(masses, xyz) -> np.ndarray:
    return (masses[:, None] * xyz).sum(0) / masses.sum()


def _principal_frame(masses, xyz):
    """Return (com, axes) where ``axes`` columns are the mass-weighted principal
    axes as a proper rotation (det = +1).  Rotating a lab vector by ``axes.T``
    expresses it in the molecule's principal-axis frame.

    For a single atom (no orientation) the identity frame is returned.
    """
    com = _com(masses, xyz)
    if xyz.shape[0] < 2:
        return com, np.eye(3)
    r = xyz - com
    # inertia tensor
    inertia = np.zeros((3, 3))
    for m, (x, y, z) in zip(masses, r):
        inertia[0, 0] += m * (y * y + z * z)
        inertia[1, 1] += m * (x * x + z * z)
        inertia[2, 2] += m * (x * x + y * y)
        inertia[0, 1] -= m * x * y
        inertia[0, 2] -= m * x * z
        inertia[1, 2] -= m * y * z
    inertia[1, 0] = inertia[0, 1]
    inertia[2, 0] = inertia[0, 2]
    inertia[2, 1] = inertia[1, 2]
    _, axes = np.linalg.eigh(inertia)  # ascending eigenvalues, orthonormal cols
    if np.linalg.det(axes) < 0:        # keep a proper rotation, not a reflection
        axes[:, 0] = -axes[:, 0]
    return com, axes


def _is_orientable(xyz) -> bool:
    """True if the fragment defines a usable first principal axis (>=2 atoms,
    not effectively a single point)."""
    if xyz.shape[0] < 2:
        return False
    return np.ptp(xyz, axis=0).max() > 1e-6


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
@dataclass
class DimerMetrics:
    n: int
    R: np.ndarray                       # (n,) COM-COM separation, Å
    min_contact: np.ndarray            # (n,) closest A-B atom pair, Å
    approach_lonlat: np.ndarray        # (n, 2) lon,lat (rad) of B, in A's frame
    approach_vec: np.ndarray           # (n, 3) unit approach vectors in A's frame
    orient_angle: np.ndarray           # (n,) angle between A,B principal axes, deg
    energy: Optional[np.ndarray]       # (n,) interaction energy, kJ/mol, or None
    contact_pairs: dict = field(default_factory=dict)  # (el,el)->distances array

    @property
    def has_energy(self) -> bool:
        return self.energy is not None and np.isfinite(self.energy).any()


def compute_metrics(ensemble, contact_cutoff: float = 6.0) -> DimerMetrics:
    """Compute distribution metrics for a list of :class:`Dimer`.

    Parameters
    ----------
    ensemble : list[Dimer]
    contact_cutoff : float
        Only intermolecular A-B distances below this (Å) are pooled into the
        contact-RDF; the per-dimer minimum contact is always exact.
    """
    n = len(ensemble)
    R = np.full(n, np.nan)
    min_contact = np.full(n, np.nan)
    lonlat = np.full((n, 2), np.nan)
    avec = np.full((n, 3), np.nan)
    orient = np.full(n, np.nan)
    energies = np.full(n, np.nan)
    have_E = False
    pairs: dict = {}

    for i, d in enumerate(ensemble):
        mA, mB = _masses(d.symbols_A), _masses(d.symbols_B)
        comA, axesA = _principal_frame(mA, d.xyz_A)
        comB, axesB = _principal_frame(mB, d.xyz_B)

        vec = comB - comA
        r = float(np.linalg.norm(vec))
        R[i] = r

        # approach direction, expressed in A's principal-axis frame
        if r > 1e-9:
            u = axesA.T @ (vec / r)
            avec[i] = u
            lonlat[i, 0] = np.arctan2(u[1], u[0])
            lonlat[i, 1] = np.arcsin(np.clip(u[2], -1.0, 1.0))

        # relative orientation: angle between first principal axes (sign-folded)
        if _is_orientable(d.xyz_A) and _is_orientable(d.xyz_B):
            c = abs(float(np.dot(axesA[:, 0], axesB[:, 0])))
            orient[i] = np.degrees(np.arccos(np.clip(c, 0.0, 1.0)))

        # intermolecular contacts
        diff = d.xyz_A[:, None, :] - d.xyz_B[None, :, :]
        dist = np.linalg.norm(diff, axis=2)          # (nA, nB)
        min_contact[i] = dist.min()
        near = dist < contact_cutoff
        for ia, sa in enumerate(d.symbols_A):
            for ib, sb in enumerate(d.symbols_B):
                if near[ia, ib]:
                    key = tuple(sorted((sa, sb)))
                    pairs.setdefault(key, []).append(dist[ia, ib])

        if d.energy is not None and np.isfinite(d.energy):
            energies[i] = float(d.energy)
            have_E = True

    contact_pairs = {k: np.asarray(v) for k, v in pairs.items()}
    return DimerMetrics(
        n=n,
        R=R,
        min_contact=min_contact,
        approach_lonlat=lonlat,
        approach_vec=avec,
        orient_angle=orient,
        energy=energies if have_E else None,
        contact_pairs=contact_pairs,
    )


def energy_histogram(energy, nbins: int = 24, clip_high: Optional[float] = None):
    """Histogram of interaction energies plus a flatness metric.

    Returns (edges, counts, flatness) where ``flatness`` is the coefficient of
    variation of the non-empty bin counts -- 0 is perfectly flat (ideal for the
    energy-stratified sampler); large values mean the samples pile up in a few
    energy bins (what uniform / Boltzmann sampling does).
    """
    e = np.asarray(energy)
    e = e[np.isfinite(e)]
    if clip_high is not None:
        e = e[e <= clip_high]
    if e.size == 0:
        return np.array([0.0, 1.0]), np.array([0]), 0.0
    edges = np.linspace(e.min(), e.max(), nbins + 1)
    counts, _ = np.histogram(e, edges)
    nz = counts[counts > 0].astype(float)
    flatness = float(nz.std() / nz.mean()) if nz.size > 1 else 0.0
    return edges, counts, flatness


def approach_concentration(approach_vec) -> float:
    """Mean resultant length of the approach unit vectors (in A's frame).

    0 = isotropic (uniform orientation), 1 = all partners approach from one
    direction.  A simple scalar for "is the orientational sampling biased?".
    (Multi-lobe patterns -- e.g. several symmetry-equivalent H-bond directions --
    partially cancel, so read this together with the sphere map.)
    """
    v = np.asarray(approach_vec)
    v = v[np.isfinite(v).all(axis=1)]
    if v.size == 0:
        return 0.0
    return float(np.linalg.norm(v.mean(axis=0)))


def summarize(metrics: DimerMetrics) -> dict:
    """A compact dict of scalar diagnostics (for logging / assertions / the
    step's text summary)."""
    m = metrics
    out = {
        "n": m.n,
        "R_min": float(np.nanmin(m.R)),
        "R_max": float(np.nanmax(m.R)),
        "min_contact_p05": float(np.nanpercentile(m.min_contact, 5)),
        "min_contact_median": float(np.nanmedian(m.min_contact)),
        "approach_concentration": approach_concentration(m.approach_vec),
    }
    if m.has_energy:
        e = m.energy[np.isfinite(m.energy)]
        _, _, flat = energy_histogram(e)
        out.update(
            energy_min=float(e.min()),
            energy_mean=float(e.mean()),
            energy_frac_attractive=float((e < 0).mean()),
            energy_flatness=flat,
        )
    return out


# --------------------------------------------------------------------------- #
# Plotting (plotly imported lazily) -- returns a go.Figure
# --------------------------------------------------------------------------- #
def _mollweide_xy(lon, lat):
    """Project (lon, lat) in radians to Mollweide x, y (equal-area).

    Boundary is the ellipse x in [-2, 2], y in [-1, 1] (scaled by 1/sqrt2).
    """
    lat = np.clip(lat, -np.pi / 2 + 1e-6, np.pi / 2 - 1e-6)
    theta = lat.copy()
    for _ in range(6):  # Newton solve 2θ + sin2θ = π sinφ
        theta -= (2 * theta + np.sin(2 * theta) - np.pi * np.sin(lat)) / (
            2 + 2 * np.cos(2 * theta)
        )
    x = (2 * np.sqrt(2) / np.pi) * lon * np.cos(theta) / np.sqrt(2)
    y = np.sqrt(2) * np.sin(theta) / np.sqrt(2)
    return x, y


# A brand-neutral, light/dark-neutral palette echoing the notebook pages.
_INK, _CC, _DFT, _OK, _BENCH, _MUTED = (
    "#16161a", "#4a3aa7", "#2160c4", "#1f5130", "#c0392b", "#6b6b72",
)


def _axis_id(row, col, ncols):
    """Axis reference (e.g. 'x3') for the subplot at (row, col) in an ncols grid.
    Panel (1,1) is axis 'x'."""
    idx = (row - 1) * ncols + col
    return "x" if idx == 1 else f"x{idx}"


# --- per-panel trace builders: the single source shared by the combined ----- #
#     dashboard and the standalone panels, so the two views cannot drift.      #
def _add_separation(fig, m, row, col):
    import plotly.graph_objects as go
    R = m.R[np.isfinite(m.R)]
    fig.add_trace(go.Histogram(x=R, nbinsx=24, marker_color=_CC,
                               showlegend=False), row=row, col=col)
    fig.update_xaxes(title_text="COM separation R (Å)", row=row, col=col)
    fig.update_yaxes(title_text="count", row=row, col=col)


def _add_contact(fig, m, row, col):
    import plotly.graph_objects as go
    mc = m.min_contact[np.isfinite(m.min_contact)]
    fig.add_trace(go.Histogram(x=mc, nbinsx=24, marker_color=_DFT, opacity=0.5,
                               name="closest contact"), row=row, col=col)
    if m.contact_pairs:
        top = sorted(m.contact_pairs.items(), key=lambda kv: -kv[1].size)[:3]
        lo = min(mc.min(), min(arr.min() for _, arr in top))
        hi = max(mc.max(), max(arr.max() for _, arr in top))
        edges = np.linspace(lo, hi, 31)
        centers = 0.5 * (edges[:-1] + edges[1:])
        for (a, b), arr in top:
            counts, _ = np.histogram(arr, edges)
            fig.add_trace(go.Scatter(x=centers, y=counts, mode="lines",
                                     line_shape="hvh", name=f"{a}–{b} (all)"),
                          row=row, col=col)
    fig.update_xaxes(title_text="intermolecular distance (Å)", row=row, col=col)
    fig.update_yaxes(title_text="count", row=row, col=col)


def _add_approach(fig, m, row, col, ncols, colorbar=None):
    import plotly.graph_objects as go
    good = np.isfinite(m.approach_lonlat).all(axis=1)
    mx, my = _mollweide_xy(m.approach_lonlat[good, 0], m.approach_lonlat[good, 1])
    # bounding ellipse + a few graticule parallels for orientation
    u = np.linspace(0, 2 * np.pi, 181)
    fig.add_trace(go.Scatter(x=2 * np.cos(u), y=np.sin(u), mode="lines",
                             line=dict(color=_MUTED, width=1), hoverinfo="skip",
                             showlegend=False), row=row, col=col)
    for lat0 in (-np.pi / 4, 0.0, np.pi / 4):
        gl = np.linspace(-np.pi, np.pi, 121)
        gx, gy = _mollweide_xy(gl, np.full_like(gl, lat0))
        fig.add_trace(go.Scatter(x=gx, y=gy, mode="lines",
                                 line=dict(color=_MUTED, width=0.4),
                                 hoverinfo="skip", showlegend=False),
                      row=row, col=col)
    if m.has_energy:
        e = m.energy[good]
        marker = dict(size=5, color=e, colorscale="Viridis", showscale=True,
                      colorbar=colorbar or dict(title="ΔE", thickness=12))
    else:
        marker = dict(size=4, color=_CC, opacity=0.6)
    fig.add_trace(go.Scatter(x=mx, y=my, mode="markers", marker=marker,
                             hoverinfo="skip", showlegend=False), row=row, col=col)
    # equal aspect, no ticks -- it's a map, not a plot
    fig.update_xaxes(visible=False, range=[-2.15, 2.15], row=row, col=col)
    fig.update_yaxes(visible=False, range=[-1.1, 1.1],
                     scaleanchor=_axis_id(row, col, ncols), scaleratio=1,
                     row=row, col=col)


def _add_orientation(fig, m, row, col):
    import plotly.graph_objects as go
    oa = m.orient_angle[np.isfinite(m.orient_angle)]
    if oa.size:
        fig.add_trace(go.Histogram(x=oa, xbins=dict(start=0, end=90, size=5),
                                   marker_color=_OK, showlegend=False),
                      row=row, col=col)
        fig.update_xaxes(title_text="axis–axis angle (deg)", range=[0, 90],
                         row=row, col=col)
    else:
        fig.add_annotation(text="orientation n/a<br>(monatomic fragment)",
                           showarrow=False, font=dict(color=_MUTED),
                           xref="x domain", yref="y domain", x=0.5, y=0.5,
                           row=row, col=col)
    fig.update_yaxes(title_text="count", row=row, col=col)


def _add_energy_hist(fig, m, row, col):
    import plotly.graph_objects as go
    e = m.energy[np.isfinite(m.energy)]
    edges, counts, _ = energy_histogram(e)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig.add_trace(go.Bar(x=centers, y=counts, marker_color=_BENCH,
                         showlegend=False), row=row, col=col)
    fig.add_vline(x=0.0, line=dict(color=_INK, width=0.8, dash="dot"),
                  row=row, col=col)
    fig.update_xaxes(title_text="ΔE interaction (kJ/mol)", row=row, col=col)
    fig.update_yaxes(title_text="count", row=row, col=col)


def _add_energy_vs_R(fig, m, row, col, colorbar=None):
    import plotly.graph_objects as go
    gd = np.isfinite(m.R) & np.isfinite(m.energy)
    fig.add_trace(go.Histogram2d(
        x=m.R[gd], y=m.energy[gd], nbinsx=30, nbinsy=30, colorscale="Magma_r",
        colorbar=colorbar or dict(title="count", thickness=12)),
        row=row, col=col)
    fig.add_hline(y=0.0, line=dict(color=_INK, width=0.8, dash="dot"),
                  row=row, col=col)
    fig.update_xaxes(title_text="COM separation R (Å)", row=row, col=col)
    fig.update_yaxes(title_text="ΔE (kJ/mol)", row=row, col=col)


def _subtitle(m: DimerMetrics) -> str:
    s = summarize(m)
    sub = (f"N = {s['n']}   R∈[{s['R_min']:.2f}, {s['R_max']:.2f}] Å   "
           f"approach concentration = {s['approach_concentration']:.2f}")
    if m.has_energy:
        sub += (f"   ΔE_min = {s['energy_min']:.1f}   "
                f"attractive = {100 * s['energy_frac_attractive']:.0f}%")
    return sub


def make_dashboard(metrics: DimerMetrics, title: str = "Dimer sampling"):
    """Assemble all diagnostics as a single combined plotly ``go.Figure``.

    Panels (energy panels dropped automatically when no energies are present):
      separation histogram | contact distances | approach-direction sphere map |
      relative-orientation histogram | energy histogram + flatness | ΔE vs R.

    See :func:`make_panels` for the same panels as individual figures.
    """
    from plotly.subplots import make_subplots

    m = metrics
    if m.has_energy:
        _, _, flat = energy_histogram(m.energy[np.isfinite(m.energy)])
        titles = [
            "Separation coverage", "Contact distribution",
            "Approach direction (A frame)", "Relative orientation",
            f"Energy distribution (flatness CV = {flat:.2f})",
            "Binding-curve envelope",
        ]
        fig = make_subplots(rows=2, cols=3, subplot_titles=titles,
                            horizontal_spacing=0.085, vertical_spacing=0.16)
        ncols = 3
    else:
        titles = ["Separation coverage", "Contact distribution",
                  "Approach direction (A frame)", "Relative orientation"]
        fig = make_subplots(rows=1, cols=4, subplot_titles=titles,
                            horizontal_spacing=0.06)
        ncols = 4

    _add_separation(fig, m, 1, 1)
    _add_contact(fig, m, 1, 2)
    _add_approach(fig, m, 1, 3, ncols,
                  colorbar=(dict(title="ΔE", len=0.42, y=0.80, x=1.005,
                                 thickness=12) if m.has_energy else None))
    ori_rc = (2, 1) if m.has_energy else (1, 4)
    _add_orientation(fig, m, ori_rc[0], ori_rc[1])
    if m.has_energy:
        _add_energy_hist(fig, m, 2, 2)
        _add_energy_vs_R(fig, m, 2, 3,
                         colorbar=dict(title="count", len=0.42, y=0.20, x=1.005,
                                       thickness=12))

    fig.update_layout(
        template="simple_white",
        title=dict(text=f"<b>{title}</b><br><sup>{_subtitle(m)}</sup>", x=0.045,
                   xanchor="left", font=dict(size=17)),
        height=760 if m.has_energy else 400,
        width=1280,
        bargap=0.05,
        margin=dict(l=70, r=90, t=95, b=60),
        legend=dict(font=dict(size=10), yanchor="top", y=1.0, xanchor="right",
                    x=1.0),
    )
    for ann in fig.layout.annotations:  # shrink the subplot titles a touch
        if ann.text in titles:
            ann.font = dict(size=12)
    return fig


# --------------------------------------------------------------------------- #
# Individual panels -- same trace builders, each on its own go.Figure
# --------------------------------------------------------------------------- #
def _single(add_fn, m, title, width=560, height=400, show_legend=False):
    """Wrap one panel builder in a standalone 1x1 figure with a shared header."""
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=1)
    add_fn(fig)
    fig.update_layout(
        template="simple_white",
        title=dict(text=f"<b>{title}</b><br><sup>{_subtitle(m)}</sup>", x=0.03,
                   xanchor="left", font=dict(size=15)),
        width=width, height=height, bargap=0.05,
        margin=dict(l=65, r=80, t=80, b=55),
        showlegend=show_legend,
        legend=dict(font=dict(size=10), yanchor="top", y=0.98, xanchor="right",
                    x=0.99),
    )
    return fig


def panel_separation(m: DimerMetrics, title: str = "Separation coverage"):
    return _single(lambda f: _add_separation(f, m, 1, 1), m, title)


def panel_contact(m: DimerMetrics, title: str = "Contact distribution"):
    return _single(lambda f: _add_contact(f, m, 1, 1), m, title, show_legend=True)


def panel_approach(m: DimerMetrics, title: str = "Approach direction (A frame)"):
    return _single(lambda f: _add_approach(f, m, 1, 1, 1,
                   colorbar=dict(title="ΔE", thickness=12) if m.has_energy else None),
                   m, title, width=560, height=470)


def panel_orientation(m: DimerMetrics, title: str = "Relative orientation"):
    return _single(lambda f: _add_orientation(f, m, 1, 1), m, title)


def panel_energy(m: DimerMetrics, title: Optional[str] = None):
    """Energy-distribution panel, or ``None`` if the ensemble has no energies."""
    if not m.has_energy:
        return None
    if title is None:
        _, _, flat = energy_histogram(m.energy[np.isfinite(m.energy)])
        title = f"Energy distribution (flatness CV = {flat:.2f})"
    return _single(lambda f: _add_energy_hist(f, m, 1, 1), m, title)


def panel_energy_vs_R(m: DimerMetrics, title: str = "Binding-curve envelope"):
    """ΔE-vs-R panel, or ``None`` if the ensemble has no energies."""
    if not m.has_energy:
        return None
    return _single(lambda f: _add_energy_vs_R(f, m, 1, 1,
                   colorbar=dict(title="count", thickness=12)),
                   m, title, width=600)


def make_panels(metrics: DimerMetrics) -> dict:
    """Return each diagnostic as its own ``go.Figure``.

    Keys: ``separation``, ``contact``, ``approach``, ``orientation`` always;
    ``energy`` and ``energy_vs_R`` only when the ensemble carries energies.
    Uses the same trace builders as :func:`make_dashboard`.
    """
    m = metrics
    panels = {
        "separation": panel_separation(m),
        "contact": panel_contact(m),
        "approach": panel_approach(m),
        "orientation": panel_orientation(m),
    }
    if m.has_energy:
        panels["energy"] = panel_energy(m)
        panels["energy_vs_R"] = panel_energy_vs_R(m)
    return panels
