# -*- coding: utf-8 -*-

"""Tests for the headless dimer-building logic in DimerBuilder."""

import math

import numpy as np
import pytest

from molsystem.system_db import SystemDB
from seamm_util import Q_

import dimer_builder_step
from dimer_builder_step.dimer_builder import DimerBuilder, vdw_radii


def test_mdi_method_and_basis():
    """The engine method/basis picked from a model chemistry: MOPAC/xTB get a
    method alone (no basis); ORCA gets the real (un-aliased) keyword + basis."""
    assert DimerBuilder._mdi_method_and_basis(
        {"method": "PM6-ORG"}, {"mdi_method_arg": "PM6-ORG"}
    ) == ("PM6-ORG", None)
    assert DimerBuilder._mdi_method_and_basis(
        {"method": "REVDSD-PBEP86-D4_2021"},
        {"mdi_method_arg": "REVDSD-PBEP86-D4/2021", "mdi_basis_arg": "def2-TZVP"},
    ) == ("REVDSD-PBEP86-D4/2021", "def2-TZVP")
    # Falls back to the parsed method when no mdi_method_arg is given.
    assert DimerBuilder._mdi_method_and_basis({"method": "PM7"}, {}) == ("PM7", None)


def _add_water(configuration):
    r0 = 0.9572
    theta0 = 104.52
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))
    ids = configuration.atoms.append(
        x=[0.0, x, -x], y=[0.0, 0.0, 0.0], z=[0.0, z, z], atno=[8, 1, 1]
    )
    configuration.bonds.append(i=[ids[0], ids[0]], j=[ids[1], ids[2]], bondorder=[1, 1])
    return configuration


_db_counter = [0]


@pytest.fixture()
def db_two_waters():
    _db_counter[0] += 1
    db = SystemDB(filename=f"file:dimer_test_{_db_counter[0]}?mode=memory&cache=shared")
    a = db.create_system(name="A").create_configuration(name="w1")
    _add_water(a)
    b = db.create_system(name="B").create_configuration(name="w1")
    _add_water(b)

    yield db

    db.close()


def _P(**overrides):
    """A full parameter dict with sensible defaults for the build."""
    P = {
        "input mode": "two monomer sets",
        "monomer A": "A",
        "monomer A configurations": "all",
        "monomer A configuration name": "",
        "monomer B": "B",
        "monomer B configurations": "all",
        "monomer B configuration name": "",
        "number of orientations": 5,
        "random seed": "1",
        "contact method": "van der Waals radii",
        "innermost gap": Q_(-0.5, "Å"),
        "maximum separation": Q_(10.0, "Å"),
        "spacing": "geometric",
        "number of separations": 8,
        "separations": "",
        "energy levels": "-De, -De/2, 0, kBT, 5*kBT",
        "sampling temperature": Q_(300.0, "K"),
        "orientation weighting": "reject shallow orientations",
        "minimum well depth": Q_(1.0, "kJ/mol"),
        "system name": "from monomers",
        "configuration name": "orientation,distance",
        "save scan variables as properties": "yes",
        "analysis plots": "none",
    }
    P.update(overrides)
    return P


def _subset_atom_ids(db, conf, name):
    """All atom ids in the named subsets of a configuration (or None)."""
    if not db.templates.exists(name, "general"):
        return None
    template = db.templates.get(name, "general")
    ids = []
    for subset in conf.subsets.get(template):
        ids.extend(subset.atoms.ids)
    return ids


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #


def test_vdw_radii_in_angstrom():
    r = vdw_radii(["O", "H", "C"])
    assert np.allclose(r, [1.52, 1.10, 1.70], atol=0.02)


def test_contact_distance_two_atoms_along_z():
    node = dimer_builder_step.DimerBuilder()
    A = np.array([[0.0, 0.0, 0.0]])
    B = np.array([[0.0, 0.0, 0.0]])
    axis = np.array([0.0, 0.0, 1.0])
    contact = node._contact_distance(A, np.array([1.5]), B, np.array([1.2]), axis)
    assert math.isclose(contact, 2.7, abs_tol=1.0e-9)


def test_contact_distance_offset_atoms():
    """Lateral offset reduces the along-axis contact distance."""
    node = dimer_builder_step.DimerBuilder()
    A = np.array([[0.0, 0.0, 0.0]])
    B = np.array([[1.0, 0.0, 0.0]])  # 1 Å lateral offset
    axis = np.array([0.0, 0.0, 1.0])
    R = 2.7
    contact = node._contact_distance(A, np.array([1.5]), B, np.array([1.2]), axis)
    assert math.isclose(contact, math.sqrt(R**2 - 1.0), abs_tol=1.0e-9)


def test_separation_schedule_geometric_range():
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    contact = 3.0
    d = node._separation_schedule(contact, P)
    assert len(d) == 8
    assert math.isclose(d.min(), contact - 0.5, abs_tol=1.0e-6)
    assert math.isclose(d.max(), 10.0, abs_tol=1.0e-6)
    steps = np.diff(d)
    assert np.all(steps > 0)  # sorted, strictly increasing
    # Steps grow smoothly with separation: no oversized first step.
    assert np.all(np.diff(steps) >= -1.0e-9)


def test_separation_schedule_explicit():
    node = dimer_builder_step.DimerBuilder()
    P = _P(spacing="explicit", separations="0.0, 1.0, 2.0")
    d = node._separation_schedule(5.0, P)
    assert np.allclose(d, [5.0, 6.0, 7.0])


# --------------------------------------------------------------------------- #
# Full build — Mode A
# --------------------------------------------------------------------------- #


def test_build_mode_a_counts_and_atomset(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    rng = np.random.default_rng(1)

    system, stats = node._build(db, P, rng)

    # 5 orientations x 8 separations
    assert stats["n_configurations"] == 40
    assert len(system.configurations) == 40
    # Every configuration is a dimer of 6 atoms ...
    assert all(c.n_atoms == 6 for c in system.configurations)
    # ... and they are all conformers sharing a single atomset.
    assert len({c.atomset for c in system.configurations}) == 1
    assert system.name == "A + B"


def test_build_mode_a_orientation_distance_names(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    # 3 orientations x 4 separations, default 'orientation/distance' naming.
    P = _P(spacing="explicit", separations="0.0, 1.0, 2.0, 3.0")
    P["number of orientations"] = 3
    system, stats = node._build(db, P, np.random.default_rng(12))

    names = [c.name for c in system.configurations]
    assert names[:4] == ["1,1", "1,2", "1,3", "1,4"]
    assert names[4:8] == ["2,1", "2,2", "2,3", "2,4"]
    assert names[-1] == "3,4"


def test_build_mode_a_separation_range(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    system, stats = node._build(db, P, np.random.default_rng(2))

    # Largest center-to-center separation reaches the maximum; smallest is a
    # slight overlap inside contact (but still positive).
    assert math.isclose(stats["max_separation"], 10.0, abs_tol=1.0e-6)
    assert 1.0 < stats["min_separation"] < 4.0


def test_build_mode_a_preserves_monomer_geometry(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(3))

    conf = system.configurations[0]
    xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
    # Atoms 0-2 are monomer A (water): O-H bonds ~0.9572 Å.
    oh = np.linalg.norm(xyz[1:3] - xyz[0], axis=1)
    assert np.allclose(oh, 0.9572, atol=1.0e-3)


def test_build_mode_a_monomer_a_is_fixed(db_two_waters):
    """With a single A conformer, monomer A is identical in every frame."""
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(7))

    reference = None
    for conf in system.configurations:
        xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
        A = xyz[:3]
        if reference is None:
            reference = A
            # A is centered at its center of mass at the origin.
            masses = np.asarray(conf.atoms.atomic_masses)[:3]
            com = (masses[:, None] * A).sum(axis=0) / masses.sum()
            assert np.allclose(com, [0.0, 0.0, 0.0], atol=1.0e-9)
        else:
            assert np.allclose(A, reference, atol=1.0e-9)


class _HarmonicEngine:
    """A fake engine whose energy is harmonic in the two fragments' separation.

    Duck-types the bits of seamm_mdi.MDIEngine that _energy_anchor uses, so the
    energy-contact machinery can be tested without MDI or a real code.
    """

    def __init__(self, nA, r0):
        self.nA = nA
        self.r0 = r0
        self._xyz = None

    def set_coordinates(self, xyz, units="bohr"):
        self._xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)

    def energy(self, units="hartree"):
        a = self._xyz[: self.nA].mean(axis=0)
        b = self._xyz[self.nA :].mean(axis=0)
        r = np.linalg.norm(b - a)
        return float((r - self.r0) ** 2)


def test_minimize_on_grid_parabola():
    node = dimer_builder_step.DimerBuilder()
    f = lambda d: (d - 2.7) ** 2 + 1.0  # noqa: E731
    d_min, k, n = node._minimize_on_grid(f, 1.0, 5.0, 9)
    assert d_min == pytest.approx(2.7, abs=0.05)
    assert 0 < k < n - 1  # interior minimum


def _assemble_along_z(A, B):
    def assemble(d):
        return np.vstack([A, B + np.array([0.0, 0.0, d])])

    return assemble


def test_energy_anchor_finds_minimum():
    node = dimer_builder_step.DimerBuilder()
    engine = _HarmonicEngine(nA=3, r0=3.2)
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    anchor = node._energy_anchor(engine, assemble, seed=2.5, P=_P())
    assert anchor == pytest.approx(3.2, abs=0.05)


def test_energy_anchor_falls_back_when_no_well():
    # Monotonically decreasing energy (no binding well) -> anchor at the seed.
    class _Repulsive:
        def set_coordinates(self, xyz, units="bohr"):
            self._xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)

        def energy(self, units="hartree"):
            a = self._xyz[:3].mean(axis=0)
            b = self._xyz[3:].mean(axis=0)
            return -float(np.linalg.norm(b - a))  # keeps falling as they separate

    node = dimer_builder_step.DimerBuilder()
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    anchor = node._energy_anchor(_Repulsive(), assemble, seed=2.9, P=_P())
    assert anchor == pytest.approx(2.9)


# --------------------------------------------------------------------------- #
# Energy-stratified radial sampling
# --------------------------------------------------------------------------- #


class _LJEngine:
    """A fake engine with a Lennard-Jones energy in the fragment separation.

    Duck-types the MDIEngine bits used by _energy_profile, giving a real well
    (depth eps at r = sigma·2^(1/6)) and a ΔE that decays to ~0 at large R.
    ``eps`` is in hartree, so the profile's kJ/mol well depth is eps·2625.5.
    """

    def __init__(self, nA, sigma, eps):
        self.nA = nA
        self.sigma = sigma
        self.eps = eps
        self._xyz = None

    def set_coordinates(self, xyz, units="bohr"):
        self._xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)

    def energy(self, units="hartree"):
        a = self._xyz[: self.nA].mean(axis=0)
        b = self._xyz[self.nA :].mean(axis=0)
        r = np.linalg.norm(b - a)
        sr6 = (self.sigma / r) ** 6
        return float(4.0 * self.eps * (sr6**2 - sr6))


def test_kBT_at_300K():
    node = dimer_builder_step.DimerBuilder()
    assert node._kBT(_P()) == pytest.approx(2.4943, abs=1.0e-3)


def test_energy_levels_parser():
    node = dimer_builder_step.DimerBuilder()
    levels = node._energy_levels(20.0, 2.5, _P())
    assert levels == pytest.approx([-20.0, -10.0, 0.0, 2.5, 12.5])


def test_energy_levels_rejects_bad_token():
    node = dimer_builder_step.DimerBuilder()
    with pytest.raises(ValueError):
        node._energy_levels(20.0, 2.5, _P(**{"energy levels": "De, __import__"}))


def test_energy_profile_finds_lj_well():
    node = dimer_builder_step.DimerBuilder()
    engine = _LJEngine(nA=3, sigma=3.0, eps=0.01)
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    ds, dE, d_min, De = node._energy_profile(engine, assemble, seed=3.0, P=_P())
    # The LJ minimum is at sigma * 2**(1/6) and the well depth is eps (in kJ/mol).
    assert d_min == pytest.approx(3.0 * 2.0 ** (1.0 / 6.0), abs=0.5)
    assert De == pytest.approx(0.01 * 2625.4996, rel=0.1)
    assert dE[-1] == pytest.approx(0.0, abs=1.0e-6)  # far-point reference


def test_interaction_energies_for_geometric_energy_scan():
    """A geometric scan with an energy engine still records per-point ΔE."""
    node = dimer_builder_step.DimerBuilder()
    engine = _LJEngine(nA=3, sigma=3.0, eps=0.01)
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    P = _P(**{"contact method": "energy", "spacing": "geometric"})
    distances, gap_ref, dE_at, De = node._plan_scan(engine, assemble, 3.2, P)
    assert dE_at is not None  # geometric+energy now carries ΔE
    assert De == pytest.approx(0.01 * 2625.4996, rel=0.15)  # LJ well depth (kJ/mol)
    # ΔE ~ 0 at the far point, negative in the well.
    assert dE_at(9.0) == pytest.approx(0.0, abs=0.5)
    assert min(dE_at(float(d)) for d in distances) < 0.0


def test_stratified_separations_hits_levels_twice():
    """A negative ΔE level is crossed once on the wall and once on the tail."""
    node = dimer_builder_step.DimerBuilder()
    ds = np.linspace(2.0, 10.0, 81)
    De = 20.0
    d_min = 3.2
    # A smooth well: -De at d_min, positive (wall) below, -> 0 above.
    dE = De * (((3.0 / ds) ** 12) - 2.0 * ((3.0 / ds) ** 6))
    dE = dE / (-dE.min()) * De  # normalize the depth to exactly De
    d_min = float(ds[np.argmin(dE)])
    distances = node._stratified_separations(ds, dE, d_min, De, _P())
    # The wall point, the well bottom, and both roots of -De/2 are present.
    assert distances[0] == pytest.approx(ds[0])
    assert any(abs(d - d_min) < 0.2 for d in distances)
    half = [d for d in distances if d < d_min], [d for d in distances if d > d_min]
    assert len(half[0]) >= 1 and len(half[1]) >= 1


def test_accept_orientation_reject_and_none():
    node = dimer_builder_step.DimerBuilder()
    rng = np.random.default_rng(0)
    P_reject = _P()  # reject shallow, min depth 1.0 kJ/mol
    assert node._accept_orientation(5.0, P_reject, rng) is True
    assert node._accept_orientation(0.2, P_reject, rng) is False
    P_none = _P(**{"orientation weighting": "none"})
    assert node._accept_orientation(0.0, P_none, rng) is True


def test_accept_orientation_downweight_probability():
    node = dimer_builder_step.DimerBuilder()
    P = _P(**{"orientation weighting": "downweight by depth"})
    rng = np.random.default_rng(42)
    # Deep well (10x the half-weight depth) -> kept nearly always.
    kept = sum(node._accept_orientation(10.0, P, rng) for _ in range(200))
    assert kept > 180


# --------------------------------------------------------------------------- #
# Sampling diagnostics (vendored dimer_analysis module, plotly)
# --------------------------------------------------------------------------- #


def test_dimer_analysis_metrics_and_summary():
    from dimer_builder_step import dimer_analysis

    rng = np.random.default_rng(0)
    ensemble = []
    for i in range(20):
        z = 3.0 + 0.1 * i
        ensemble.append(
            dimer_analysis.Dimer(
                symbols_A=["O", "H", "H"],
                xyz_A=np.array(
                    [[0.0, 0.0, 0.0], [0.76, 0.59, 0.0], [-0.76, 0.59, 0.0]]
                ),
                symbols_B=["O", "H", "H"],
                xyz_B=np.array([[0.0, 0.0, z], [0.76, 0.59, z], [-0.76, 0.59, z]]),
                energy=float(rng.normal()),
                separation=z,
                orientation=1,
            )
        )
    metrics = dimer_analysis.compute_metrics(ensemble)
    s = dimer_analysis.summarize(metrics)
    assert s["n"] == 20
    assert metrics.has_energy
    assert "energy_flatness" in s


def test_build_collects_ensemble_when_requested(db_two_waters):
    node = dimer_builder_step.DimerBuilder()
    P = _P(**{"analysis plots": "basic"})
    _, stats = node._build(db_two_waters, P, np.random.default_rng(1))

    ensemble = stats["ensemble"]
    assert len(ensemble) == stats["n_configurations"]
    d = ensemble[0]
    assert d.symbols_A == ["O", "H", "H"] and d.symbols_B == ["O", "H", "H"]
    assert d.xyz_A.shape == (3, 3) and d.xyz_B.shape == (3, 3)
    # van der Waals contact method -> no interaction energies collected.
    assert d.energy is None


def test_build_skips_ensemble_when_none(db_two_waters):
    node = dimer_builder_step.DimerBuilder()
    _, stats = node._build(db_two_waters, _P(), np.random.default_rng(1))
    assert stats["ensemble"] == []


def test_make_dashboard_returns_plotly_figure(db_two_waters):
    import json

    from dimer_builder_step import dimer_analysis

    node = dimer_builder_step.DimerBuilder()
    P = _P(**{"analysis plots": "basic"})
    _, stats = node._build(db_two_waters, P, np.random.default_rng(2))
    metrics = dimer_analysis.compute_metrics(stats["ensemble"])
    figure = dimer_analysis.make_dashboard(metrics, title="test")
    # A native plotly go.Figure that serializes to the SEAMM .graph format.
    assert len(figure.data) > 0
    payload = json.loads(figure.to_json())
    assert "data" in payload and "layout" in payload


def test_make_panels_individual_figures(db_two_waters):
    from dimer_builder_step import dimer_analysis

    node = dimer_builder_step.DimerBuilder()
    _, stats = node._build(
        db_two_waters, _P(**{"analysis plots": "detailed"}), np.random.default_rng(2)
    )
    metrics = dimer_analysis.compute_metrics(stats["ensemble"])
    panels = dimer_analysis.make_panels(metrics)
    # vdW contact -> no energies -> the four geometry panels, no energy panels.
    assert set(panels) == {"separation", "contact", "approach", "orientation"}
    assert all(len(fig.data) > 0 for fig in panels.values())


def test_detailed_writes_panel_graphs(db_two_waters, tmp_path):
    from unittest import mock

    node = dimer_builder_step.DimerBuilder()
    _, stats = node._build(
        db_two_waters, _P(**{"analysis plots": "detailed"}), np.random.default_rng(2)
    )
    # 'directory' is a read-only property (flowchart root + node id); patch it.
    with mock.patch.object(
        type(node),
        "directory",
        new_callable=mock.PropertyMock,
        return_value=str(tmp_path),
    ):
        node._run_diagnostics(stats["ensemble"], "detailed", stats["system"])
    written = {p.name for p in tmp_path.glob("dimer_sampling*.graph")}
    assert "dimer_sampling.graph" in written  # the combined dashboard
    assert "dimer_sampling_separation.graph" in written  # a per-panel graph
    assert "dimer_sampling_contact.graph" in written


def test_direction_angles_known():
    node = dimer_builder_step.DimerBuilder()
    assert node._direction_angles([0.0, 0.0, 1.0]) == (0.0, 0.0)
    th, ph = node._direction_angles([1.0, 0.0, 0.0])
    assert math.isclose(th, 90.0) and math.isclose(ph, 0.0)
    th, ph = node._direction_angles([0.0, 1.0, 0.0])
    assert math.isclose(th, 90.0) and math.isclose(ph, 90.0)


def test_euler_zyz_identity_and_roundtrip():
    node = dimer_builder_step.DimerBuilder()
    # Identity -> all zero.
    assert node._euler_zyz(np.eye(3)) == (0.0, 0.0, 0.0)

    # Round-trip: build R = Rz(a) Ry(b) Rz(c), extract, rebuild, compare.
    def Rz(t):
        c, s = math.cos(t), math.sin(t)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

    def Ry(t):
        c, s = math.cos(t), math.sin(t)
        return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])

    a, b, c = math.radians(40), math.radians(70), math.radians(20)
    R = Rz(a) @ Ry(b) @ Rz(c)
    alpha, beta, gamma = node._euler_zyz(R)
    R2 = Rz(math.radians(alpha)) @ Ry(math.radians(beta)) @ Rz(math.radians(gamma))
    assert np.allclose(R, R2, atol=1.0e-9)


def test_build_mode_a_tags_geometry_properties(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(11))

    conf = system.configurations[0]
    for name in (
        "dimer separation",
        "dimer gap",
        "dimer orientation",
        "approach theta",
        "approach phi",
        "movable alpha",
        "movable beta",
        "movable gamma",
    ):
        assert conf.properties.exists(f"{name}#DimerBuilder#scan")


def test_orient_to_principal_axes_diagonalizes_inertia():
    """After reorientation the inertia tensor is diagonal (axes on x/y/z)."""
    node = dimer_builder_step.DimerBuilder()
    # A tilted water, arbitrary placement.
    r0, theta0 = 0.9572, 104.52
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))
    xyz = np.array([[0.0, 0.0, 0.0], [x, 0.0, z], [-x, 0.0, z]]) + 5.0
    R = node._orient_to_principal_axes(xyz, [15.999, 1.008, 1.008])

    masses = np.array([15.999, 1.008, 1.008])
    com = (masses[:, None] * R).sum(axis=0) / masses.sum()
    inertia = np.zeros((3, 3))
    for m, r in zip(masses, R - com):
        inertia += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    off_diag = inertia - np.diag(np.diag(inertia))
    assert np.allclose(off_diag, 0.0, atol=1.0e-9)
    # COM at the origin, and a proper rotation preserves bond lengths.
    assert np.allclose(com, 0.0, atol=1.0e-9)
    assert np.allclose(np.linalg.norm(R[1:] - R[0], axis=1), 0.9572, atol=1.0e-6)


def test_build_mode_a_tags_properties(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(4))

    conf = system.configurations[0]
    assert conf.properties.exists("dimer separation#DimerBuilder#scan")
    assert conf.properties.exists("dimer gap#DimerBuilder#scan")
    assert conf.properties.exists("dimer orientation#DimerBuilder#scan")


def test_build_mode_a_no_severe_overlap(db_two_waters):
    """No atom pair across the two monomers is grossly inside vdW contact."""
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(5))

    rA = vdw_radii(["O", "H", "H"])
    for conf in system.configurations:
        xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
        A = xyz[:3]
        B = xyz[3:]
        D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
        Rsum = rA[:, None] + rA[None, :]
        # Innermost point allows a slight overlap; require no pair closer than
        # 70% of the vdW sum.
        assert np.all(D > 0.7 * Rsum)


def test_build_mode_a_creates_fixed_movable_subsets(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(8))

    conf = system.configurations[0]
    atom_ids = conf.atoms.ids
    # A is 'fixed' (first 3 atoms), B is 'movable' (last 3).
    assert set(_subset_atom_ids(db, conf, "fixed")) == set(atom_ids[:3])
    assert set(_subset_atom_ids(db, conf, "movable")) == set(atom_ids[3:])
    # Every configuration carries both subsets.
    for c in system.configurations:
        assert len(_subset_atom_ids(db, c, "fixed")) == 3
        assert len(_subset_atom_ids(db, c, "movable")) == 3


# --------------------------------------------------------------------------- #
# Full build — Mode B (prepared dimers)
# --------------------------------------------------------------------------- #


def test_build_mode_b_from_prepared_dimer(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    # A prepared water dimer, B displaced 3 Å along z.
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"

    system, stats = node._build(db, P, np.random.default_rng(6))

    assert stats["n_configurations"] == 8  # 1 dimer x 8 separations
    assert len(system.configurations) == 8
    assert all(c.n_atoms == 6 for c in system.configurations)
    assert len({c.atomset for c in system.configurations}) == 1


def test_build_mode_b_default_last_molecule_movable(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"
    system, stats = node._build(db, P, np.random.default_rng(9))

    conf = system.configurations[0]
    out_ids = conf.atoms.ids
    # No user subsets -> last molecule (atoms 3-5) is movable, rest fixed.
    assert set(_subset_atom_ids(db, conf, "movable")) == set(out_ids[3:])
    assert set(_subset_atom_ids(db, conf, "fixed")) == set(out_ids[:3])


def test_build_mode_b_honors_user_subsets(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    dimer_sys = db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")
    d0 = dimer_sys.configuration

    # The user designates the FIRST molecule as movable (overriding the default).
    for nm in ("fixed", "movable"):
        if not db.templates.exists(nm, "general"):
            db.templates.create(name=nm, category="general")
    ids = d0.atoms.ids
    d0.subsets.create(template=db.templates.get("movable", "general"), atoms=ids[:3])
    d0.subsets.create(template=db.templates.get("fixed", "general"), atoms=ids[3:])

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"
    system, stats = node._build(db, P, np.random.default_rng(10))

    conf = system.configurations[0]
    out_ids = conf.atoms.ids
    # Movable follows the user's choice: the first molecule.
    assert set(_subset_atom_ids(db, conf, "movable")) == set(out_ids[:3])
