.. _user-guide:

**********
User Guide
**********

The Dimer Builder step generates sets of configurations of a pair (or small
cluster) of molecules across a range of separations and relative orientations.
The configurations are stored as conformers of a new system, ready to be fed to
a quantum-chemistry step for computing interaction energies, or exported as a
training set for machine-learned force fields (in the spirit of the DES370K
dataset).

Input
=====

Set **Input** to choose how the structures are provided:

**two monomer sets**
    Assemble dimers from conformers of *monomer A* and *monomer B* at random
    relative orientations. Each of the two monomers is given as:

    * ``current`` -- the current system;
    * a **system name** -- use that system's configurations as the conformer
      pool, selected by the *using configurations* control (``all``, ``last``,
      ``first``, or a name pattern); or
    * a ``$variable`` holding a list of configurations (all are used).

**prepared dimers**
    Take already-assembled complexes (from *monomer A*) and scan each one,
    sliding a "movable" group relative to a "fixed" group. The two groups come
    from ``fixed`` / ``movable`` subsets on the input if present (so you can
    group several molecules per side); otherwise the last molecule is movable
    and the rest are fixed.

Orientation sampling
====================

For "two monomer sets", **Number of orientations** samples are generated. Each
sample independently draws a random conformer of A, a random conformer of B, and
a random relative orientation. Monomer A is held fixed (centered on its
principal axes) so it is a stable visual anchor across the whole run; monomer B
carries the randomness. Set **Random seed** to an integer for reproducible
orientations, or leave it as ``random``.

The radial scan
===============

For each orientation the center-to-center separation is scanned from just inside
contact out to a maximum:

* **Find contact using** -- ``van der Waals radii`` (a fast geometric estimate)
  or ``energy`` (see below).
* **Innermost gap** -- how far inside the contact distance to start (a small
  negative value gives a slight overlap).
* **Maximum separation** -- the largest center-to-center distance.
* **Spacing** -- ``geometric`` (the default; points cluster near contact and
  thin out in the tail), ``linear``, ``explicit`` (a list you provide in
  **Separations**, as gaps beyond contact), or ``energy-stratified`` (see
  below; requires the ``energy`` contact method).
* **Number of separations** -- how many points per scan (for geometric/linear;
  for energy-stratified it sets the resolution of the ΔE profile).

Energy-based contact
--------------------

With **Find contact using** = ``energy``, the step drives a quantum-chemistry
engine to find the energy **minimum** along each approach direction and anchors
the scan there; orientations with no binding well fall back to the van der Waals
contact. This requires a **Model Chemistry step before the Dimer Builder step**
to define the engine and method; the dialog reminds you if one is missing. The
engine is driven over `MDI <https://molssi-mdi.github.io/MDI_Library/>`_, so any
MDI-capable model chemistry works -- **MOPAC** or **xTB** (semiempirical), or
**ORCA** (HF, MP2, or an analytic-gradient DFT functional; ORCA methods without
an analytic gradient, such as DLPNO-CCSD(T), are not offered for this). A cheap
method is usually the right choice: the energy only *places* the scan, it does
not produce the final data set, and ORCA in particular runs its binary once per
geometry (reusing orbitals between points) so it is heavier than the
semiempirical engines. The step reports which model chemistry was used and how
many times it was called.

Energy-stratified sampling
--------------------------

With **Spacing** = ``energy-stratified`` (which needs the ``energy`` contact
method) the step builds a set of dimers that is **flat in interaction energy**
rather than evenly spaced in distance. This is usually what you want for a
machine-learned force field: the model must be accurate everywhere it will be
queried -- the repulsive wall, the bottom of the attractive well, and the
long-range tail -- so the training set should cover that whole energy range
evenly, instead of piling most configurations in the shallow, nearly
non-interacting region that uniform sampling produces.

It works by pooling candidate configurations from *all* orientations and then
down-selecting about **Target configurations** of them. **Down-select by**
chooses how:

* ``energy bins + diversity`` (default) -- sort the candidates into interaction-
  energy bins (flat in energy) and, within each bin, keep a *geometrically
  diverse, de-duplicated* subset by clustering the collective variables
  (separation, approach direction, relative orientation, closest contact). This
  gives a set that is flat in energy, reaches deep into the attractive well, and
  is not dominated by near-identical geometries.
* ``descriptor diversity`` -- a single global clustering (the DIRECT method) over
  those collective variables **plus** the interaction energy, keeping one per
  cluster. **Energy weight** sets how strongly ΔE counts relative to each
  geometric variable (larger = flatter in energy, less geometric spread). This
  maximizes geometric diversity but the energy flatness depends on the weight.
* ``energy bins`` -- energy bins with a plain random pick per bin (flat in
  energy, no geometric de-duplication).

Other controls:

* **Number of energy bins** -- how many ΔE bins to spread the kept
  configurations across (the two binned methods).
* **Target configurations** -- the approximate total to keep. Deeply bound
  geometries are rare, so the total may come out smaller; raise **Number of
  orientations** to find more deep configurations.
* **ΔE levels** -- sets the energy window. The most repulsive value (default
  ``+5*kBT``) caps the wall, so no configuration is pushed to an absurd
  repulsive energy; the symbols ``De`` (well depth) and ``kBT`` (thermal energy
  at the **Sampling temperature**) may be used.
* **Weight orientations by well depth** -- an optional pre-filter. The default
  ``none`` keeps every orientation; ``reject shallow orientations`` /
  ``downweight by depth`` bias toward the more strongly bound orientations first
  (using **Minimum well depth**).

Because the interaction energy varies almost entirely at short range, a
flat-in-energy set naturally has most of its configurations at short
separations; that is expected.

What is stored
==============

All the generated configurations become conformers of a **new system** (named
from the two monomers by default). On each configuration the step:

* records the scan geometry as properties -- the separation and the gap beyond
  contact, the approach-direction angles, and the movable group's orientation
  (Euler) angles, plus the interaction energy and well depth when the ``energy``
  contact method was used -- so configurations can be filtered and analyzed
  later (toggle with **Save scan variables as properties**); and
* marks the two pieces as ``fixed`` and ``movable`` **subsets**, so a downstream
  interaction-energy or counterpoise calculation can find them.

**Name the configurations** controls the configuration names:
``orientation,distance`` (e.g. ``2,1`` -- orientation 2, point 1; the points of
one scan group together), ``separation``, or ``sequential``.

Sampling diagnostics
====================

Set **Sampling diagnostics** to check *what kind of coverage a run actually
produced* -- how the configurations are spread over separation, contact
distance, approach direction, relative orientation, and (with the ``energy``
contact method) interaction energy:

* ``none`` -- skip the diagnostics.
* ``basic`` -- print a short numerical summary and write a single combined
  dashboard, ``dimer_sampling.graph``, that opens interactively in the SEAMM
  Dashboard.
* ``detailed`` -- also write each panel as its own graph
  (``dimer_sampling_<panel>.graph``) for a closer look.

To also save the graphs as static images, add ``graph-formats`` to the Dimer
Builder section of ``~/SEAMM/dimer_builder.ini`` (or pass ``--graph-formats``),
e.g. ``graph-formats = png pdf svg``.

Reading the dashboard
---------------------

The header line summarizes the run: **N** (number of configurations), the
separation range **R**, the **approach concentration**, and -- when energies
were computed -- the minimum interaction energy **ΔE_min** and the percentage of
configurations that are attractive (ΔE < 0). The panels are:

**Separation coverage**
    Histogram of the center-to-center separation *R*. Shows how far apart the
    molecules are across the set. With energy-stratified sampling most points
    sit at short *R* (that is where the interaction changes), with a thin tail
    to the maximum separation.

**Contact distribution**
    The shaded histogram is the *closest* atom-atom distance between the two
    molecules in each configuration -- the practical measure of how tightly they
    touch. The overlaid lines are the most common element-pair distances (e.g.
    O--H, H--H); a peak at short O--H, for example, is the fingerprint of
    hydrogen bonding. Watch for a build-up at very short distances, which means
    hard atomic clashes.

**Approach direction (A frame)**
    Where molecule B sits relative to molecule A, drawn on a flattened sphere
    (a Mollweide map) in A's own frame; with energies, each point is colored by
    its ΔE. Points spread evenly over the whole map mean B approaches from all
    directions equally (isotropic); points clustered in one region mean the
    sampling favors particular docking directions. The single-number version of
    this is the **approach concentration**: **0 = perfectly isotropic**, **1 =
    all from one direction**. A low value with a directional interaction like a
    hydrogen bond signals that the orientations are *not* concentrated on the
    physically important motif.

**Relative orientation**
    Histogram of the angle between the two molecules' long axes -- how the
    partners are turned relative to each other, independent of which side they
    approach from.

**Energy distribution** *(energy contact method only)*
    Histogram of the interaction energy ΔE. For energy-stratified sampling this
    should be roughly **flat** from the repulsive side (positive) through zero to
    the bottom of the well (negative). The title reports a **flatness CV**
    (coefficient of variation of the bar heights): **0 is perfectly flat**;
    values well below ~0.5 are good, while a value near ~2 means the sample is
    piled up in one energy range (what uniform sampling gives).

**Binding-curve envelope** *(energy contact method only)*
    A 2-D density of ΔE versus separation *R* -- the interaction energy as the
    molecules come together, pooled over all orientations. The lower edge traces
    the best (deepest) binding curve; the spread at each *R* shows how much the
    energy depends on orientation.

A "good" energy-stratified run therefore shows a flat **Energy distribution**
(low flatness CV) that reaches well into the attractive region, a clean
**Binding-curve envelope**, and no pile-up of hard clashes in the **Contact
distribution**. The **approach concentration** tells you separately whether the
*orientational* coverage is biased or isotropic.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
