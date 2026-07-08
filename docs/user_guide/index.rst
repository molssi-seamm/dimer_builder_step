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
  thin out in the tail), ``linear``, or ``explicit`` (a list you provide in
  **Separations**, as gaps beyond contact).
* **Number of separations** -- how many points per scan (for geometric/linear).

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

What is stored
==============

All the generated configurations become conformers of a **new system** (named
from the two monomers by default). On each configuration the step:

* records the scan geometry as properties -- the separation and the gap beyond
  contact, the approach-direction angles, and the movable group's orientation
  (Euler) angles -- so configurations can be filtered and analyzed later
  (toggle with **Save scan variables as properties**); and
* marks the two pieces as ``fixed`` and ``movable`` **subsets**, so a downstream
  interaction-energy or counterpoise calculation can find them.

**Name the configurations** controls the configuration names:
``orientation,distance`` (e.g. ``2,1`` -- orientation 2, point 1; the points of
one scan group together), ``separation``, or ``sequential``.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
