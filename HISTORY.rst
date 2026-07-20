=======
History
=======

2026.7.20 -- Energy-stratified sampling and sampling diagnostics
    * Added an **energy-stratified** spacing for the radial scan (with the energy
      contact method): instead of a fixed geometric ladder, the scan points are
      placed to hit a spread of interaction-energy levels (e.g. the well depth,
      half the well depth, zero, and a few multiples of kBT) along the energy
      profile. This gives a sample that is flat in energy — from the repulsive
      wall through the attractive well to the long-range tail — which is what a
      machine-learned force field needs, rather than piling most points in the
      shallow near-zero-interaction region.
    * Orientations can be weighted by how deeply they bind: shallow, glancing
      orientations can be rejected, or kept with a probability that grows with
      the well depth, so effort concentrates on the physically important
      geometries.
    * The energy-stratified scan is anchored on the actual energy minimum found
      from the engine, with the repulsive wall sampled by energy, so it no longer
      relies on the van der Waals estimate being accurate (the "innermost gap"
      setting therefore does not apply to this spacing).
    * The interaction energy is now recorded for every configuration whenever the
      energy contact method is used (any spacing), and saved as a property.
    * Added **sampling diagnostics**: after a build, the step reports a short
      summary and, with "Sampling diagnostics" set to ``basic`` or ``detailed``,
      writes interactive graphs for the Dashboard — separation coverage, contact
      distances, approach direction, relative orientation, and (when energies
      were computed) the interaction-energy distribution and binding-curve
      envelope. ``detailed`` also writes each panel as its own graph. Extra image
      formats (PDF, PNG, SVG, ...) can be requested with ``graph-formats`` in
      ``seamm.ini``.

2026.7.9 -- ORCA engine for energy contact, and consistent CoM separation
    * The energy-based contact search can now use **ORCA** (via a Model
      Chemistry step), in addition to MOPAC and xTB. Any MDI-capable model
      chemistry works; the model chemistry's method and, for ORCA, its basis set
      are passed to the engine automatically.
    * ``seamm_mdi`` is now a declared dependency (it was previously imported
      lazily), so the energy contact method works out of the box.
    * The center-to-center separation for prepared dimers is now measured between
      the centers of mass of the "fixed" and "movable" groups, matching the
      two-monomer-sets path (which was already mass-weighted). Previously the
      prepared-dimers path used the unweighted geometric center, so the reported
      ``dimer separation`` / ``dimer gap`` properties now have a single,
      consistent meaning across both input modes.

2026.7.6 -- Energy-based contact and GUI refinements
    * Added an ``energy`` contact method: with a Model Chemistry step before the
      Dimer Builder step, the contact distance is found from the energy minimum
      along each approach direction (falling back to the van der Waals estimate
      for orientations with no binding well). It uses the new ``seamm_mdi`` MDI
      driver to evaluate the energy, and reports which model chemistry was used
      and how many times it was called.
    * Configuration names now use a comma (e.g. ``2,1``) rather than a slash,
      which is reserved for separating system and configuration names.
    * The dialog no longer shows the unused Flowchart tab, and reminds you to add
      a Model Chemistry step when the energy contact method is chosen and none
      precedes this step.

2026.6.30 -- Initial release
    * Generates sets of dimer (molecule-pair) configurations across a range of
      separations and relative orientations, for building interaction-energy data
      sets and training sets for machine-learned force fields.
    * Two ways to provide the input structures:

        - two sets of monomer conformers, which are assembled into dimers at
          random relative orientations; or
        - prepared complexes, each scanned along the axis between a "fixed" and a
          "movable" group (taken from subsets if present, otherwise the last
          molecule is movable and the rest are fixed).
    * For each orientation, scans the center-to-center separation from just inside
      the van der Waals contact distance out to a chosen maximum, with geometric
      (default), linear, or explicit spacing.
    * Records the scan geometry on every configuration as properties -- separation,
      gap beyond contact, the approach-direction angles, and the movable group's
      orientation angles -- and marks the two molecules as "fixed" and "movable"
      subsets, so the pieces are easy to find downstream.
    * Stores all the generated structures as conformers of a new system, named by
      orientation and point (1/1, 1/2, ...) by default.
