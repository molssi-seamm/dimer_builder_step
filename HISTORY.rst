=======
History
=======

2026.7.20 -- Energy-stratified sampling and sampling diagnostics
    * Added an **energy-stratified** spacing (with the energy contact method) that
      produces a set of dimers **flat in interaction energy** — evenly covering
      the repulsive wall, the attractive well, and the long-range tail — which is
      what a machine-learned force field needs, rather than the pile-up near zero
      interaction that uniform sampling gives. It works globally: candidate
      configurations are pooled across all orientations, sorted into
      interaction-energy bins, and each bin is capped at the same count. Two
      controls, "Number of energy bins" and "Target configurations", set the
      resolution and the approximate total (deep-well energies are rare, so those
      bins — and the total — may come out smaller; raise the number of
      orientations to fill them).
    * The scan is anchored on the actual energy minimum found from the engine, and
      the repulsive side is capped at the largest "ΔE level" (default +5·kBT), so
      configurations are never pushed to absurd repulsive energies and the
      van der Waals estimate is only a starting guess (the "innermost gap" setting
      does not apply to this spacing).
    * Orientations may optionally be pre-filtered by how deeply they bind (reject
      shallow, or downweight by depth); the default keeps every orientation and
      lets the global stratification balance the set.
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
