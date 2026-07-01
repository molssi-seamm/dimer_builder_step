=======
History
=======

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
