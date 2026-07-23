# -*- coding: utf-8 -*-
"""
Control parameters for the Dimer Builder step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class DimerBuilderParameters(seamm.Parameters):
    """
    The control parameters for Dimer Builder.

    The keys are the parameters for this plug-in; each value is a dictionary
    describing that parameter (default, kind, units, enumeration, format,
    description, help). This is the Tier-1 set (engine-free generation); the
    "energy" contact method and any further options are layered on later.

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
    """

    parameters = {
        # ------------------------------------------------------------------ #
        # Input: where the monomers (or prepared dimers) come from
        # ------------------------------------------------------------------ #
        "input mode": {
            "default": "two monomer sets",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("two monomer sets", "prepared dimers"),
            "format_string": "",
            "description": "Input:",
            "help_text": (
                "'two monomer sets' assembles dimers from conformers of monomer "
                "A and monomer B at random relative orientations. 'prepared "
                "dimers' takes already-assembled dimers (from 'monomer A') and "
                "scans each, splitting it into its two molecules by connectivity."
            ),
        },
        "monomer A": {
            "default": "current",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current",),
            "format_string": "",
            "description": "Monomer A:",
            "help_text": (
                "The source of monomer A: 'current' for the current system, a "
                "system name, or a variable ($name) holding a list of "
                "configurations. In 'prepared dimers' mode this is the source of "
                "the dimers."
            ),
        },
        "monomer A configurations": {
            "default": "all",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "all",
                "last",
                "first",
                "name is",
                "name matches",
                "name regexp",
            ),
            "format_string": "",
            "description": "using configurations:",
            "help_text": (
                "Which configurations of the monomer A system to use as the "
                "conformer pool. Ignored when monomer A is a variable holding a "
                "list of configurations (all of them are used)."
            ),
        },
        "monomer A configuration name": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "matching:",
            "help_text": (
                "The configuration name or pattern, used with 'name is', 'name "
                "matches', or 'name regexp'."
            ),
        },
        "monomer B": {
            "default": "current",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current",),
            "format_string": "",
            "description": "Monomer B:",
            "help_text": (
                "The source of monomer B: 'current' for the current system, a "
                "system name, or a variable ($name) holding a list of "
                "configurations. Ignored in 'prepared dimers' mode."
            ),
        },
        "monomer B configurations": {
            "default": "all",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "all",
                "last",
                "first",
                "name is",
                "name matches",
                "name regexp",
            ),
            "format_string": "",
            "description": "using configurations:",
            "help_text": (
                "Which configurations of the monomer B system to use as the "
                "conformer pool. Ignored when monomer B is a variable holding a "
                "list of configurations (all of them are used)."
            ),
        },
        "monomer B configuration name": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "matching:",
            "help_text": (
                "The configuration name or pattern, used with 'name is', 'name "
                "matches', or 'name regexp'."
            ),
        },
        # ------------------------------------------------------------------ #
        # Orientation sampling (N random pairings)
        # ------------------------------------------------------------------ #
        "number of orientations": {
            "default": 10,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of orientations:",
            "help_text": (
                "How many samples to generate. Each sample draws a random "
                "conformer of monomer A, a random conformer of monomer B, and a "
                "random relative orientation, then is scanned radially."
            ),
        },
        "random seed": {
            "default": "random",
            "kind": "string",
            "default_units": "",
            "enumeration": ("random",),
            "format_string": "",
            "description": "Random seed:",
            "help_text": (
                "The seed for the random-number generator. Use 'random' for a "
                "fresh, non-reproducible seed, or an integer for reproducible "
                "orientations."
            ),
        },
        # ------------------------------------------------------------------ #
        # Radial scan: surface-gap coordinate, geometric spacing
        # ------------------------------------------------------------------ #
        "contact method": {
            "default": "van der Waals radii",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("van der Waals radii", "energy"),
            "format_string": "",
            "description": "Find contact using:",
            "help_text": (
                "How to locate the contact distance that anchors the scan. 'van "
                "der Waals radii' is a fast geometric estimate (Tier 1). 'energy' "
                "uses the energy from the sub-flowchart to find the onset of "
                "repulsion (requires an engine in the sub-flowchart)."
            ),
        },
        "innermost gap": {
            "default": -0.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Innermost gap:",
            "help_text": (
                "The closest point of the scan, measured as the gap beyond "
                "contact (0 = touching). A negative value starts inside contact, "
                "i.e. with a slight overlap. Not used for 'energy-stratified' "
                "spacing, which anchors on the energy minimum and sets its inner "
                "bound by energy (the repulsive-wall target level)."
            ),
        },
        "maximum separation": {
            "default": 10.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Maximum separation:",
            "help_text": ("The largest center-of-mass separation to scan out to."),
        },
        "spacing": {
            "default": "geometric",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("geometric", "linear", "explicit", "energy-stratified"),
            "format_string": "",
            "description": "Spacing:",
            "help_text": (
                "How the scan points are distributed along the gap coordinate. "
                "'geometric' clusters points near contact and thins them in the "
                "tail; 'linear' spaces them evenly; 'explicit' uses the list in "
                "'separations'. 'energy-stratified' places points at target "
                "interaction-energy levels along the ΔE(R) profile (requires the "
                "'energy' contact method), giving a flat-in-energy sample from the "
                "repulsive wall through the well to the asymptote."
            ),
        },
        "number of separations": {
            "default": 15,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of separations:",
            "help_text": (
                "How many points to place along each radial profile, for "
                "'geometric' or 'linear' spacing. For 'energy-stratified' spacing "
                "this is the resolution of the ΔE(R) profile (the number of "
                "energy evaluations), not the number of output points -- those are "
                "set by the 'energy levels'."
            ),
        },
        # ------------------------------------------------------------------ #
        # Energy-stratified spacing (needs the 'energy' contact method)
        # ------------------------------------------------------------------ #
        "energy levels": {
            "default": "-De, -De/2, 0, kBT, 5*kBT",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "ΔE levels:",
            "help_text": (
                "Sets the interaction-energy window for 'energy-stratified' "
                "spacing, as a comma-separated list of linear expressions in the "
                "well depth 'De' and thermal energy 'kBT' (kJ/mol) -- e.g. "
                "'-De, -De/2, 0, kBT, 5*kBT'. The largest positive value caps the "
                "repulsive wall (points above it are dropped); the candidates are "
                "then globally binned across this range and capped per bin."
            ),
        },
        "sampling temperature": {
            "default": 300.0,
            "kind": "float",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Sampling temperature:",
            "help_text": (
                "The temperature that sets kBT for the thermal 'energy levels' in "
                "'energy-stratified' spacing."
            ),
        },
        "orientation weighting": {
            "default": "none",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "none",
                "reject shallow orientations",
                "downweight by depth",
            ),
            "format_string": "",
            "description": "Weight orientations by well depth:",
            "help_text": (
                "Optional per-orientation pre-filter for 'energy-stratified' "
                "spacing ('two monomer sets' mode). The default 'none' keeps every "
                "orientation and lets the global energy-stratification balance the "
                "ensemble. 'reject shallow orientations' skips any orientation whose "
                "well is shallower than the 'minimum well depth'; 'downweight by "
                "depth' keeps each with probability 1 - 2^(-De/minimum well depth). "
                "Prepared dimers are always kept."
            ),
        },
        "minimum well depth": {
            "default": 2.5,
            "kind": "float",
            "default_units": "kJ/mol",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Minimum well depth:",
            "help_text": (
                "For 'reject shallow orientations', the smallest ΔE well depth an "
                "orientation must have to be kept (defaults to ~kBT at 300 K). For "
                "'downweight by depth', the half-weight depth. Ignored when "
                "orientation weighting is 'none'."
            ),
        },
        "selection method": {
            "default": "energy bins + diversity",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "energy bins + diversity",
                "descriptor diversity",
                "energy bins",
            ),
            "format_string": "",
            "description": "Down-select by:",
            "help_text": (
                "How to down-select the pooled candidate configurations for "
                "'energy-stratified' spacing (all target about 'target "
                "configurations'). 'energy bins + diversity' bins by interaction "
                "energy (flat in energy) and, within each bin, keeps a "
                "geometrically diverse, de-duplicated subset (DIRECT clustering of "
                "the geometric collective variables). 'descriptor diversity' is a "
                "single global DIRECT clustering over those variables PLUS the "
                "interaction energy (scaled by 'energy weight') -- maximally "
                "diverse, but flatness depends on the weight. 'energy bins' caps "
                "each energy bin equally with a random pick (flat in energy, no "
                "geometric de-duplication)."
            ),
        },
        "energy weight": {
            "default": 8.0,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Energy weight:",
            "help_text": (
                "For the 'descriptor diversity' selection: how strongly the "
                "interaction energy counts relative to each geometric collective "
                "variable in the clustering (1 = equal; larger keeps the sample "
                "flatter in energy at the cost of some geometric diversity)."
            ),
        },
        "number of energy bins": {
            "default": 12,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of energy bins:",
            "help_text": (
                "For 'energy-stratified' spacing: the number of interaction-energy "
                "bins the pooled candidate configurations are sorted into for the "
                "global stratification. Each bin is capped at the same count, so the "
                "kept ensemble is flat in energy from the repulsive wall through the "
                "well to the asymptote."
            ),
        },
        "target configurations": {
            "default": 300,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Target configurations:",
            "help_text": (
                "For 'energy-stratified' spacing: the approximate total number of "
                "configurations to keep. The per-bin cap is this divided by the "
                "number of energy bins; sparsely-populated energy ranges (e.g. the "
                "deep well) may yield fewer, so the actual total can be smaller. "
                "Increase the number of orientations to fill the deep bins."
            ),
        },
        "tail coverage": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Add long-range distance coverage:",
            "help_text": (
                "Energy-stratified selection is flat in energy, which starves the "
                "weak long-range tail (5+ Å maps to ≈0 ΔE, so it gets only ~one "
                "bin's worth) even though molecular dynamics traverses it "
                "constantly. With this on, the kept set is supplemented with a "
                "distance-coverage floor and a few far-separation anchors so the "
                "tail is sampled independently of interaction strength: energy-flat "
                "where the interaction is strong, distance-dense in the weak tail."
            ),
        },
        "tail minimum separation": {
            "default": 4.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Tail coverage from:",
            "help_text": (
                "The start of the distance-coverage floor. From here out to the "
                "maximum separation the kept set is guaranteed a minimum number of "
                "configurations in each separation bin (over a spread of "
                "orientations), on top of the energy-based selection."
            ),
        },
        "tail spacing": {
            "default": 0.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Tail coverage spacing:",
            "help_text": (
                "The separation-bin width of the distance-coverage floor "
                "(smaller = denser tail coverage)."
            ),
        },
        "tail configurations per bin": {
            "default": 2,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Tail configurations per bin:",
            "help_text": (
                "The minimum number of configurations to keep in each separation "
                "bin of the distance-coverage floor (chosen for geometric diversity "
                "when the energy selection did not already supply them)."
            ),
        },
        "asymptote anchors": {
            "default": 20,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Asymptote anchors:",
            "help_text": (
                "How many far-separation configurations (between the maximum "
                "separation and the 'anchor separation', over a spread of "
                "orientations) to add to pin the interaction energy to zero at "
                "large distance. 0 disables them."
            ),
        },
        "anchor separation": {
            "default": 15.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Anchor separation:",
            "help_text": (
                "The largest center-to-center separation sampled for the asymptote "
                "anchors (the profile is extended to here so the ≈0 tail is "
                "represented out to this distance)."
            ),
        },
        "separations": {
            "default": "",
            "kind": "string",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Separations:",
            "help_text": (
                "An explicit list of gaps beyond contact (e.g. "
                "'-0.5, 0.0, 0.5, 1:5:0.5'), used when spacing is 'explicit'."
            ),
        },
        # ------------------------------------------------------------------ #
        # Output: a new system of dimer conformers
        # ------------------------------------------------------------------ #
        "system name": {
            "default": "from monomers",
            "kind": "string",
            "default_units": "",
            "enumeration": ("from monomers",),
            "format_string": "",
            "description": "Name the dimer system:",
            "help_text": (
                "The name for the new system holding the dimer configurations. "
                "'from monomers' builds it from the two monomer names; otherwise "
                "the literal text is used."
            ),
        },
        "configuration name": {
            "default": "orientation,distance",
            "kind": "string",
            "default_units": "",
            "enumeration": ("orientation,distance", "separation", "sequential"),
            "format_string": "",
            "description": "Name the configurations:",
            "help_text": (
                "How to name each generated configuration. 'orientation,distance' "
                "labels them by orientation and point index (1,1  1,2  ...  2,1  "
                "...) so the points of one scan group together; 'separation' "
                "labels them by their separation; 'sequential' numbers them 1, 2, "
                "... across all configurations. A comma (not '/') is used since "
                "'/' separates system and configuration names elsewhere."
            ),
        },
        "save scan variables as properties": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Save scan variables as properties:",
            "help_text": (
                "Whether to store the separation, gap, and orientation index on "
                "each configuration as properties, for downstream filtering."
            ),
        },
        "analysis plots": {
            "default": "basic",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("none", "basic", "detailed"),
            "format_string": "",
            "description": "Sampling diagnostics:",
            "help_text": (
                "Distribution diagnostics for the generated ensemble (separation "
                "coverage, contact distances, approach direction, relative "
                "orientation, and -- if interaction energies were computed -- the "
                "ΔE distribution and flatness). 'none' skips them. 'basic' prints "
                "the scalar summary and writes the combined interactive "
                "'dimer_sampling.graph' for the Dashboard. 'detailed' additionally "
                "writes each panel as its own 'dimer_sampling_<panel>.graph' for "
                "closer inspection. Any extra image formats (png/pdf/svg/...) are "
                "controlled by graph-formats in seamm.ini."
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("DimerBuilderParameters.__init__")

        super().__init__(
            defaults={
                **DimerBuilderParameters.parameters,
                **defaults,
            },
            data=data,
        )
