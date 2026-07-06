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
                "i.e. with a slight overlap."
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
            "enumeration": ("geometric", "linear", "explicit"),
            "format_string": "",
            "description": "Spacing:",
            "help_text": (
                "How the scan points are distributed along the gap coordinate. "
                "'geometric' clusters points near contact and thins them in the "
                "tail; 'linear' spaces them evenly; 'explicit' uses the list in "
                "'separations'."
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
                "'geometric' or 'linear' spacing."
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
