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

    This is an initial, foundational set of parameters; the full set (input
    mode, separation schedule, contact-distance handling, engine on/off, output
    options, ...) is fleshed out as the step is implemented. The keys are the
    parameters for this plug-in; each value is a dictionary describing that
    parameter (default, kind, units, enumeration, format, description, help).

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
    """

    parameters = {
        "maximum separation": {
            "default": 10.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Maximum separation:",
            "help_text": ("The largest center-of-mass separation to scan out to."),
        },
        "number of orientations": {
            "default": 10,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of orientations:",
            "help_text": (
                "How many relative orientations of the two molecules to sample."
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
