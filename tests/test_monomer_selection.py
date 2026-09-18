# -*- coding: utf-8 -*-

"""The monomer sources use the standard SEAMM structure selection, prefixed per
monomer, and old flowcharts are translated on load."""

import dimer_builder_step
import seamm
from dimer_builder_step.dimer_builder_parameters import monomer_selection_keys


def test_prefixed_standard_block():
    P = dimer_builder_step.DimerBuilderParameters()
    for prefix in ("monomer A", "monomer B"):
        for std_key, key in monomer_selection_keys(prefix).items():
            assert key in P
        assert P[f"{prefix} systems"].value == "current"
        assert P[f"{prefix} configurations"].value == "all"
        assert prefix not in P
    std = seamm.standard_parameters.structure_selection_parameters
    assert set(P["monomer A systems"].enumeration) == set(
        std["source systems"]["enumeration"]
    )


def test_legacy_translation():
    P = dimer_builder_step.DimerBuilderParameters()
    P.from_dict(
        {
            "monomer A": {"value": "current", "units": None},
            "monomer B": {"value": "waters", "units": None},
            "monomer B configurations": {"value": "last", "units": None},
        }
    )
    assert P["monomer A systems"].value == "current"
    assert P["monomer B systems"].value == "name is"
    assert P["monomer B system name"].value == "waters"
    assert P["monomer B configurations"].value == "last"
    P.from_dict({"monomer A": {"value": "$dimers", "units": None}})
    assert P["monomer A systems"].value == "$dimers"


def test_description_names_the_sources():
    node = dimer_builder_step.DimerBuilder()
    node._id = (1,)
    P = node.parameters.values_to_dict()
    P["monomer B systems"] = "name is"
    P["monomer B system name"] = "waters"
    P["monomer B configurations"] = "last"
    text = " ".join(node.description_text(P).split())
    assert "monomer A (all configurations of the current system)" in text
    assert "monomer B (the last configuration of the system 'waters')" in text
