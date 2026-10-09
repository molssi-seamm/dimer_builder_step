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


def _db_with(systems):
    """A system database holding the named systems, each with conformers given
    as lists of atomic numbers."""
    from molsystem.system_db import SystemDB

    db = SystemDB(filename="file:monomer_pools?mode=memory&cache=shared")
    for name, molecules in systems:
        system = db.create_system(name=name)
        for numbers in molecules:
            configuration = system.create_configuration()
            configuration.atoms.append(
                x=[float(i) for i in range(len(numbers))],
                y=[0.0] * len(numbers),
                z=[0.0] * len(numbers),
                atno=numbers,
            )
    return db


def _P(**monomer_a):
    node = dimer_builder_step.DimerBuilder()
    P = node.parameters.values_to_dict()
    P.update({f"monomer A {k}": v for k, v in monomer_a.items()})
    return node, P


def test_name_is_must_match_one_system():
    """Two systems with the same name made the pool ambiguous, and the build
    died in molsystem with an IndexError (science, ChemAI 5238/5240)."""
    import pytest

    water = [8, 1, 1]
    db = _db_with([("dup", [water]), ("dup", [water]), ("only", [water, water])])
    node, P = _P(systems="name is", **{"system name": "dup"})
    with pytest.raises(ValueError, match="2 systems are named 'dup'"):
        node._resolve_pool(P, "monomer A", db)
    node, P = _P(systems="name is", **{"system name": "absent"})
    with pytest.raises(ValueError, match="no system is named 'absent'"):
        node._resolve_pool(P, "monomer A", db)
    node, P = _P(systems="name is", **{"system name": "only"})
    assert len(node._resolve_pool(P, "monomer A", db)) == 2
    db.close()


def test_a_pool_is_one_molecule():
    import pytest

    db = _db_with([("water", [[8, 1, 1]]), ("methane", [[6, 1, 1, 1, 1]])])
    node, P = _P(systems="all")
    with pytest.raises(ValueError, match="not all the same molecule"):
        node._resolve_pool(P, "monomer A", db)
    db.close()
