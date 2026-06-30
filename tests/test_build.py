# -*- coding: utf-8 -*-

"""Tests for the headless dimer-building logic in DimerBuilder."""

import math

import numpy as np
import pytest

from molsystem.system_db import SystemDB
from seamm_util import Q_

import dimer_builder_step
from dimer_builder_step.dimer_builder import vdw_radii


def _add_water(configuration):
    r0 = 0.9572
    theta0 = 104.52
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))
    ids = configuration.atoms.append(
        x=[0.0, x, -x], y=[0.0, 0.0, 0.0], z=[0.0, z, z], atno=[8, 1, 1]
    )
    configuration.bonds.append(i=[ids[0], ids[0]], j=[ids[1], ids[2]], bondorder=[1, 1])
    return configuration


_db_counter = [0]


@pytest.fixture()
def db_two_waters():
    _db_counter[0] += 1
    db = SystemDB(filename=f"file:dimer_test_{_db_counter[0]}?mode=memory&cache=shared")
    a = db.create_system(name="A").create_configuration(name="w1")
    _add_water(a)
    b = db.create_system(name="B").create_configuration(name="w1")
    _add_water(b)

    yield db

    db.close()


def _P(**overrides):
    """A full parameter dict with sensible defaults for the build."""
    P = {
        "input mode": "two monomer sets",
        "monomer A": "A",
        "monomer A configurations": "all",
        "monomer A configuration name": "",
        "monomer B": "B",
        "monomer B configurations": "all",
        "monomer B configuration name": "",
        "number of orientations": 5,
        "random seed": "1",
        "contact method": "van der Waals radii",
        "innermost gap": Q_(-0.5, "Å"),
        "maximum separation": Q_(10.0, "Å"),
        "spacing": "geometric",
        "number of separations": 8,
        "separations": "",
        "system name": "from monomers",
        "configuration name": "sequential",
        "save scan variables as properties": "yes",
    }
    P.update(overrides)
    return P


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #


def test_vdw_radii_in_angstrom():
    r = vdw_radii(["O", "H", "C"])
    assert np.allclose(r, [1.52, 1.10, 1.70], atol=0.02)


def test_contact_distance_two_atoms_along_z():
    node = dimer_builder_step.DimerBuilder()
    A = np.array([[0.0, 0.0, 0.0]])
    B = np.array([[0.0, 0.0, 0.0]])
    axis = np.array([0.0, 0.0, 1.0])
    contact = node._contact_distance(A, np.array([1.5]), B, np.array([1.2]), axis)
    assert math.isclose(contact, 2.7, abs_tol=1.0e-9)


def test_contact_distance_offset_atoms():
    """Lateral offset reduces the along-axis contact distance."""
    node = dimer_builder_step.DimerBuilder()
    A = np.array([[0.0, 0.0, 0.0]])
    B = np.array([[1.0, 0.0, 0.0]])  # 1 Å lateral offset
    axis = np.array([0.0, 0.0, 1.0])
    R = 2.7
    contact = node._contact_distance(A, np.array([1.5]), B, np.array([1.2]), axis)
    assert math.isclose(contact, math.sqrt(R**2 - 1.0), abs_tol=1.0e-9)


def test_separation_schedule_geometric_range():
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    contact = 3.0
    d = node._separation_schedule(contact, P)
    assert len(d) == 8
    assert math.isclose(d.min(), contact - 0.5, abs_tol=1.0e-6)
    assert math.isclose(d.max(), 10.0, abs_tol=1.0e-6)
    assert np.all(np.diff(d) > 0)  # sorted, strictly increasing


def test_separation_schedule_explicit():
    node = dimer_builder_step.DimerBuilder()
    P = _P(spacing="explicit", separations="0.0, 1.0, 2.0")
    d = node._separation_schedule(5.0, P)
    assert np.allclose(d, [5.0, 6.0, 7.0])


# --------------------------------------------------------------------------- #
# Full build — Mode A
# --------------------------------------------------------------------------- #


def test_build_mode_a_counts_and_atomset(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    rng = np.random.default_rng(1)

    system, stats = node._build(db, P, rng)

    # 5 orientations x 8 separations
    assert stats["n_configurations"] == 40
    assert len(system.configurations) == 40
    # Every configuration is a dimer of 6 atoms ...
    assert all(c.n_atoms == 6 for c in system.configurations)
    # ... and they are all conformers sharing a single atomset.
    assert len({c.atomset for c in system.configurations}) == 1
    assert system.name == "A + B"


def test_build_mode_a_separation_range(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    P = _P()
    system, stats = node._build(db, P, np.random.default_rng(2))

    # Largest center-to-center separation reaches the maximum; smallest is a
    # slight overlap inside contact (but still positive).
    assert math.isclose(stats["max_separation"], 10.0, abs_tol=1.0e-6)
    assert 1.0 < stats["min_separation"] < 4.0


def test_build_mode_a_preserves_monomer_geometry(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(3))

    conf = system.configurations[0]
    xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
    # Atoms 0-2 are monomer A (water): O-H bonds ~0.9572 Å.
    oh = np.linalg.norm(xyz[1:3] - xyz[0], axis=1)
    assert np.allclose(oh, 0.9572, atol=1.0e-3)


def test_build_mode_a_monomer_a_is_fixed(db_two_waters):
    """With a single A conformer, monomer A is identical in every frame."""
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(7))

    reference = None
    for conf in system.configurations:
        xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
        A = xyz[:3]
        if reference is None:
            reference = A
            # A is centered at the origin.
            assert np.allclose(A.mean(axis=0), [0.0, 0.0, 0.0], atol=1.0e-9)
        else:
            assert np.allclose(A, reference, atol=1.0e-9)


def test_build_mode_a_tags_properties(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(4))

    conf = system.configurations[0]
    assert conf.properties.exists("dimer separation")
    assert conf.properties.exists("dimer gap")
    assert conf.properties.exists("dimer orientation")


def test_build_mode_a_no_severe_overlap(db_two_waters):
    """No atom pair across the two monomers is grossly inside vdW contact."""
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(5))

    rA = vdw_radii(["O", "H", "H"])
    for conf in system.configurations:
        xyz = np.asarray(conf.atoms.get_coordinates(fractionals=False, as_array=True))
        A = xyz[:3]
        B = xyz[3:]
        D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
        Rsum = rA[:, None] + rA[None, :]
        # Innermost point allows a slight overlap; require no pair closer than
        # 70% of the vdW sum.
        assert np.all(D > 0.7 * Rsum)


# --------------------------------------------------------------------------- #
# Full build — Mode B (prepared dimers)
# --------------------------------------------------------------------------- #


def test_build_mode_b_from_prepared_dimer(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    # A prepared water dimer, B displaced 3 Å along z.
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"

    system, stats = node._build(db, P, np.random.default_rng(6))

    assert stats["n_configurations"] == 8  # 1 dimer x 8 separations
    assert len(system.configurations) == 8
    assert all(c.n_atoms == 6 for c in system.configurations)
    assert len({c.atomset for c in system.configurations}) == 1
