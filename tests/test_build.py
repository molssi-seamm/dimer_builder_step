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
        "configuration name": "orientation/distance",
        "save scan variables as properties": "yes",
    }
    P.update(overrides)
    return P


def _subset_atom_ids(db, conf, name):
    """All atom ids in the named subsets of a configuration (or None)."""
    if not db.templates.exists(name, "general"):
        return None
    template = db.templates.get(name, "general")
    ids = []
    for subset in conf.subsets.get(template):
        ids.extend(subset.atoms.ids)
    return ids


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
    steps = np.diff(d)
    assert np.all(steps > 0)  # sorted, strictly increasing
    # Steps grow smoothly with separation: no oversized first step.
    assert np.all(np.diff(steps) >= -1.0e-9)


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


def test_build_mode_a_orientation_distance_names(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    # 3 orientations x 4 separations, default 'orientation/distance' naming.
    P = _P(spacing="explicit", separations="0.0, 1.0, 2.0, 3.0")
    P["number of orientations"] = 3
    system, stats = node._build(db, P, np.random.default_rng(12))

    names = [c.name for c in system.configurations]
    assert names[:4] == ["1/1", "1/2", "1/3", "1/4"]
    assert names[4:8] == ["2/1", "2/2", "2/3", "2/4"]
    assert names[-1] == "3/4"


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
            # A is centered at its center of mass at the origin.
            masses = np.asarray(conf.atoms.atomic_masses)[:3]
            com = (masses[:, None] * A).sum(axis=0) / masses.sum()
            assert np.allclose(com, [0.0, 0.0, 0.0], atol=1.0e-9)
        else:
            assert np.allclose(A, reference, atol=1.0e-9)


class _HarmonicEngine:
    """A fake engine whose energy is harmonic in the two fragments' separation.

    Duck-types the bits of seamm_mdi.MDIEngine that _energy_anchor uses, so the
    energy-contact machinery can be tested without MDI or a real code.
    """

    def __init__(self, nA, r0):
        self.nA = nA
        self.r0 = r0
        self._xyz = None

    def set_coordinates(self, xyz, units="bohr"):
        self._xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)

    def energy(self, units="hartree"):
        a = self._xyz[: self.nA].mean(axis=0)
        b = self._xyz[self.nA :].mean(axis=0)
        r = np.linalg.norm(b - a)
        return float((r - self.r0) ** 2)


def test_minimize_on_grid_parabola():
    node = dimer_builder_step.DimerBuilder()
    f = lambda d: (d - 2.7) ** 2 + 1.0  # noqa: E731
    d_min, k, n = node._minimize_on_grid(f, 1.0, 5.0, 9)
    assert d_min == pytest.approx(2.7, abs=0.05)
    assert 0 < k < n - 1  # interior minimum


def _assemble_along_z(A, B):
    def assemble(d):
        return np.vstack([A, B + np.array([0.0, 0.0, d])])

    return assemble


def test_energy_anchor_finds_minimum():
    node = dimer_builder_step.DimerBuilder()
    engine = _HarmonicEngine(nA=3, r0=3.2)
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    anchor = node._energy_anchor(engine, assemble, seed=2.5, P=_P())
    assert anchor == pytest.approx(3.2, abs=0.05)


def test_energy_anchor_falls_back_when_no_well():
    # Monotonically decreasing energy (no binding well) -> anchor at the seed.
    class _Repulsive:
        def set_coordinates(self, xyz, units="bohr"):
            self._xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)

        def energy(self, units="hartree"):
            a = self._xyz[:3].mean(axis=0)
            b = self._xyz[3:].mean(axis=0)
            return -float(np.linalg.norm(b - a))  # keeps falling as they separate

    node = dimer_builder_step.DimerBuilder()
    assemble = _assemble_along_z(np.zeros((3, 3)), np.zeros((3, 3)))
    anchor = node._energy_anchor(_Repulsive(), assemble, seed=2.9, P=_P())
    assert anchor == pytest.approx(2.9)


def test_direction_angles_known():
    node = dimer_builder_step.DimerBuilder()
    assert node._direction_angles([0.0, 0.0, 1.0]) == (0.0, 0.0)
    th, ph = node._direction_angles([1.0, 0.0, 0.0])
    assert math.isclose(th, 90.0) and math.isclose(ph, 0.0)
    th, ph = node._direction_angles([0.0, 1.0, 0.0])
    assert math.isclose(th, 90.0) and math.isclose(ph, 90.0)


def test_euler_zyz_identity_and_roundtrip():
    node = dimer_builder_step.DimerBuilder()
    # Identity -> all zero.
    assert node._euler_zyz(np.eye(3)) == (0.0, 0.0, 0.0)

    # Round-trip: build R = Rz(a) Ry(b) Rz(c), extract, rebuild, compare.
    def Rz(t):
        c, s = math.cos(t), math.sin(t)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

    def Ry(t):
        c, s = math.cos(t), math.sin(t)
        return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])

    a, b, c = math.radians(40), math.radians(70), math.radians(20)
    R = Rz(a) @ Ry(b) @ Rz(c)
    alpha, beta, gamma = node._euler_zyz(R)
    R2 = Rz(math.radians(alpha)) @ Ry(math.radians(beta)) @ Rz(math.radians(gamma))
    assert np.allclose(R, R2, atol=1.0e-9)


def test_build_mode_a_tags_geometry_properties(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(11))

    conf = system.configurations[0]
    for name in (
        "dimer separation",
        "dimer gap",
        "dimer orientation",
        "approach theta",
        "approach phi",
        "movable alpha",
        "movable beta",
        "movable gamma",
    ):
        assert conf.properties.exists(f"{name}#DimerBuilder#scan")


def test_orient_to_principal_axes_diagonalizes_inertia():
    """After reorientation the inertia tensor is diagonal (axes on x/y/z)."""
    node = dimer_builder_step.DimerBuilder()
    # A tilted water, arbitrary placement.
    r0, theta0 = 0.9572, 104.52
    x = r0 * math.sin(math.radians(theta0 / 2))
    z = r0 * math.cos(math.radians(theta0 / 2))
    xyz = np.array([[0.0, 0.0, 0.0], [x, 0.0, z], [-x, 0.0, z]]) + 5.0
    R = node._orient_to_principal_axes(xyz, [15.999, 1.008, 1.008])

    masses = np.array([15.999, 1.008, 1.008])
    com = (masses[:, None] * R).sum(axis=0) / masses.sum()
    inertia = np.zeros((3, 3))
    for m, r in zip(masses, R - com):
        inertia += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    off_diag = inertia - np.diag(np.diag(inertia))
    assert np.allclose(off_diag, 0.0, atol=1.0e-9)
    # COM at the origin, and a proper rotation preserves bond lengths.
    assert np.allclose(com, 0.0, atol=1.0e-9)
    assert np.allclose(np.linalg.norm(R[1:] - R[0], axis=1), 0.9572, atol=1.0e-6)


def test_build_mode_a_tags_properties(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(4))

    conf = system.configurations[0]
    assert conf.properties.exists("dimer separation#DimerBuilder#scan")
    assert conf.properties.exists("dimer gap#DimerBuilder#scan")
    assert conf.properties.exists("dimer orientation#DimerBuilder#scan")


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


def test_build_mode_a_creates_fixed_movable_subsets(db_two_waters):
    db = db_two_waters
    node = dimer_builder_step.DimerBuilder()
    system, stats = node._build(db, _P(), np.random.default_rng(8))

    conf = system.configurations[0]
    atom_ids = conf.atoms.ids
    # A is 'fixed' (first 3 atoms), B is 'movable' (last 3).
    assert set(_subset_atom_ids(db, conf, "fixed")) == set(atom_ids[:3])
    assert set(_subset_atom_ids(db, conf, "movable")) == set(atom_ids[3:])
    # Every configuration carries both subsets.
    for c in system.configurations:
        assert len(_subset_atom_ids(db, c, "fixed")) == 3
        assert len(_subset_atom_ids(db, c, "movable")) == 3


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


def test_build_mode_b_default_last_molecule_movable(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"
    system, stats = node._build(db, P, np.random.default_rng(9))

    conf = system.configurations[0]
    out_ids = conf.atoms.ids
    # No user subsets -> last molecule (atoms 3-5) is movable, rest fixed.
    assert set(_subset_atom_ids(db, conf, "movable")) == set(out_ids[3:])
    assert set(_subset_atom_ids(db, conf, "fixed")) == set(out_ids[:3])


def test_build_mode_b_honors_user_subsets(db_two_waters):
    db = db_two_waters
    wA = db.get_system("A").configuration
    wB = db.get_system("B").configuration
    M = np.eye(4)
    M[:3, 3] = [0.0, 0.0, 3.0]
    dimer_sys = db.create_combined_system([wA, wB], transforms=[None, M], name="dimer")
    d0 = dimer_sys.configuration

    # The user designates the FIRST molecule as movable (overriding the default).
    for nm in ("fixed", "movable"):
        if not db.templates.exists(nm, "general"):
            db.templates.create(name=nm, category="general")
    ids = d0.atoms.ids
    d0.subsets.create(template=db.templates.get("movable", "general"), atoms=ids[:3])
    d0.subsets.create(template=db.templates.get("fixed", "general"), atoms=ids[3:])

    node = dimer_builder_step.DimerBuilder()
    P = _P()
    P["input mode"] = "prepared dimers"
    P["monomer A"] = "dimer"
    system, stats = node._build(db, P, np.random.default_rng(10))

    conf = system.configurations[0]
    out_ids = conf.atoms.ids
    # Movable follows the user's choice: the first molecule.
    assert set(_subset_atom_ids(db, conf, "movable")) == set(out_ids[:3])
