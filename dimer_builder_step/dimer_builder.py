# -*- coding: utf-8 -*-

"""Non-graphical part of the Dimer Builder step in a SEAMM flowchart"""

import logging
import importlib.resources
import pprint  # noqa: F401

import numpy as np

import dimer_builder_step
import molsystem
from molsystem import random_rotation_matrix
import seamm
from seamm_util import getParser, parse_list, ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Dimer Builder")

# Add this module's properties to the standard properties
path = importlib.resources.files("dimer_builder_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


def vdw_radii(symbols):
    """Van der Waals radii, in Å, for a list of element symbols.

    Uses the mendeleev package (radii are tabulated in picometres). Falls back
    to the Pyykkö covalent radius, then to 1.5 Å, for any element lacking a vdW
    radius.

    Parameters
    ----------
    symbols : [str]
        The element symbols.

    Returns
    -------
    numpy.ndarray
        The van der Waals radii in Å, in the same order as ``symbols``.
    """
    import mendeleev

    elements = mendeleev.element(list(symbols))
    if not isinstance(elements, (list, tuple)):
        elements = [elements]

    radii = []
    for element in elements:
        r = element.vdw_radius
        if r is None:
            r = element.covalent_radius_pyykko
        if r is None:
            r = 150.0
        radii.append(r / 100.0)  # picometres -> Å
    return np.array(radii)


class DimerBuilder(seamm.Node):
    """
    The non-graphical part of a Dimer Builder step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : DimerBuilderParameters
        The control parameters for Dimer Builder.

    See Also
    --------
    TkDimerBuilder,
    DimerBuilder, DimerBuilderParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Dimer Builder",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for Dimer Builder in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Dimer Builder {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Dimer Builder", namespace=namespace
        )  # yapf: disable

        super().__init__(
            flowchart=flowchart,
            title="Dimer Builder",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = dimer_builder_step.metadata
        self.parameters = dimer_builder_step.DimerBuilderParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return dimer_builder_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return dimer_builder_step.__git_revision__

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "dimer-builder-step"
        parser = getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for the dimer builder itself can be added here.
            pass

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool
            If True, omit the description of the sub-flowchart steps.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        if P["input mode"] == "two monomer sets":
            text = (
                f"Build dimers from monomer A ({P['monomer A']}) and monomer B "
                f"({P['monomer B']}), sampling {P['number of orientations']} random "
                "relative orientations (each a random conformer of A and of B at a "
                "random orientation)."
            )
        else:
            text = (
                f"Scan the prepared structures from {P['monomer A']}, sliding the "
                "'movable' group away from and into the 'fixed' group along their "
                "center-to-center axis. The two groups are taken from 'fixed'/"
                "'movable' subsets if present, otherwise the last molecule is "
                "movable and the rest are fixed."
            )

        text += (
            f" For each, locate the contact distance using the '{P['contact method']}' "
            f"method, then scan the center-to-center separation from an innermost gap "
            f"of {P['innermost gap']} out to {P['maximum separation']} with "
            f"{P['number of separations']} points ({P['spacing']} spacing). The "
            f"configurations are stored in a new system named '{P['system name']}'."
        )

        result = self.header + "\n\n" + str(__(text, indent=4 * " "))

        node = self.subflowchart.get_node("1").next()
        if not short and node is not None:
            result += "\n\n"
            result += str(
                __(
                    "The following sub-flowchart is available for energy-based "
                    "options (e.g. the 'energy' contact method):",
                    indent=4 * " ",
                )
            )
            result += "\n\n"
            while node is not None:
                result += str(__(node.description_text(), indent=7 * " ", wrap=False))
                result += "\n"
                node = node.next()

        return result

    def run(self):
        """Run a Dimer Builder step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Print what we are about to do
        printer.important(__(self.description_text(P, short=True), indent=self.indent))
        printer.important("")

        system_db = self.get_variable("_system_db")
        rng = self._make_rng(P["random seed"])

        new_system, stats = self._build(system_db, P, rng)

        # The contact distance uses van der Waals radii from mendeleev. (The
        # plug-in's own citation is added automatically by the base class.)
        if "mendeleev" in self._bibliography:
            self.references.cite(
                raw=self._bibliography["mendeleev"],
                alias="mendeleev",
                module="dimer_builder_step",
                level=2,
                note="Van der Waals radii used to estimate the contact distance.",
            )

        # Make the new system & its first configuration current
        system_db.system = new_system
        new_system.configuration = new_system.configurations[0].id

        self.analyze(P=P, stats=stats)

        return next_node

    def analyze(self, P=None, stats=None, indent="", **kwargs):
        """Summarize what was generated."""
        if stats is None:
            return

        text = (
            f"Created {stats['n_configurations']} configurations in the system "
            f"'{stats['system']}' from {stats['n_seeds']} "
            f"{'orientation' if stats['n_seeds'] == 1 else 'orientations'}, with "
            f"center-to-center separations from {stats['min_separation']:.2f} to "
            f"{stats['max_separation']:.2f} Å."
        )
        printer.important(__(text, indent=4 * " "))
        printer.important("")

    # ----------------------------------------------------------------- #
    # Implementation helpers
    # ----------------------------------------------------------------- #

    def _make_rng(self, seed):
        """A numpy random generator from the 'random seed' parameter."""
        if isinstance(seed, str):
            if seed.strip() == "" or seed.strip().lower() == "random":
                return np.random.default_rng()
            seed = int(seed)
        return np.random.default_rng(int(seed))

    def _resolve_pool(self, spec, configurations, name, system_db):
        """Resolve an input specification to a list of configurations.

        ``spec`` is either an already-dereferenced list of configurations (from a
        ``$variable``), or a string: 'current', a system name, or '$variable'.
        ``configurations``/``name`` select within a system (ignored for a list).
        """
        # A variable holding a list of configurations: use all of them.
        if not isinstance(spec, str):
            return list(spec)

        spec = spec.strip()
        if spec.startswith("$"):
            value = self.get_variable(spec[1:])
            return list(value)

        if spec == "" or spec.lower() == "current":
            system = system_db.system
        else:
            system = system_db.get_system(spec)

        return self._select_configurations(system, configurations, name)

    def _select_configurations(self, system, how, name):
        """Pick configurations from a system, mirroring the loop step."""
        configurations = system.configurations
        if how == "all":
            return configurations
        elif how == "last":
            return [configurations[-1]]
        elif how == "first":
            return [configurations[0]]
        elif how == "name is":
            return [c for c in configurations if c.name == name]
        elif how == "name matches":
            import fnmatch

            return [c for c in configurations if fnmatch.fnmatch(c.name, name)]
        elif how == "name regexp":
            import re

            pattern = re.compile(name)
            return [c for c in configurations if pattern.search(c.name)]
        else:
            raise ValueError(f"Unknown configuration selector '{how}'.")

    def _contact_distance(self, A_xyz, A_radii, B_xyz, B_radii, axis):
        """The center-to-center distance at which A and B first touch.

        A is fixed with its center at the origin; B (center at the origin) is slid
        along ``axis`` to a center-to-center distance ``d``. This returns the
        largest ``d`` at which any atom pair is in van der Waals contact -- i.e.
        the onset of contact as B is brought in from infinity.
        """
        axis = np.asarray(axis, dtype=float)
        W = A_xyz[:, np.newaxis, :] - B_xyz[np.newaxis, :, :]  # (nA, nB, 3)
        wu = W @ axis  # (nA, nB)
        w2 = np.einsum("ijk,ijk->ij", W, W)
        R = A_radii[:, np.newaxis] + B_radii[np.newaxis, :]  # (nA, nB)
        disc = wu**2 - w2 + R**2
        mask = disc >= 0.0
        if not mask.any():
            # No pair can touch along this axis; fall back to the closest approach.
            return float(np.max(wu))
        return float(np.max(wu[mask] + np.sqrt(disc[mask])))

    def _separation_schedule(self, contact, P):
        """Center-to-center distances for the radial scan of one orientation.

        The scan coordinate is the gap beyond contact (0 = touching). Geometric
        spacing clusters points near contact; an explicit list is given as gaps.
        """
        inner_gap = P["innermost gap"].to("Å").magnitude
        max_sep = P["maximum separation"].to("Å").magnitude
        n = P["number of separations"]
        spacing = P["spacing"]

        if spacing == "explicit":
            gaps = np.array(parse_list(P["separations"]), dtype=float)
            distances = contact + gaps
        elif spacing == "linear":
            distances = np.linspace(contact + inner_gap, max_sep, n)
        else:  # geometric
            # Geometric spacing of the center-to-center distance over the whole
            # range, so steps grow smoothly with separation and cluster near
            # contact -- with no isolated, oversized first step.
            d_min = max(contact + inner_gap, 0.1)
            d_max = max_sep if max_sep > d_min else d_min + 0.5
            distances = np.geomspace(d_min, d_max, n)

        distances = np.unique(np.clip(distances, 0.1, None))
        return distances

    def _principal_axes(self, xyz, masses):
        """The center of mass and the principal-axis rotation of a molecule.

        ``axes`` is a proper rotation (det = +1) whose columns are the principal
        axes of inertia, ordered by ascending moment; centering ``xyz`` and
        rotating by ``axes`` puts the molecule in its principal-axis frame.
        """
        xyz = np.asarray(xyz, dtype=float)
        masses = np.asarray(masses, dtype=float)
        com = (masses[:, np.newaxis] * xyz).sum(axis=0) / masses.sum()
        centered = xyz - com

        inertia = np.zeros((3, 3))
        for m, r in zip(masses, centered):
            inertia += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))

        # Eigenvectors (columns) are the principal axes, eigenvalues ascending.
        _, axes = np.linalg.eigh(inertia)
        if np.linalg.det(axes) < 0.0:
            axes[:, 0] = -axes[:, 0]  # keep a proper rotation, not a reflection
        return com, axes

    def _orient_to_principal_axes(self, xyz, masses):
        """Coordinates centered at the COM and rotated onto the principal axes.

        A symmetric molecule's rotation axis then lands on a Cartesian axis. The
        transform is a proper rotation (never a reflection) so chirality is kept.
        """
        com, axes = self._principal_axes(xyz, masses)
        return (np.asarray(xyz, dtype=float) - com) @ axes

    def _direction_angles(self, axis):
        """Polar (theta) and azimuthal (phi) angles of a vector, in degrees."""
        axis = np.asarray(axis, dtype=float)
        theta = np.degrees(np.arccos(np.clip(axis[2], -1.0, 1.0)))
        phi = np.degrees(np.arctan2(axis[1], axis[0]))
        return float(theta), float(phi)

    def _euler_zyz(self, R):
        """ZYZ Euler angles (alpha, beta, gamma) of a rotation matrix, degrees."""
        R = np.asarray(R, dtype=float)
        beta = np.arctan2(np.hypot(R[0, 2], R[1, 2]), R[2, 2])
        if np.sin(beta) > 1.0e-8:
            alpha = np.arctan2(R[1, 2], R[0, 2])
            gamma = np.arctan2(R[2, 1], -R[2, 0])
        else:  # gimbal lock: fold the two z-rotations together
            alpha = np.arctan2(-R[0, 1], R[0, 0])
            gamma = 0.0
        return (
            float(np.degrees(alpha)),
            float(np.degrees(beta)),
            float(np.degrees(gamma)),
        )

    def _ensure_templates(self, system_db):
        """Get (creating if needed) the 'fixed' and 'movable' subset templates."""
        templates = system_db.templates
        result = []
        for name in ("fixed", "movable"):
            if templates.exists(name, "general"):
                result.append(templates.get(name, "general"))
            else:
                result.append(templates.create(name=name, category="general"))
        return result

    def _add_subsets(self, configuration, t_fixed, t_movable, fixed_ids, movable_ids):
        """Tag the fixed and movable atom groups as subsets on a configuration."""
        configuration.subsets.create(template=t_fixed, atoms=fixed_ids)
        configuration.subsets.create(template=t_movable, atoms=movable_ids)

    def _fixed_movable_indices(self, configuration):
        """The fixed and movable atom groups (0-based indices) for a structure.

        Uses 'fixed'/'movable' subsets on the configuration if present; otherwise
        the last molecule (by connectivity) is movable and the rest are fixed.
        """
        system_db = configuration.system_db
        atom_ids = configuration.atoms.ids
        index_of = {atom_id: i for i, atom_id in enumerate(atom_ids)}

        def subset_indices(name):
            if not system_db.templates.exists(name, "general"):
                return None
            subsets = configuration.subsets.get(
                system_db.templates.get(name, "general")
            )
            if not subsets:
                return None
            return sorted(index_of[aid] for s in subsets for aid in s.atoms.ids)

        movable_idx = subset_indices("movable")
        if movable_idx is None:
            molecules = configuration.find_molecules(as_indices=True)
            if len(molecules) < 2:
                raise ValueError(
                    "A prepared structure must contain at least two molecules (or "
                    "'fixed'/'movable' subsets) to scan."
                )
            movable_idx = molecules[-1]

        fixed_idx = subset_indices("fixed")
        if fixed_idx is None:
            movable_set = set(movable_idx)
            fixed_idx = [i for i in range(len(atom_ids)) if i not in movable_set]

        return fixed_idx, movable_idx

    # ----------------------------------------------------------------- #
    # Energy-based contact ('contact method' = 'energy'), via seamm_mdi
    # ----------------------------------------------------------------- #

    def _open_energy_engine(self, elements, charge, multiplicity):
        """Start an MDI engine for the dimer from the upstream model chemistry.

        Returns a started ``seamm_mdi.MDIEngine`` or raises with guidance if no
        MDI-capable model chemistry is available.
        """
        from seamm_mdi import MDIEngine

        try:
            mc = self.get_variable("_model_chemistry")
        except Exception:
            raise ValueError(
                "The 'energy' contact method needs a model chemistry. Add a "
                "Model Chemistry step before the Dimer Builder step."
            )
        options = mc.get("options", {}) if isinstance(mc, dict) else {}
        if not options.get("mdi_capable", False):
            raise ValueError(
                f"The model chemistry '{mc.get('level', mc)}' is not MDI-capable; "
                "the 'energy' contact method needs an MDI engine such as MOPAC or "
                "xTB."
            )

        step = self.flowchart.plugin_manager.get(mc["step"])
        executor = self.flowchart.executor
        seamm_options = self.global_options
        method = mc["method"]
        n_atoms = len(elements)

        def build_argv(hostname, port):
            return step.get_mdi_engine_command(
                executor,
                seamm_options,
                method=method,
                port=port,
                hostname=hostname,
                charge=charge,
                multiplicity=multiplicity,
                n_atoms=n_atoms,
            )

        engine = MDIEngine(build_argv, elements=elements, name="DimerBuilder")
        engine.start()
        return engine

    def _energy_anchor(self, engine, assemble, seed, P):
        """Locate the energy minimum along the approach axis (the scan anchor).

        ``assemble(d)`` returns the full dimer coordinates (Å) with the movable
        group placed at center-to-center distance ``d``. ``seed`` is the van der
        Waals contact estimate used to bracket the search.
        """
        max_sep = P["maximum separation"].to("Å").magnitude
        lo = max(seed - 0.5, 0.5)
        hi = min(seed + 5.0, max_sep)
        if hi <= lo:
            hi = lo + 1.0

        def energy_at(d):
            engine.set_coordinates(assemble(d), units="Å")
            return engine.energy(units="hartree")

        return self._minimize_on_grid(energy_at, lo, hi, 9)

    @staticmethod
    def _minimize_on_grid(func, lo, hi, n):
        """Minimize a 1-D function on a uniform grid, refining parabolically.

        Coarse but robust and derivative-free -- enough to anchor the scan.
        """
        ds = np.linspace(lo, hi, n)
        es = np.array([func(float(d)) for d in ds])
        k = int(np.argmin(es))
        if 0 < k < n - 1:
            e0, e1, e2 = es[k - 1], es[k], es[k + 1]
            denom = e0 - 2.0 * e1 + e2
            if denom > 0.0:  # concave up -> a real interior minimum
                h = ds[1] - ds[0]
                return float(ds[k] + 0.5 * h * (e0 - e2) / denom)
        return float(ds[k])

    def _build(self, system_db, P, rng):
        """Generate the dimer configurations. Returns (system, stats)."""
        if P["input mode"] == "two monomer sets":
            return self._build_from_monomers(system_db, P, rng)
        else:
            return self._build_from_dimers(system_db, P)

    def _build_from_monomers(self, system_db, P, rng):
        """Mode A: assemble dimers from two monomer conformer pools."""
        A_pool = self._resolve_pool(
            P["monomer A"],
            P["monomer A configurations"],
            P["monomer A configuration name"],
            system_db,
        )
        B_pool = self._resolve_pool(
            P["monomer B"],
            P["monomer B configurations"],
            P["monomer B configuration name"],
            system_db,
        )
        if len(A_pool) == 0 or len(B_pool) == 0:
            raise ValueError("Both monomer A and monomer B must supply structures.")

        A0, B0 = A_pool[0], B_pool[0]

        name = self._system_name(P, A0.system.name, B0.system.name)
        dimer_sys = system_db.create_combined_system([A0, B0], name=name)
        base = dimer_sys.configuration

        # Tag molecule A as 'fixed' and molecule B as 'movable'. The atom ids are
        # shared across all conformers (one atomset), so compute them once.
        nA = A0.n_atoms
        atom_ids = base.atoms.ids
        fixed_ids = atom_ids[:nA]
        movable_ids = atom_ids[nA:]
        t_fixed, t_movable = self._ensure_templates(system_db)

        A_radii = vdw_radii(A0.atoms.symbols)
        B_radii = vdw_radii(B0.atoms.symbols)

        engine = None
        if P["contact method"] == "energy":
            elements = list(A0.atoms.atomic_numbers) + list(B0.atoms.atomic_numbers)
            charge = (A0.charge or 0) + (B0.charge or 0)
            engine = self._open_energy_engine(elements, charge, 1)

        save_props = self._truthy(P["save scan variables as properties"])
        count = 0
        separations = []
        try:
            for orientation in range(1, P["number of orientations"] + 1):
                Ac = A_pool[int(rng.integers(len(A_pool)))]
                Bc = B_pool[int(rng.integers(len(B_pool)))]

                # Monomer A is held fixed: centered at the COM on its principal
                # axes, so it is identical across the whole scan (and across
                # orientations, for a single conformer) -- a clean visual anchor.
                xyzA = self._orient_to_principal_axes(
                    np.array(
                        Ac.atoms.get_coordinates(fractionals=False, as_array=True)
                    ),
                    Ac.atoms.atomic_masses,
                )

                # Monomer B carries the randomness: its principal-axis frame is
                # rotated by R_B, and it approaches from a random direction.
                R_B = random_rotation_matrix(rng)
                xyzB = (
                    self._orient_to_principal_axes(
                        np.array(
                            Bc.atoms.get_coordinates(fractionals=False, as_array=True)
                        ),
                        Bc.atoms.atomic_masses,
                    )
                    @ R_B.T
                )
                axis = rng.standard_normal(3)
                axis = axis / np.linalg.norm(axis)

                theta, phi = self._direction_angles(axis)
                alpha, beta, gamma = self._euler_zyz(R_B)

                contact = self._contact_distance(xyzA, A_radii, xyzB, B_radii, axis)
                if engine is not None:

                    def assemble(d, xyzA=xyzA, xyzB=xyzB, axis=axis):
                        return np.vstack([xyzA, xyzB + axis * d])

                    contact = self._energy_anchor(engine, assemble, contact, P)

                for point, d in enumerate(
                    self._separation_schedule(contact, P), start=1
                ):
                    coordinates = np.vstack([xyzA, xyzB + axis * d])
                    conf = base if count == 0 else dimer_sys.copy_configuration(base)
                    conf.atoms.set_coordinates(coordinates, fractionals=False)
                    geometry = {
                        "separation": d,
                        "gap": d - contact,
                        "theta": theta,
                        "phi": phi,
                        "alpha": alpha,
                        "beta": beta,
                        "gamma": gamma,
                    }
                    self._name_and_tag(
                        conf, P, count + 1, orientation, point, geometry, save_props
                    )
                    self._add_subsets(conf, t_fixed, t_movable, fixed_ids, movable_ids)
                    separations.append(d)
                    count += 1
        finally:
            if engine is not None:
                engine.close()

        return dimer_sys, self._stats(name, P["number of orientations"], separations)

    def _build_from_dimers(self, system_db, P):
        """Mode B: radial profiles from prepared complexes (fixed orientation).

        The 'movable' group is slid out from / in toward the 'fixed' group along
        their center-to-center axis, preserving the input relative orientation.
        The two groups come from 'fixed'/'movable' subsets on the input if they
        exist; otherwise the last molecule is movable and the rest are fixed.
        """
        pool = self._resolve_pool(
            P["monomer A"],
            P["monomer A configurations"],
            P["monomer A configuration name"],
            system_db,
        )
        if len(pool) == 0:
            raise ValueError("No prepared structures were found.")

        d0 = pool[0]
        fixed_idx, movable_idx = self._fixed_movable_indices(d0)

        name = self._system_name(P, d0.system.name, None)
        out_sys = system_db.create_combined_system([d0], name=name)
        base = out_sys.configuration

        radii = vdw_radii(d0.atoms.symbols)
        fixed_radii = radii[fixed_idx]
        movable_radii = radii[movable_idx]
        masses = np.asarray(d0.atoms.atomic_masses)
        movable_masses = masses[movable_idx]

        t_fixed, t_movable = self._ensure_templates(system_db)
        atom_ids = base.atoms.ids
        fixed_ids = [atom_ids[i] for i in fixed_idx]
        movable_ids = [atom_ids[i] for i in movable_idx]

        engine = None
        if P["contact method"] == "energy":
            engine = self._open_energy_engine(
                list(d0.atoms.atomic_numbers), d0.charge or 0, 1
            )

        save_props = self._truthy(P["save scan variables as properties"])
        count = 0
        separations = []
        orientation = 0
        try:
            for orientation, structure in enumerate(pool, start=1):
                if structure.n_atoms != d0.n_atoms:
                    self.logger.warning(
                        f"Skipping '{structure.name}': it has a different number of "
                        "atoms than the first structure (mixed compositions are not "
                        "supported)."
                    )
                    continue
                xyz = np.array(
                    structure.atoms.get_coordinates(fractionals=False, as_array=True),
                    dtype=float,
                )
                fixed_xyz = xyz[fixed_idx]
                movable_xyz = xyz[movable_idx]
                fixed_center = fixed_xyz.mean(axis=0)
                movable_center = movable_xyz.mean(axis=0)
                axis = movable_center - fixed_center
                distance0 = np.linalg.norm(axis)
                if distance0 < 1.0e-6:
                    continue
                axis = axis / distance0

                fixed_centered = fixed_xyz - fixed_center
                movable_centered = movable_xyz - movable_center
                contact = self._contact_distance(
                    fixed_centered, fixed_radii, movable_centered, movable_radii, axis
                )

                theta, phi = self._direction_angles(axis)
                _, movable_axes = self._principal_axes(movable_xyz, movable_masses)
                alpha, beta, gamma = self._euler_zyz(movable_axes)

                def assemble(
                    d,
                    fixed_centered=fixed_centered,
                    movable_centered=movable_centered,
                    axis=axis,
                    fixed_idx=fixed_idx,
                    movable_idx=movable_idx,
                    template=xyz,
                ):
                    c = np.empty_like(template)
                    c[fixed_idx] = fixed_centered
                    c[movable_idx] = movable_centered + axis * d
                    return c

                if engine is not None:
                    contact = self._energy_anchor(engine, assemble, contact, P)

                for point, d in enumerate(
                    self._separation_schedule(contact, P), start=1
                ):
                    conf = base if count == 0 else out_sys.copy_configuration(base)
                    conf.atoms.set_coordinates(assemble(d), fractionals=False)
                    geometry = {
                        "separation": d,
                        "gap": d - contact,
                        "theta": theta,
                        "phi": phi,
                        "alpha": alpha,
                        "beta": beta,
                        "gamma": gamma,
                    }
                    self._name_and_tag(
                        conf, P, count + 1, orientation, point, geometry, save_props
                    )
                    self._add_subsets(conf, t_fixed, t_movable, fixed_ids, movable_ids)
                    separations.append(d)
                    count += 1
        finally:
            if engine is not None:
                engine.close()

        return out_sys, self._stats(name, orientation, separations)

    def _system_name(self, P, name_a, name_b):
        """Resolve the output system name."""
        requested = P["system name"]
        if requested == "from monomers":
            if name_b is None:
                return f"{name_a} (scan)"
            return f"{name_a} + {name_b}"
        return requested

    def _name_and_tag(self, conf, P, index, orientation, point, geometry, save_props):
        """Name a generated configuration and optionally tag its geometry.

        ``index`` is the global running count, ``orientation`` the orientation
        number, and ``point`` the 1-based distance index within that orientation.
        ``geometry`` holds the scan coordinates for this configuration:
        ``separation`` and ``gap`` (Å), the approach-direction angles ``theta``
        and ``phi`` (degrees), and the movable group's ZYZ Euler angles
        ``alpha``/``beta``/``gamma`` (degrees).
        """
        naming = P["configuration name"]
        if naming == "separation":
            conf.name = f"{geometry['separation']:.2f} Å"
        elif naming == "sequential":
            conf.name = str(index)
        else:  # orientation/distance
            conf.name = f"{orientation}/{point}"

        if save_props:
            p = "#DimerBuilder#scan"
            self._put_property(
                conf, f"dimer separation{p}", geometry["separation"], "Å"
            )
            self._put_property(conf, f"dimer gap{p}", geometry["gap"], "Å")
            self._put_property(
                conf, f"dimer orientation{p}", int(orientation), None, _type="int"
            )
            self._put_property(conf, f"approach theta{p}", geometry["theta"], "degree")
            self._put_property(conf, f"approach phi{p}", geometry["phi"], "degree")
            self._put_property(conf, f"movable alpha{p}", geometry["alpha"], "degree")
            self._put_property(conf, f"movable beta{p}", geometry["beta"], "degree")
            self._put_property(conf, f"movable gamma{p}", geometry["gamma"], "degree")

    def _put_property(self, conf, name, value, units, _type="float"):
        """Store a property on a configuration (defining it if not registered).

        The scan properties are normally pre-registered from data/properties.csv;
        the fallback definition here covers a configuration whose database lacks
        that registration.
        """
        properties = conf.properties
        if not properties.exists(name):
            properties.add(name, _type, units=units, noerror=True)
        properties.put(name, value)

    def _stats(self, name, n_seeds, separations):
        return {
            "system": name,
            "n_seeds": n_seeds,
            "n_configurations": len(separations),
            "min_separation": min(separations) if separations else 0.0,
            "max_separation": max(separations) if separations else 0.0,
        }

    @staticmethod
    def _truthy(value):
        return value is True or (isinstance(value, str) and value.lower() == "yes")
