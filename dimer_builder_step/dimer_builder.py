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
                f"Scan the prepared dimers from {P['monomer A']}, splitting each into "
                "its two molecules by connectivity and scanning along their "
                "center-to-center axis."
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
            g_max = max_sep - contact
            if g_max <= 0.1:
                distances = np.array([contact + inner_gap, contact + 0.5])
            else:
                positive = np.geomspace(0.05, g_max, max(n - 1, 1))
                gaps = np.concatenate(([inner_gap], positive))
                distances = contact + gaps

        distances = np.unique(np.clip(distances, 0.1, None))
        return distances

    def _orient_to_principal_axes(self, xyz, masses):
        """Center at the center of mass and rotate onto the principal axes.

        Returns coordinates whose center of mass is at the origin and whose
        principal axes of inertia are aligned with x/y/z, so a symmetric
        molecule's rotation axis lands on a Cartesian axis. The transform is a
        proper rotation (never a reflection), so chirality is preserved.
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
        return centered @ axes

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

        A_radii = vdw_radii(A0.atoms.symbols)
        B_radii = vdw_radii(B0.atoms.symbols)

        save_props = self._truthy(P["save scan variables as properties"])
        count = 0
        separations = []
        for orientation in range(1, P["number of orientations"] + 1):
            Ac = A_pool[int(rng.integers(len(A_pool)))]
            Bc = B_pool[int(rng.integers(len(B_pool)))]

            # Monomer A is held fixed: centered at the origin in its input
            # orientation, so it is identical across the whole scan (and across
            # orientations, for a single conformer) -- a clean visual anchor.
            xyzA = np.array(
                Ac.atoms.get_coordinates(fractionals=False, as_array=True), dtype=float
            )
            xyzA = self._orient_to_principal_axes(xyzA, Ac.atoms.atomic_masses)

            # Monomer B carries all the randomness: a random orientation, and a
            # random approach direction along which it is scanned in and out.
            xyzB = np.array(
                Bc.atoms.get_coordinates(fractionals=False, as_array=True), dtype=float
            )
            xyzB = (xyzB - xyzB.mean(axis=0)) @ random_rotation_matrix(rng).T
            axis = rng.standard_normal(3)
            axis = axis / np.linalg.norm(axis)

            contact = self._contact_distance(xyzA, A_radii, xyzB, B_radii, axis)
            for d in self._separation_schedule(contact, P):
                coordinates = np.vstack([xyzA, xyzB + axis * d])
                conf = base if count == 0 else dimer_sys.copy_configuration(base)
                conf.atoms.set_coordinates(coordinates, fractionals=False)
                self._name_and_tag(
                    conf, P, count + 1, orientation, d, d - contact, save_props
                )
                separations.append(d)
                count += 1

        return dimer_sys, self._stats(name, P["number of orientations"], separations)

    def _build_from_dimers(self, system_db, P):
        """Mode B: radial profiles from prepared dimers (fixed orientation)."""
        pool = self._resolve_pool(
            P["monomer A"],
            P["monomer A configurations"],
            P["monomer A configuration name"],
            system_db,
        )
        if len(pool) == 0:
            raise ValueError("No prepared dimers were found.")

        d0 = pool[0]
        molecules = d0.find_molecules(as_indices=True)
        if len(molecules) != 2:
            raise ValueError(
                "Each prepared dimer must contain exactly two molecules; found "
                f"{len(molecules)}."
            )
        idxA, idxB = molecules[0], molecules[1]

        name = self._system_name(P, d0.system.name, None)
        out_sys = system_db.create_combined_system([d0], name=name)
        base = out_sys.configuration

        radii = vdw_radii(d0.atoms.symbols)
        A_radii, B_radii = radii[idxA], radii[idxB]

        save_props = self._truthy(P["save scan variables as properties"])
        count = 0
        separations = []
        for orientation, dimer in enumerate(pool, start=1):
            if dimer.n_atoms != d0.n_atoms:
                self.logger.warning(
                    f"Skipping '{dimer.name}': it has a different number of atoms "
                    "than the first dimer (mixed compositions are not supported)."
                )
                continue
            xyz = np.array(
                dimer.atoms.get_coordinates(fractionals=False, as_array=True),
                dtype=float,
            )
            fragA = xyz[idxA]
            fragB = xyz[idxB]
            comA = fragA.mean(axis=0)
            comB = fragB.mean(axis=0)
            axis = comB - comA
            distance0 = np.linalg.norm(axis)
            if distance0 < 1.0e-6:
                continue
            axis = axis / distance0

            A_centered = fragA - comA
            B_centered = fragB - comB
            contact = self._contact_distance(
                A_centered, A_radii, B_centered, B_radii, axis
            )

            for d in self._separation_schedule(contact, P):
                coordinates = np.empty_like(xyz)
                coordinates[idxA] = A_centered
                coordinates[idxB] = B_centered + axis * d
                conf = base if count == 0 else out_sys.copy_configuration(base)
                conf.atoms.set_coordinates(coordinates, fractionals=False)
                self._name_and_tag(
                    conf, P, count + 1, orientation, d, d - contact, save_props
                )
                separations.append(d)
                count += 1

        return out_sys, self._stats(name, orientation, separations)

    def _system_name(self, P, name_a, name_b):
        """Resolve the output system name."""
        requested = P["system name"]
        if requested == "from monomers":
            if name_b is None:
                return f"{name_a} (scan)"
            return f"{name_a} + {name_b}"
        return requested

    def _name_and_tag(self, conf, P, index, orientation, separation, gap, save_props):
        """Name a generated configuration and optionally tag scan variables."""
        if P["configuration name"] == "separation":
            conf.name = f"{separation:.2f} Å"
        else:
            conf.name = str(index)

        if save_props:
            self._put_property(conf, "dimer separation", separation, "Å")
            self._put_property(conf, "dimer gap", gap, "Å")
            self._put_property(conf, "dimer orientation", float(orientation), None)

    def _put_property(self, conf, name, value, units):
        """Define (once) and store a float property on a configuration."""
        properties = conf.properties
        if not properties.exists(name):
            properties.add(name, "float", units=units, noerror=True)
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
