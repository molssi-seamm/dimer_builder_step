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
            parser.add_argument(
                parser_name,
                "--graph-formats",
                default=tuple(),
                choices=("html", "png", "jpeg", "webp", "svg", "pdf"),
                nargs="+",
                help="extra formats to write for the sampling-diagnostics graph",
            )
            parser.add_argument(
                parser_name,
                "--graph-width",
                default=1024,
                help="width of graphs in formats that support it, defaults to 1024",
            )
            parser.add_argument(
                parser_name,
                "--graph-height",
                default=1024,
                help="height of graphs in formats that support it, defaults to 1024",
            )

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

        if P["spacing"] == "energy-stratified":
            text += (
                f" For each, evaluate the interaction energy ΔE(R) along the "
                f"approach using the '{P['contact method']}' method out to "
                f"{P['maximum separation']}, then pool the candidate points across "
                f"all orientations and keep about {P['target configurations']} of "
                f"them"
            )
            method = P.get("selection method", "energy bins")
            if method == "descriptor diversity":
                text += (
                    " by clustering in a space of geometric collective variables "
                    "plus the interaction energy and keeping one per cluster (a "
                    "DIRECT-style, de-duplicated, diverse set)."
                )
            elif method == "energy bins + diversity":
                text += (
                    f" by sorting them into {P['number of energy bins']} ΔE bins "
                    "(flat in interaction energy) and keeping a geometrically "
                    "diverse subset of each bin (DIRECT clustering)."
                )
            else:
                text += (
                    f" by sorting them into {P['number of energy bins']} ΔE bins "
                    f"over the range set by '{P['energy levels']}' and capping each "
                    f"bin equally (flat in interaction energy)."
                )
            if self._truthy(P.get("tail coverage", "no")):
                text += (
                    f" A long-range distance-coverage floor from "
                    f"{P['tail minimum separation']} to {P['maximum separation']} "
                    f"and {P['asymptote anchors']} asymptote anchors out to "
                    f"{P['anchor separation']} are added to sample the weak tail "
                    f"and pin the energy to zero at large separation."
                )
        else:
            text += (
                f" For each, locate the contact distance using the "
                f"'{P['contact method']}' method, then scan the center-to-center "
                f"separation from an innermost gap of {P['innermost gap']} out to "
                f"{P['maximum separation']} with {P['number of separations']} points "
                f"({P['spacing']} spacing)."
            )
        text += (
            f" The configurations are stored in a new system named "
            f"'{P['system name']}'."
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

        if stats.get("model_chemistry"):
            text = (
                f"The contact distances were found from the energy using the model "
                f"chemistry '{stats['model_chemistry']}', which was evaluated "
                f"{stats['n_energy_calls']} times."
            )
            printer.important(__(text, indent=4 * " "))

        level = P.get("analysis plots", "none") if P else "none"
        ensemble = stats.get("ensemble")
        if level != "none" and ensemble:
            self._run_diagnostics(ensemble, level, stats.get("system", "dimers"))

        printer.important("")

    # ----------------------------------------------------------------- #
    # Sampling diagnostics (vendored dimer_analysis module)
    # ----------------------------------------------------------------- #

    @staticmethod
    def _collect_ensemble(P):
        """Whether the build should accumulate the ensemble for diagnostics."""
        return (P.get("analysis plots", "none") if P else "none") != "none"

    @staticmethod
    def _make_dimer_record(symbols_A, xyz_A, symbols_B, xyz_B, geometry, label):
        """Build a ``dimer_analysis.Dimer`` for one generated configuration.

        The interaction energy (kJ/mol, ΔE = E_dimer − (E_A + E_B), attractive
        negative) is carried only when it was computed (the 'energy' contact
        method); the diagnostics drop the energy panels when it is absent.
        """
        from dimer_builder_step import dimer_analysis

        return dimer_analysis.Dimer(
            symbols_A=list(symbols_A),
            xyz_A=np.asarray(xyz_A, dtype=float),
            symbols_B=list(symbols_B),
            xyz_B=np.asarray(xyz_B, dtype=float),
            energy=geometry.get("interaction_energy"),
            separation=geometry.get("separation"),
            orientation=geometry.get("orientation"),
            label=label,
        )

    def _run_diagnostics(self, ensemble, level, title):
        """Print the scalar summary and write the diagnostics graph(s).

        Uses the vendored, framework-free ``dimer_analysis`` module: the metrics
        are numpy-only and the plotting helpers return plotly ``go.Figure``
        objects, written as SEAMM ``.graph`` files (plotly JSON) for the
        Dashboard, plus any extra image formats requested via ``graph-formats``
        in seamm.ini. 'basic' writes the combined dashboard; 'detailed' also
        writes each panel as its own graph for closer inspection. Diagnostics
        are best-effort and never abort the run.
        """
        from dimer_builder_step import dimer_analysis

        try:
            metrics = dimer_analysis.compute_metrics(ensemble)
            s = dimer_analysis.summarize(metrics)
        except Exception as e:
            printer.important(
                __(f"Could not compute the sampling diagnostics: {e}", indent=4 * " ")
            )
            return

        text = (
            f"Sampling diagnostics ({s['n']} dimers): COM separation "
            f"{s['R_min']:.2f}-{s['R_max']:.2f} Å, closest-contact median "
            f"{s['min_contact_median']:.2f} Å, approach concentration "
            f"{s['approach_concentration']:.2f} (0 = isotropic)."
        )
        printer.important(__(text, indent=4 * " "))
        if "energy_flatness" in s:
            text = (
                f"ΔE from {s['energy_min']:.1f} kJ/mol, "
                f"{100 * s['energy_frac_attractive']:.0f}% attractive, "
                f"energy-flatness CV {s['energy_flatness']:.2f} "
                f"(0 = perfectly flat-in-energy)."
            )
            printer.important(__(text, indent=4 * " "))

        import os

        os.makedirs(self.directory, exist_ok=True)

        # The combined dashboard (always). A .graph is plotly JSON ({data,
        # layout}); the Dashboard renders it interactively.
        try:
            figure = dimer_analysis.make_dashboard(metrics, title=str(title))
        except Exception as e:  # plotting is best-effort, never fatal
            printer.important(
                __(f"Could not build the diagnostics dashboard: {e}", indent=4 * " ")
            )
            return
        base = os.path.join(self.directory, "dimer_sampling")
        with open(base + ".graph", "w") as fd:
            fd.write(self._graph_json(figure))
        self._write_extra_graph_formats(figure, base)

        if level != "detailed":
            printer.important(
                __(
                    "Wrote the sampling dashboard 'dimer_sampling.graph'.",
                    indent=4 * " ",
                )
            )
            return

        # 'detailed': each panel as its own graph, from the same trace builders
        # as the dashboard (so combined and separated views never diverge).
        try:
            panels = dimer_analysis.make_panels(metrics)
        except Exception as e:
            printer.important(
                __(f"Could not build the individual panels: {e}", indent=4 * " ")
            )
            panels = {}
        for name, panel in panels.items():
            if panel is None:
                continue
            panel_base = f"{base}_{name}"
            with open(panel_base + ".graph", "w") as fd:
                fd.write(self._graph_json(panel))
            self._write_extra_graph_formats(panel, panel_base)
        printer.important(
            __(
                f"Wrote the sampling dashboard and {len(panels)} panel graphs "
                "('dimer_sampling*.graph').",
                indent=4 * " ",
            )
        )

    def _write_extra_graph_formats(self, figure, base):
        """Write the extra image formats requested by 'graph-formats' in seamm.ini.

        Mirrors the LAMMPS step. Each format is best-effort: a missing 'kaleido'
        (needed for static images) is reported once, not fatal.
        """
        import shlex

        options = getattr(self, "options", {}) or {}
        formats = options.get("graph_formats", ())
        if isinstance(formats, str):
            formats = shlex.split(formats)
        if not formats:
            return
        width = int(options.get("graph_width", 1024))
        height = int(options.get("graph_height", 1024))
        for fmt in formats:
            path = f"{base}.{fmt}"
            try:
                if fmt in ("html", "htm"):
                    figure.write_html(path)
                else:
                    figure.write_image(path, format=fmt, width=width, height=height)
            except Exception as e:
                printer.important(
                    __(
                        f"Could not write the diagnostics graph as '{fmt}': {e}",
                        indent=4 * " ",
                    )
                )

    @staticmethod
    def _graph_json(figure):
        """Serialize a plotly figure to SEAMM .graph JSON with plain-list arrays.

        plotly's ``to_json`` base64-encodes numpy arrays as ``{"dtype","bdata"}``
        typed arrays; the Dashboard's bundled plotly.js does not decode that form,
        so the axes render but the traces are empty. Expand any typed array back
        to a plain JSON list (the .graph is the only consumer that goes through
        the Dashboard's plotly.js; the png/pdf/html exports use plotly directly
        and handle their own encoding).
        """
        import base64
        import json

        def expand(obj):
            if isinstance(obj, dict):
                if "bdata" in obj and "dtype" in obj and isinstance(obj["bdata"], str):
                    arr = np.frombuffer(
                        base64.b64decode(obj["bdata"]), dtype=np.dtype(obj["dtype"])
                    )
                    shape = obj.get("shape")
                    if shape:
                        arr = arr.reshape([int(s) for s in str(shape).split(",")])
                    return arr.tolist()
                return {k: expand(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [expand(v) for v in obj]
            return obj

        return json.dumps(expand(json.loads(figure.to_json())))

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

    @staticmethod
    def _mdi_method_and_basis(mc, options):
        """The (method, basis) to launch the MDI engine with.

        ``method`` is the engine's real keyword -- ``options['mdi_method_arg']``,
        which un-aliases an ORCA functional whose '/' was replaced in the
        model-chemistry string -- falling back to the parsed ``mc['method']``.
        ``basis`` is ``options['mdi_basis_arg']``, or ``None`` for engines that
        take a method alone (MOPAC, xTB) and whose ``get_mdi_engine_command`` has
        no basis argument, so a basis must not be passed to them.
        """
        method = options.get("mdi_method_arg") or mc.get("method")
        basis = options.get("mdi_basis_arg")
        return method, basis

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

        self._energy_model = mc.get("level", mc.get("method", "the model chemistry"))
        step = self.flowchart.plugin_manager.get(mc["step"])
        executor = self.flowchart.executor
        seamm_options = self.global_options
        n_atoms = len(elements)
        method, basis = self._mdi_method_and_basis(mc, options)

        def build_argv(hostname, port):
            kwargs = {
                "method": method,
                "port": port,
                "hostname": hostname,
                "charge": charge,
                "multiplicity": multiplicity,
                "n_atoms": n_atoms,
            }
            if basis is not None:
                kwargs["basis"] = basis
            return step.get_mdi_engine_command(executor, seamm_options, **kwargs)

        engine = MDIEngine(build_argv, elements=elements, name="DimerBuilder")
        engine.start()
        return engine

    def _energy_anchor(self, engine, assemble, seed, P):
        """Locate the energy minimum along the approach axis (the scan anchor).

        ``assemble(d)`` returns the full dimer coordinates (Å) with the movable
        group placed at center-to-center distance ``d``. ``seed`` is the van der
        Waals contact estimate used to bracket the search.
        """
        # The energy minimum for a bound orientation sits at or just beyond the
        # vdW contact, so search a modest window around the seed at a fine step;
        # a minimum pinned at the outer edge means "no well in range".
        max_sep = P["maximum separation"].to("Å").magnitude
        lo = max(seed - 0.5, 0.5)
        hi = min(seed + 3.0, max_sep)
        if hi <= lo:
            hi = lo + 1.0

        def energy_at(d):
            engine.set_coordinates(assemble(d), units="Å")
            return engine.energy(units="hartree")

        d_min, k, n = self._minimize_on_grid(energy_at, lo, hi, 11)
        # If the minimum sits at the outer edge the energy is still falling as the
        # molecules separate -- this orientation has no binding well in range, so
        # anchor the scan at the geometric (vdW) contact instead.
        if k >= n - 1:
            return seed
        return d_min

    @staticmethod
    def _minimize_on_grid(func, lo, hi, n):
        """Minimize a 1-D function on a uniform grid, refining parabolically.

        Coarse but robust and derivative-free -- enough to anchor the scan.
        Returns (d_min, k, n) where k is the index of the lowest grid point (so
        the caller can tell an interior minimum from one pinned at an edge).
        """
        ds = np.linspace(lo, hi, n)
        es = np.array([func(float(d)) for d in ds])
        k = int(np.argmin(es))
        d_min = float(ds[k])
        if 0 < k < n - 1:
            e0, e1, e2 = es[k - 1], es[k], es[k + 1]
            denom = e0 - 2.0 * e1 + e2
            if denom > 0.0:  # concave up -> a real interior minimum
                h = ds[1] - ds[0]
                d_min = float(ds[k] + 0.5 * h * (e0 - e2) / denom)
        return d_min, k, n

    # ----------------------------------------------------------------- #
    # Energy-stratified radial sampling ('spacing' = 'energy-stratified')
    # ----------------------------------------------------------------- #

    def _energy_profile(self, engine, assemble, seed, P):
        """Sample the interaction energy ΔE(R) along the approach axis.

        Energy-anchored: the van der Waals ``seed`` is only a starting guess to
        bracket the minimum search -- it does not set the scan bounds. The
        profile is built relative to the located energy minimum:

        * **outward** from the minimum to the maximum separation -- the
          attractive tail plus the far-point reference (E at max separation ≈
          E_A + E_B, so ΔE → 0 there and ΔE tracks only the intermolecular
          physics); and
        * **inward** from the minimum, one step at a time, until ΔE rises past
          the largest positive target level -- so the repulsive wall is covered
          *by energy*, not by a vdW-relative offset (down to a physical floor).

        Returns ``(ds, dE, d_min, De)``: the center-to-center grid ``ds`` (Å),
        the interaction energy ``dE`` (kJ/mol), the ΔE minimum ``d_min`` (Å), and
        the well depth ``De`` = −min(ΔE) (kJ/mol, 0 if there is no well). Note
        'innermost gap' is not used here -- it applies only to the vdW-contact /
        geometric spacings.
        """
        # Extend the outward scan to the anchor separation when tail coverage is
        # on, so the ≈0 far tail (and its far-point reference) reaches out to the
        # asymptote anchors. Same number of points -> no extra engine calls.
        max_sep = self._outer_separation(P)
        n = max(int(P["number of separations"]), 11)
        floor = max(seed - 2.0, 1.0)  # keep separations physically sane
        hartree_to_kJmol = Q_(1.0, "hartree").to("kJ/mol").magnitude
        kBT = self._kBT(P)

        def energy_at(d):
            engine.set_coordinates(assemble(float(d)), units="Å")
            return engine.energy(units="hartree")

        # 1. Locate the energy minimum; the vdW seed only brackets the search.
        anchor = self._energy_anchor(engine, assemble, seed, P)

        # 2. Outward branch: anchor -> outer separation (tail + far reference).
        d_hi = max_sep if max_sep > anchor + 0.5 else anchor + 3.0
        ds = list(np.linspace(anchor, d_hi, n))
        es = [energy_at(d) for d in ds]
        e_ref = es[-1]  # E at maximum separation ≈ E_A + E_B
        step = (d_hi - anchor) / (n - 1)

        # 3. Inward branch: walk in from the anchor until ΔE reaches the largest
        #    positive target level (the repulsive wall), down to the floor.
        d = anchor - step
        while d >= floor:
            e = energy_at(d)
            ds.insert(0, d)
            es.insert(0, e)
            dE = (np.array(es) - e_ref) * hartree_to_kJmol
            De = float(-dE.min()) if dE.min() < 0.0 else 0.0
            targets = self._energy_levels(De, kBT, P)
            max_pos = max([t for t in targets if t > 0.0], default=0.0)
            if (e - e_ref) * hartree_to_kJmol >= max_pos:
                break
            d -= step

        ds = np.array(ds)
        dE = (np.array(es) - e_ref) * hartree_to_kJmol
        k = int(np.argmin(dE))
        d_min = float(ds[k])
        De = float(-dE[k]) if dE[k] < 0.0 else 0.0
        return ds, dE, d_min, De

    def _kBT(self, P):
        """The thermal energy k_B·T (kJ/mol) for the sampling temperature."""
        T = P["sampling temperature"].to("K").magnitude
        return float((Q_(T, "K") * ureg.molar_gas_constant).to("kJ/mol").magnitude)

    def _energy_levels(self, De, kBT, P):
        """The target ΔE levels (kJ/mol) from the 'energy levels' spec.

        Each comma-separated token is a linear expression in the symbols ``De``
        (well depth) and ``kBT`` (thermal energy), both in kJ/mol -- e.g.
        ``-De``, ``-De/2``, ``0``, ``kBT``, ``5*kBT``. Only arithmetic on those
        two names is allowed (evaluated with no builtins and a whitelisted
        character set), since this is the user's own flowchart parameter.
        """
        import re

        namespace = {"De": float(De), "kBT": float(kBT), "__builtins__": {}}
        levels = []
        for token in P["energy levels"].split(","):
            token = token.strip()
            if token == "":
                continue
            if not re.fullmatch(r"[0-9DekBT.+\-*/() ]+", token):
                raise ValueError(
                    f"Invalid ΔE level '{token}': only numbers, the symbols 'De' "
                    "and 'kBT', and + - * / ( ) are allowed."
                )
            try:
                levels.append(float(eval(token, namespace)))  # noqa: S307
            except Exception as e:
                raise ValueError(f"Could not evaluate the ΔE level '{token}': {e}")
        return levels

    def _repulsive_cap(self, P):
        """The largest positive ΔE level (kJ/mol) -- the repulsive-wall cap.

        Points more repulsive than this are dropped, so the innermost gap can
        never push a configuration to an absurd energy (e.g. +130 kJ/mol). The
        positive levels are kBT-based, so De is irrelevant here.
        """
        levels = self._energy_levels(0.0, self._kBT(P), P)
        positive = [t for t in levels if t > 0.0]
        return max(positive) if positive else None

    def _outer_separation(self, P):
        """The outer scan bound: the anchor separation when tail coverage is on,
        else the maximum separation, so the ≈0 far tail is sampled out to where
        the asymptote anchors live."""
        max_sep = P["maximum separation"].to("Å").magnitude
        if self._truthy(P.get("tail coverage", "no")):
            return max(max_sep, P["anchor separation"].to("Å").magnitude)
        return max_sep

    def _energy_candidates(self, ds, dE, P):
        """Dense (distance, ΔE) candidate points from one orientation's profile.

        Interpolates ΔE on a fine radial grid spanning the sampled profile and
        drops any point above the repulsive cap. These candidates from every
        orientation are pooled and then globally stratified by energy; a dense
        grid (rather than the ~handful of profile points) lets the global
        binning fill each energy band smoothly.
        """
        interp = self._make_interpolator(ds, dE)
        n = max(int(P["number of separations"]) * 3, 40)
        grid = np.linspace(float(ds[0]), float(ds[-1]), n)
        cap = self._repulsive_cap(P)
        out = []
        for d in grid:
            e = interp(float(d))
            if cap is None or e <= cap:
                out.append((float(d), e))
        return out

    def _global_stratify(self, dE_values, P, rng):
        """Indices of a flat-in-energy subset of the pooled candidates.

        Bins the pooled interaction energies across [min, repulsive cap] into
        'number of energy bins' bins and keeps up to (target / n_bins) per bin,
        randomly subsampling over-full bins (sparser bins -- e.g. the deep well
        -- contribute all they have). This is what makes the ensemble flat in
        energy: per-orientation radial targeting cannot, because uniform random
        orientations pile up near ΔE = 0.
        """
        dE = np.asarray(dE_values, dtype=float)
        if dE.size == 0:
            return []
        n_bins = max(int(P["number of energy bins"]), 1)
        per_bin = max(int(P["target configurations"]) // n_bins, 1)
        lo = float(dE.min())
        cap = self._repulsive_cap(P)
        hi = float(cap) if (cap is not None and cap > lo) else float(dE.max())
        if hi <= lo:
            hi = lo + 1.0
        edges = np.linspace(lo, hi, n_bins + 1)
        which = np.clip(np.digitize(dE, edges) - 1, 0, n_bins - 1)
        selected = []
        for b in range(n_bins):
            members = np.where(which == b)[0]
            if len(members) <= per_bin:
                selected.extend(members.tolist())
            else:
                pick = rng.choice(members, size=per_bin, replace=False)
                selected.extend(int(i) for i in pick)
        return sorted(selected)

    def _candidate_dimers(
        self, candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx
    ):
        """Build ``dimer_analysis.Dimer`` objects for a list of candidates."""
        from dimer_builder_step import dimer_analysis

        out = []
        for o, d, _e in candidates:
            coords = orient_data[o]["rebuild"](d)
            out.append(
                dimer_analysis.Dimer(
                    symbols_A=symbols_A,
                    xyz_A=coords[a_idx],
                    symbols_B=symbols_B,
                    xyz_B=coords[b_idx],
                )
            )
        return out

    def _cluster_pick(self, dimers, dE, target, energy_weight):
        """DIRECT featurize → cluster → nearest-centroid pick; local indices.

        Featurizes each dimer by the collective variables ``dimer_analysis``
        computes (separation, approach unit vector, relative orientation,
        closest contact); if ``dE`` is given, the interaction energy is appended
        as an extra feature scaled by ``energy_weight`` (so energy can be made to
        count for more than one geometric dimension). Follows DIRECT (Qi et al.,
        npj Comput. Mater. 2024): standardize → PCA (variance-weighted scores) →
        BIRCH into ``target`` clusters → keep the member nearest each centroid.
        Returns the kept indices into ``dimers``.
        """
        import warnings

        from dimer_builder_step import dimer_analysis
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        from sklearn.cluster import Birch
        from sklearn.exceptions import ConvergenceWarning

        n = len(dimers)
        if n <= target:
            return list(range(n))

        m = dimer_analysis.compute_metrics(dimers)
        cols = [
            m.R,
            m.approach_vec[:, 0],
            m.approach_vec[:, 1],
            m.approach_vec[:, 2],
            m.orient_angle,
            m.min_contact,
        ]
        if dE is not None:
            cols.append(np.asarray(dE, dtype=float))
        feats = np.column_stack(cols)
        # Fill non-finite entries (e.g. orientation angle for a monatomic
        # fragment) with the column mean so every candidate is clusterable.
        for j in range(feats.shape[1]):
            col = feats[:, j]
            bad = ~np.isfinite(col)
            if bad.any():
                col[bad] = np.nanmean(col[~bad]) if (~bad).any() else 0.0

        X = StandardScaler().fit_transform(feats)
        if dE is not None:
            X[:, -1] *= float(energy_weight)  # up-weight the energy axis

        Z = PCA(n_components=min(X.shape)).fit_transform(X)
        # Normalize the (variance-weighted) PCA scores to an overall unit scale so
        # the BIRCH threshold behaves consistently, then shrink the threshold
        # until BIRCH resolves at least 'target' subclusters (it warns and
        # returns fewer otherwise -- the classic BIRCH pitfall).
        scale = float(Z.std())
        if scale > 0.0:
            Z = Z / scale
        threshold = 0.5
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            for _ in range(10):
                model = Birch(threshold=threshold, n_clusters=target).fit(Z)
                if len(model.subcluster_centers_) >= target:
                    break
                threshold *= 0.5
        labels = model.labels_

        keep = []
        for lab in np.unique(labels):
            members = np.where(labels == lab)[0]
            centroid = Z[members].mean(axis=0)
            nearest = members[np.argmin(((Z[members] - centroid) ** 2).sum(axis=1))]
            keep.append(int(nearest))
        return sorted(keep)

    def _direct_select(
        self, candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
    ):
        """Method A: one global DIRECT clustering over the CVs **plus ΔE**.

        ΔE is included as a feature scaled by the 'energy weight' parameter, so a
        single clustering covers geometry and energy together; near-duplicate
        geometries collapse (de-duplication) and the kept set spreads over both
        axes. Returns the kept indices into ``candidates``.
        """
        target = max(int(P["target configurations"]), 1)
        if len(candidates) <= target:
            return list(range(len(candidates)))
        dimers = self._candidate_dimers(
            candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx
        )
        dE = [c[2] for c in candidates]
        weight = float(P.get("energy weight", 1.0))
        return self._cluster_pick(dimers, dE, target, weight)

    def _select_within_bins(
        self, candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
    ):
        """Method B: flat-in-energy bins, DIRECT-diversify the geometry per bin.

        Bins the candidates by ΔE (guaranteeing flat energy, like 'energy bins'),
        then within each over-full bin keeps a geometrically diverse subset via
        DIRECT on the geometric CVs alone (ΔE excluded -- energy is already fixed
        by the bin). Returns the kept indices into ``candidates``.
        """
        dE = np.array([c[2] for c in candidates], dtype=float)
        if dE.size == 0:
            return []
        n_bins = max(int(P["number of energy bins"]), 1)
        per_bin = max(int(P["target configurations"]) // n_bins, 1)
        lo = float(dE.min())
        cap = self._repulsive_cap(P)
        hi = float(cap) if (cap is not None and cap > lo) else float(dE.max())
        if hi <= lo:
            hi = lo + 1.0
        edges = np.linspace(lo, hi, n_bins + 1)
        which = np.clip(np.digitize(dE, edges) - 1, 0, n_bins - 1)

        keep = []
        for b in range(n_bins):
            members = np.where(which == b)[0]
            if len(members) <= per_bin:
                keep.extend(int(i) for i in members)
                continue
            subset = [candidates[i] for i in members]
            dimers = self._candidate_dimers(
                subset, orient_data, symbols_A, symbols_B, a_idx, b_idx
            )
            local = self._cluster_pick(dimers, None, per_bin, 1.0)
            keep.extend(int(members[j]) for j in local)
        return sorted(keep)

    def _select_pool(
        self, candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P, rng
    ):
        """Down-select the pooled candidates by the chosen 'selection method'.

        Returns the kept subset of ``candidates`` (both methods target about
        'target configurations').
        """
        method = P.get("selection method", "energy bins")
        if method == "descriptor diversity":
            keep = self._direct_select(
                candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
            )
        elif method == "energy bins + diversity":
            keep = self._select_within_bins(
                candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
            )
        else:
            keep = self._global_stratify([c[2] for c in candidates], P, rng)

        if self._truthy(P.get("tail coverage", "no")):
            keep = self._add_distance_coverage(
                candidates, keep, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
            )
        return [candidates[i] for i in sorted(set(keep))]

    def _add_distance_coverage(
        self, candidates, keep, orient_data, symbols_A, symbols_B, a_idx, b_idx, P
    ):
        """Supplement the energy-selected core with long-range distance coverage.

        Energy-flat selection starves the weak tail (5+ Å ≈ 0 ΔE), so this adds,
        independent of interaction strength: (1) a **distance-coverage floor** --
        at least 'tail configurations per bin' in every 'tail spacing' bin from
        'tail minimum separation' to the maximum separation; and (2) a few
        **asymptote anchors** between the maximum and 'anchor separation' to pin
        E → 0 at large distance. Extra configurations are chosen for geometric
        diversity (over orientations). Returns the augmented set of indices.
        """
        d = np.array([c[1] for c in candidates], dtype=float)
        keep = set(keep)

        def fill(pool_idx, count):
            """Add up to `count` geometrically-diverse candidates from pool_idx."""
            avail = [i for i in pool_idx if i not in keep]
            if not avail or count <= 0:
                return
            if len(avail) <= count:
                keep.update(avail)
                return
            sub = [candidates[i] for i in avail]
            dimers = self._candidate_dimers(
                sub, orient_data, symbols_A, symbols_B, a_idx, b_idx
            )
            local = self._cluster_pick(dimers, None, count, 1.0)
            keep.update(avail[j] for j in local)

        # (1) distance-coverage floor over [tail_min, maximum separation]
        tail_min = P["tail minimum separation"].to("Å").magnitude
        near_max = P["maximum separation"].to("Å").magnitude
        spacing = max(P["tail spacing"].to("Å").magnitude, 1.0e-3)
        per_bin = max(int(P["tail configurations per bin"]), 0)
        edge = tail_min
        while per_bin > 0 and edge < near_max - 1.0e-6:
            hi = min(edge + spacing, near_max)
            members = np.where((d >= edge) & (d < hi))[0]
            have = sum(1 for i in members if i in keep)
            fill(members, per_bin - have)
            edge = hi

        # (2) asymptote anchors over (maximum separation, anchor separation]
        n_anchors = max(int(P["asymptote anchors"]), 0)
        if n_anchors > 0:
            anchor_max = P["anchor separation"].to("Å").magnitude
            members = np.where((d > near_max) & (d <= anchor_max))[0]
            already = sum(1 for i in members if i in keep)
            fill(members, n_anchors - already)

        return sorted(keep)

    @staticmethod
    def _make_interpolator(ds, dE):
        """A function mapping a distance (Å) to its interpolated ΔE (kJ/mol)."""

        def dE_at(d):
            return float(np.interp(d, ds, dE))

        return dE_at

    def _accept_orientation(self, De, P, rng):
        """Whether to keep an orientation given its ΔE well depth (kJ/mol).

        Applies the 'orientation weighting' rule; a purely repulsive/glancing
        orientation (De == 0) is always shallow. 'downweight by depth' needs the
        random generator; without one it degrades to the reject rule.
        """
        mode = P["orientation weighting"]
        min_depth = P["minimum well depth"].to("kJ/mol").magnitude
        if mode == "none":
            return True
        if mode == "downweight by depth" and rng is not None and min_depth > 0.0:
            keep_probability = 1.0 - 2.0 ** (-De / min_depth)
            return bool(rng.random() < keep_probability)
        return De >= min_depth  # 'reject shallow orientations' (and the fallback)

    def _plan_scan(self, engine, assemble, contact, P, rng=None, weight=False):
        """Plan one orientation's radial scan.

        Returns ``(distances, gap_reference, dE_at, De)`` -- the center-to-center
        distances to place, the reference distance for the reported 'gap' (the
        ΔE minimum for energy-stratified/energy scans, else the vdW contact), a
        function giving ΔE at a distance (or ``None``), and the well depth (or
        ``None``). Returns ``None`` when the orientation is rejected by the
        weighting rule.
        """
        if engine is not None and P["spacing"] == "energy-stratified":
            ds, dE, d_min, De = self._energy_profile(engine, assemble, contact, P)
            if weight and not self._accept_orientation(De, P, rng):
                return None
            distances = self._stratified_separations(ds, dE, d_min, De, P)
            return distances, d_min, self._make_interpolator(ds, dE), De

        if engine is not None:
            # Anchor the ladder at the energy minimum, then record the
            # interaction energy at each scan point so the diagnostics carry ΔE
            # (and the interaction-energy property is stored) even though the
            # spacing itself is geometric/linear/explicit.
            contact = self._energy_anchor(engine, assemble, contact, P)
            distances = self._separation_schedule(contact, P)
            dE_at, De = self._interaction_energies(engine, assemble, distances, P)
            return distances, contact, dE_at, De
        return self._separation_schedule(contact, P), contact, None, None

    def _interaction_energies(self, engine, assemble, distances, P):
        """ΔE (kJ/mol) at each scan distance, referenced to the far-point asymptote.

        Used by the non-stratified energy paths so every energy-contact run
        carries per-configuration interaction energies (the stratified path
        already gets them from its ΔE(R) profile). The reference is the energy at
        the maximum separation (≈ E_A + E_B). Returns ``(dE_at, De)``.
        """
        max_sep = P["maximum separation"].to("Å").magnitude
        hartree_to_kJmol = Q_(1.0, "hartree").to("kJ/mol").magnitude

        def energy_at(d):
            engine.set_coordinates(assemble(float(d)), units="Å")
            return engine.energy(units="hartree")

        ds = np.unique(np.asarray(distances, dtype=float))
        e_ref = energy_at(max(max_sep, float(ds[-1])))
        es = np.array([energy_at(d) for d in ds])
        dE = (es - e_ref) * hartree_to_kJmol
        De = float(-dE.min()) if dE.min() < 0.0 else 0.0
        return self._make_interpolator(ds, dE), De

    def _emit_configurations(
        self,
        out_sys,
        base,
        candidates,
        orient_data,
        P,
        *,
        save_props,
        collect,
        templates,
        fixed_ids,
        movable_ids,
        symbols_A,
        symbols_B,
        a_idx,
        b_idx,
    ):
        """Build the selected candidates into configurations.

        ``candidates`` is a list of ``(orientation_index, distance, ΔE-or-None)``
        (already selected/ordered); ``orient_data[i]`` holds that orientation's
        ``rebuild(d) -> coords`` closure, gap reference, well depth, and angles.
        Shared by both input modes. Returns ``(separations, ensemble, count)``.
        """
        t_fixed, t_movable = templates
        count = 0
        separations = []
        ensemble = []
        point_of = {}
        for orient_index, d, dE_value in candidates:
            o = orient_data[orient_index]
            coords = o["rebuild"](d)
            point_of[orient_index] = point_of.get(orient_index, 0) + 1
            geometry = {
                "separation": d,
                "gap": d - o["gap_ref"],
                "orientation": o["orientation"],
                "theta": o["theta"],
                "phi": o["phi"],
                "alpha": o["alpha"],
                "beta": o["beta"],
                "gamma": o["gamma"],
            }
            if dE_value is not None:
                geometry["interaction_energy"] = dE_value
                geometry["well_depth"] = o["De"]
            conf = base if count == 0 else out_sys.copy_configuration(base)
            conf.atoms.set_coordinates(coords, fractionals=False)
            self._name_and_tag(
                conf,
                P,
                count + 1,
                o["orientation"],
                point_of[orient_index],
                geometry,
                save_props,
            )
            self._add_subsets(conf, t_fixed, t_movable, fixed_ids, movable_ids)
            if collect:
                ensemble.append(
                    self._make_dimer_record(
                        symbols_A,
                        coords[a_idx],
                        symbols_B,
                        coords[b_idx],
                        geometry,
                        conf.name,
                    )
                )
            separations.append(d)
            count += 1
        return separations, ensemble, count

    def _collect_candidates(
        self,
        engine,
        assemble,
        contact,
        P,
        rng,
        orient_meta,
        orient_data,
        candidates,
        weight=True,
    ):
        """Plan one orientation and append its candidate points to the pool.

        For energy-stratified spacing, computes the ΔE(R) profile and emits a
        dense, repulsive-capped candidate set (globally stratified later); for
        the other spacings, uses the per-orientation radial schedule directly.
        ``weight`` enables the orientation pre-filter (off for prepared dimers,
        which are always kept). Returns True if the orientation contributed.
        """
        if engine is not None and P["spacing"] == "energy-stratified":
            ds, dE, d_min, De = self._energy_profile(engine, assemble, contact, P)
            if weight and not self._accept_orientation(De, P, rng):
                return False
            dE_at = self._make_interpolator(ds, dE)
            cand = self._energy_candidates(ds, dE, P)
            gap_ref = d_min
        else:
            plan = self._plan_scan(engine, assemble, contact, P, rng=rng, weight=weight)
            if plan is None:
                return False
            distances, gap_ref, dE_at, De = plan
            cand = [(float(d), dE_at(float(d)) if dE_at else None) for d in distances]

        o_idx = len(orient_data)
        orient_data.append({**orient_meta, "gap_ref": gap_ref, "De": De})
        for d, e in cand:
            candidates.append((o_idx, d, e))
        return True

    def _build(self, system_db, P, rng):
        """Generate the dimer configurations. Returns (system, stats)."""
        # Backstop for hand-edited/scripted flowcharts (the GUI prevents this):
        # energy-stratified spacing has no meaning without an energy engine.
        if P["spacing"] == "energy-stratified" and P["contact method"] != "energy":
            raise ValueError(
                "'energy-stratified' spacing requires the 'energy' contact method "
                "(an MDI engine supplied by a Model Chemistry step)."
            )
        if P["input mode"] == "two monomer sets":
            return self._build_from_monomers(system_db, P, rng)
        else:
            return self._build_from_dimers(system_db, P, rng)

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
        symbols_A = list(A0.atoms.symbols)
        symbols_B = list(B0.atoms.symbols)

        engine = None
        if P["contact method"] == "energy":
            elements = list(A0.atoms.atomic_numbers) + list(B0.atoms.atomic_numbers)
            charge = (A0.charge or 0) + (B0.charge or 0)
            engine = self._open_energy_engine(elements, charge, 1)

        save_props = self._truthy(P["save scan variables as properties"])
        collect = self._collect_ensemble(P)
        global_strat = engine is not None and P["spacing"] == "energy-stratified"

        # Phase 1: for each orientation, pool its candidate scan points.
        orient_data = []
        candidates = []
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

                def assemble(d, xyzA=xyzA, xyzB=xyzB, axis=axis):
                    return np.vstack([xyzA, xyzB + axis * d])

                self._collect_candidates(
                    engine,
                    assemble,
                    contact,
                    P,
                    rng,
                    {
                        "rebuild": assemble,
                        "orientation": orientation,
                        "theta": theta,
                        "phi": phi,
                        "alpha": alpha,
                        "beta": beta,
                        "gamma": gamma,
                    },
                    orient_data,
                    candidates,
                )
        finally:
            if engine is not None:
                engine.close()

        # Phase 2: globally down-select the pooled candidates.
        a_idx = np.arange(nA)
        b_idx = np.arange(nA, nA + B0.n_atoms)
        if global_strat:
            candidates = self._select_pool(
                candidates, orient_data, symbols_A, symbols_B, a_idx, b_idx, P, rng
            )
        candidates.sort(key=lambda c: (c[0], c[1]))

        # Phase 3: build the selected candidates.
        separations, ensemble, count = self._emit_configurations(
            dimer_sys,
            base,
            candidates,
            orient_data,
            P,
            save_props=save_props,
            collect=collect,
            templates=(t_fixed, t_movable),
            fixed_ids=fixed_ids,
            movable_ids=movable_ids,
            symbols_A=symbols_A,
            symbols_B=symbols_B,
            a_idx=a_idx,
            b_idx=b_idx,
        )

        stats = self._stats(name, P["number of orientations"], separations)
        stats["ensemble"] = ensemble
        if engine is not None:
            stats["model_chemistry"] = self._energy_model
            stats["n_energy_calls"] = engine.n_energy_calls
        return dimer_sys, stats

    def _build_from_dimers(self, system_db, P, rng):
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
        fixed_masses = masses[fixed_idx]
        movable_masses = masses[movable_idx]

        t_fixed, t_movable = self._ensure_templates(system_db)
        atom_ids = base.atoms.ids
        fixed_ids = [atom_ids[i] for i in fixed_idx]
        movable_ids = [atom_ids[i] for i in movable_idx]
        all_symbols = list(d0.atoms.symbols)
        symbols_fixed = [all_symbols[i] for i in fixed_idx]
        symbols_movable = [all_symbols[i] for i in movable_idx]

        engine = None
        if P["contact method"] == "energy":
            engine = self._open_energy_engine(
                list(d0.atoms.atomic_numbers), d0.charge or 0, 1
            )

        save_props = self._truthy(P["save scan variables as properties"])
        collect = self._collect_ensemble(P)
        global_strat = engine is not None and P["spacing"] == "energy-stratified"

        # Phase 1: pool candidate scan points from each prepared dimer.
        orient_data = []
        candidates = []
        fixed_idx = np.asarray(fixed_idx)
        movable_idx = np.asarray(movable_idx)
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
                fixed_center = np.average(fixed_xyz, axis=0, weights=fixed_masses)
                movable_center = np.average(movable_xyz, axis=0, weights=movable_masses)
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

                # Prepared dimers are user-provided orientations, so they are
                # always kept (weight=False); the pool is still energy-stratified.
                self._collect_candidates(
                    engine,
                    assemble,
                    contact,
                    P,
                    rng,
                    {
                        "rebuild": assemble,
                        "orientation": orientation,
                        "theta": theta,
                        "phi": phi,
                        "alpha": alpha,
                        "beta": beta,
                        "gamma": gamma,
                    },
                    orient_data,
                    candidates,
                    weight=False,
                )
        finally:
            if engine is not None:
                engine.close()

        # Phase 2: globally down-select the pooled candidates.
        if global_strat:
            candidates = self._select_pool(
                candidates,
                orient_data,
                symbols_fixed,
                symbols_movable,
                fixed_idx,
                movable_idx,
                P,
                rng,
            )
        candidates.sort(key=lambda c: (c[0], c[1]))

        # Phase 3: build the selected candidates.
        separations, ensemble, count = self._emit_configurations(
            out_sys,
            base,
            candidates,
            orient_data,
            P,
            save_props=save_props,
            collect=collect,
            templates=(t_fixed, t_movable),
            fixed_ids=fixed_ids,
            movable_ids=movable_ids,
            symbols_A=symbols_fixed,
            symbols_B=symbols_movable,
            a_idx=fixed_idx,
            b_idx=movable_idx,
        )

        stats = self._stats(name, len(orient_data), separations)
        stats["ensemble"] = ensemble
        if engine is not None:
            stats["model_chemistry"] = self._energy_model
            stats["n_energy_calls"] = engine.n_energy_calls
        return out_sys, stats

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
        else:  # orientation,distance
            conf.name = f"{orientation},{point}"

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
            if "interaction_energy" in geometry:
                self._put_property(
                    conf,
                    f"interaction energy{p}",
                    geometry["interaction_energy"],
                    "kJ/mol",
                )
                self._put_property(
                    conf, f"well depth{p}", geometry["well_depth"], "kJ/mol"
                )

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
