# -*- coding: utf-8 -*-

"""The graphical part of a Dimer Builder step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

from .dimer_builder_parameters import DimerBuilderParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkDimerBuilder(seamm.TkNode):
    """
    The graphical part of a Dimer Builder step in a flowchart.

    See Also
    --------
    DimerBuilder, DimerBuilderParameters
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """Initialize a graphical node."""
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """Create the dialog for the Dimer Builder parameters.

        The base class builds a notebook with 'Parameters' and 'Results' tabs
        (the step has a 'results' parameter). We deliberately do NOT add a
        'Flowchart' tab: the sub-flowchart is unused now that the 'energy'
        contact method drives an engine through a model chemistry over MDI. If a
        sub-flowchart-based fallback is added later, the tab can be reinstated
        (behind an advanced option).
        """
        super().create_dialog(title="Dimer Builder", widget="notebook")

        # Shortcut for parameters
        P = self.node.parameters

        frame = self["parameters frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Dimer Builder Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in DimerBuilderParameters.parameters:
            if key not in ("results",):
                self[key] = P[key].widget(frame)

        # Shown only when the energy contact method is selected.
        self["energy note"] = ttk.Label(
            frame,
            text=(
                "The 'energy' contact method needs a Model Chemistry step before "
                "this step, to define the engine and method (e.g. MOPAC PM6-ORG)."
            ),
            foreground="blue",
            wraplength=500,
            justify=tk.LEFT,
        )

        # Comboboxes whose value changes the layout re-lay out the dialog.
        for key in (
            "input mode",
            "spacing",
            "contact method",
            "selection method",
            "tail coverage",
            "orientation weighting",
            "monomer A configurations",
            "monomer B configurations",
        ):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].combobox.bind("<Return>", self.reset_dialog)
            self[key].combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the parameters frame in the Parameters tab."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        self["parameters frame"].grid(row=0, column=0, sticky=tk.EW, pady=10)
        frame.columnconfigure(0, weight=1)

        self.reset_parameters_frame()
        return 1

    def reset_parameters_frame(self):
        """Lay out the control parameters for the current choices."""
        mode = self["input mode"].get()
        energy = self["contact method"].get() == "energy"

        # Energy-stratified spacing is only meaningful with an energy engine, so
        # narrow the 'spacing' choices to exclude it otherwise (prevent the
        # invalid combination rather than only catching it at run time).
        spacings = ["geometric", "linear", "explicit"]
        if energy:
            spacings.append("energy-stratified")
        self["spacing"].combobox.config(values=spacings)
        if self["spacing"].get() not in spacings:
            self["spacing"].set("geometric")
        spacing = self["spacing"].get()

        frame = self["parameters frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        widgets = []

        def add(key):
            nonlocal row
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        add("input mode")

        if mode == "two monomer sets":
            sources = ("monomer A", "monomer B")
        else:
            sources = ("monomer A",)
        for prefix in sources:
            add(prefix)
            add(f"{prefix} configurations")
            if self[f"{prefix} configurations"].get() in (
                "name is",
                "name matches",
                "name regexp",
            ):
                add(f"{prefix} configuration name")

        if mode == "two monomer sets":
            add("number of orientations")
            add("random seed")

        add("contact method")
        if energy and not self._upstream_has_model_chemistry():
            self["energy note"].grid(
                row=row, column=0, columnspan=2, sticky=tk.W, pady=(0, 6)
            )
            row += 1
        for key in ("innermost gap", "maximum separation", "spacing"):
            # Energy-stratified spacing anchors on the energy minimum and sets
            # its inner bound by energy, so the innermost gap does not apply.
            if key == "innermost gap" and spacing == "energy-stratified":
                continue
            add(key)
        if spacing == "explicit":
            add("separations")
        elif spacing == "energy-stratified":
            add("number of separations")  # ΔE(R) profile resolution
            add("energy levels")
            add("sampling temperature")
            add("selection method")
            add("target configurations")
            selection = self["selection method"].get()
            if selection in ("energy bins", "energy bins + diversity"):
                add("number of energy bins")
            if selection == "descriptor diversity":
                add("energy weight")
            if mode == "two monomer sets":
                add("orientation weighting")
                if self["orientation weighting"].get() != "none":
                    add("minimum well depth")
            add("tail coverage")
            if self["tail coverage"].get() == "yes":
                add("tail minimum separation")
                add("tail spacing")
                add("tail configurations per bin")
                add("asymptote anchors")
                add("anchor separation")
        else:
            add("number of separations")

        for key in (
            "system name",
            "configuration name",
            "save scan variables as properties",
            "analysis plots",
        ):
            add(key)

        sw.align_labels(widgets, sticky=tk.E)
        frame.columnconfigure(1, weight=1)

    def _upstream_has_model_chemistry(self):
        """True if a Model Chemistry step precedes this one in the flowchart.

        Uses the shared ``previous_nodes()`` helper and checks the Python type by
        name + module (no import dependency on model_chemistry_step). On any error
        (e.g. the node is not yet linked into the flowchart) returns False, so the
        reminder is shown -- the safe default.
        """
        try:
            return any(
                type(node).__name__ == "ModelChemistry"
                and type(node).__module__.startswith("model_chemistry_step")
                for node in self.previous_nodes()
            )
        except Exception:
            return False

    def right_click(self, event):
        """Handle a right-click: add the Edit... item."""
        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)
        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
