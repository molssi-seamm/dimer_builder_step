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

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in dimer_builder_parameters.py

    See Also
    --------
    DimerBuilder, TkDimerBuilder,
    DimerBuilderParameters,
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
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
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
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the dimer_builder_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkDimerBuilder.reset_dialog
        """

        super().create_dialog(title="Dimer Builder", widget="notebook")

        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a tab for the sub-flowchart (used by the 'energy' contact method)
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = flowchart_frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()

        # Shortcut for parameters
        P = self.node.parameters

        # A frame to hold the control parameters
        parameters_frame = self["parameters frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Dimer Builder Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in DimerBuilderParameters.parameters:
            if key not in ("results",):
                self[key] = P[key].widget(parameters_frame)

        # Comboboxes whose value changes the layout re-lay out the dialog
        for key in (
            "input mode",
            "spacing",
            "monomer A configurations",
            "monomer B configurations",
        ):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].combobox.bind("<Return>", self.reset_dialog)
            self[key].combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog as needed for the current state.

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None
        """
        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        self["parameters frame"].grid(row=0, column=0, sticky=tk.EW, pady=10)
        frame.columnconfigure(0, weight=1)

        self.reset_parameters_frame()

        return 1

    def reset_parameters_frame(self):
        """Lay out the control parameters according to the current choices."""
        mode = self["input mode"].get()
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

        # Input sources
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

        # Orientation sampling only applies when assembling from monomers
        if mode == "two monomer sets":
            add("number of orientations")
            add("random seed")

        # Radial scan
        for key in ("contact method", "innermost gap", "maximum separation", "spacing"):
            add(key)
        if spacing == "explicit":
            add("separations")
        else:
            add("number of separations")

        # Output
        for key in (
            "system name",
            "configuration name",
            "save scan variables as properties",
        ):
            add(key)

        sw.align_labels(widgets, sticky=tk.E)
        frame.columnconfigure(1, weight=1)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkDimerBuilder.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
