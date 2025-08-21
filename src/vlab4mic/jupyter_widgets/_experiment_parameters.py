"""
_experiment_parameters
----------------------

This module provides EZInput-based widget interfaces for selecting and configuring experiment parameters
in the virtual microscopy simulation workflow. It enables interactive selection of structures, probes, sample parameters,
imaging modalities, and running experiments within Jupyter notebooks.

Functions
---------

- ui_select_structure(experiment):
    Returns a widget for selecting the experiment structure.

- update_widgets_visibility(ezwidget, visibility_dictionary):
    Utility to show/hide widgets in an EZInput widget based on a visibility dictionary.

- ui_select_probe(experiment, **kwargs):
    Returns a widget for selecting and adding probes to the experiment.

- ui_select_sample_parameters(experiment):
    Returns a widget for configuring sample parameters (e.g., number of particles, random orientations).

- ui_select_modality(experiment):
    Returns a widget for selecting and previewing imaging modalities.

- ui_run_experiment(experiment):
    Returns a widget for running the experiment and saving results.


Each function returns an EZInput-based widget for use in a Jupyter notebook.
"""

from ezinput import EZInput
import matplotlib.pyplot as plt
import copy
from IPython.display import display, clear_output
import ipywidgets as widgets
import io
from ipyfilechooser import FileChooser
import copy
from vlab4mic.utils.visualisation.matplotlib_plots import slider_normalised
import numpy as np
import tifffile as tif
from IPython.utils import io
from ..experiments import build_virtual_microscope

select_colour = "#4daf4ac7"
remove_colour = "#ff8000da"
update_colour = "#00bfffda"

select_icon = "fa-check"
add_icon = "fa-plus"
remove_icon = "fa-minus"
clear_icon = "fa-trash"
loding_icon = "fa-spinner fa-spin"
update_icon = "fa-wrench"  # create
toggle_icon = "fa-eye-slash"
upload_icon = "fa-upload"


def ui_select_structure(experiment):
    """
    Create a widget for selecting the experiment structure.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing available structures.

    Returns
    -------
    EZInput
        Widget for structure selection.
    """
    gui = EZInput("Select_structure")

    def select_structure(elements):
        elements["label_1"].value = "Current structure selected: Loading..."
        elements["select_structure"].icon = loding_icon
        # elements["select_structure"].disabled = True
        for wgt in elements.keys():
            elements[wgt].disabled = True
        experiment.structure_id = experiment.structures_info_list[
            elements["structures"].value
        ]
        with io.capture_output() as captured:
            experiment.build(modules="structure")
        update_structure_list()
        # elements["select_structure"].disabled = False
        elements["select_structure"].icon = select_icon
        for wgt in elements.keys():
            elements[wgt].disabled = False

    def select_structure_from_file(b):
        gui["label_1"].value = "Current structure selected: Loading..."
        gui["select_structure"].icon = loding_icon
        for wgt_name in gui.elements.keys():
            gui[wgt_name].disabled = True
        filepath = gui["File"].selected
        filename = gui["File"].selected_filename
        experiment.select_structure(
            structure_id=filename.split(".")[0], 
            structure_path=filepath
        )
        update_structure_list()
        gui["select_structure"].icon = select_icon
        for wgt_name in gui.elements.keys():
            gui[wgt_name].disabled = False


    def update_structure_list():
        if experiment.structure_id is not None:
            gui["label_1"].value = (
                "Current structure selected: " + experiment.structure_id
            )
        else:
            gui["label_1"].value = (
                "Current structure selected: " + "No structure selected yet."
            )

    def toggle_advanced_parameters(b):
        widgets_visibility["advanced_param_header"] = not widgets_visibility[
            "advanced_param_header"
        ]
        widgets_visibility["File"] = not widgets_visibility["File"]
        widgets_visibility["select_structure_from_file"] = not widgets_visibility[
            "select_structure_from_file"
        ]
        update_widgets_visibility(gui, widgets_visibility)

    if experiment.structure_id is not None:
        gui.add_label(
            value="Current structure selected: " + experiment.structure_id
        )
    else:
        gui.add_label(
            value="Current structure selected: No structure selected yet"
        )
    gui.add_dropdown(
        "structures",
        description="Select Structure:",
        options=experiment.structures_info_list.keys(),
    )
    gui.add_label(
        "Note: Time for structure loading varies depending on the size of the structure"
    )
    gui.add_button(
        "toggle_advanced_parameters",
        description="Toggle advanced parameters",
        icon=toggle_icon,
    )
    # advanced parameters
    gui.add_HTML(
        "advanced_param_header",
        "<b>Upload a PDB/CIF file</b>",
        style=dict(font_size="15px"),
    )
    gui.add_file_upload(
        "File",
        description="Select from file",
        accept=["*.pdb", "*.cif"],
        save_settings=False,
    )
    gui.add_button(
        "select_structure_from_file",
        description="Select structure from file",
        icon=select_icon,
        style={"button_color": select_colour},
    )        
    gui.add_callback(
        "select_structure",
        select_structure,
        gui.elements,
        description="Select structure",
        icon=select_icon,
        style={"button_color": select_colour},
    )
    gui["toggle_advanced_parameters"].on_click(toggle_advanced_parameters)
    gui["select_structure_from_file"].on_click(select_structure_from_file)
    widgets_visibility = {}
    _unstyle_widgets(gui, widgets_visibility)
    toggle_advanced_parameters(True)
    return gui


def update_widgets_visibility(ezwidget, visibility_dictionary):
    """
    Show or hide widgets in an EZInput widget based on a visibility dictionary.

    Parameters
    ----------
    ezwidget : EZInput
        The EZInput widget containing elements to show/hide.
    visibility_dictionary : dict
        Dictionary mapping widget names to booleans (True to show, False to hide).

    Returns
    -------
    None
    """
    for widgetname in visibility_dictionary.keys():
        if visibility_dictionary[widgetname]:
            ezwidget[widgetname].layout.display = "inline-flex"
        else:
            ezwidget[widgetname].layout.display = "None"


def _unstyle_widgets(ezwidget, visibility_dictionary):
    for wgt in ezwidget.elements.keys():
        visibility_dictionary[wgt] = True
        if isinstance(ezwidget[wgt], widgets.Button):
            ezwidget.elements[wgt].layout = widgets.Layout(
                width="50%",
                display="inline-flex",
                align_items="center",
                justify_content="center",
            )
        else:
            ezwidget.elements[wgt].layout = widgets.Layout(
                width="50%", display="inline-flex"
            )


def ui_select_probe(experiment, **kwargs):
    """
    Create a widget for selecting and adding probes to the experiment.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing probe configuration.

    Returns
    -------
    EZInput
        Widget for probe selection and addition.
    """
    experiment.remove_probes()
    probes_gui = EZInput(title="Labels")
    visibility_widgets = dict()
    probe_options = []
    if (
        experiment.structure_id
        in experiment.config_probe_per_structure_names.keys()
    ):
        probe_list = experiment.config_probe_per_structure_names[
            experiment.structure_id
        ]
        probe_options.extend(copy.copy(probe_list))
    # add probes with no targets
    for probe_name in experiment.config_global_probes_names:
        if experiment.config_probe_params[probe_name]["target"]["type"]:
            probe_options.append(probe_name)
    for probe_name in experiment.config_probe_models_names:
        probe_options.append(probe_name)

    # methods
    def select_probe(values):
        experiment.add_probe(
            probe_template=values["select_probe_template"].value,
        )

        # Set defect parameters when using simple probe selection
        defect_fraction = probes_gui["defect_fraction"].value
        defect_small = probes_gui["defect_small_cluster"].value
        defect_large = probes_gui["defect_large_cluster"].value

        if defect_fraction > 0 and defect_small > 0 and defect_large > 0:
            experiment.defect_eps["defect"] = float(defect_fraction)
            experiment.defect_eps["eps1"] = float(defect_small)
            experiment.defect_eps["eps2"] = float(defect_large)
            experiment.defect_eps["use_defects"] = True
        else:
            experiment.defect_eps["defect"] = 0.0
            experiment.defect_eps["eps1"] = 100.0
            experiment.defect_eps["eps2"] = 200.0
            experiment.defec_eps["use_defects"] = False

        probes_gui["create_particle"].disabled = False
        update_probe_list()

    def select_custom_probe(b):
        probe_template = probes_gui["select_probe_template"].value
        probe_name = probes_gui["probe_name"].value
        labelling_efficiency = probes_gui["labelling_efficiency"].value
        probe_distance_to_epitope = probes_gui["distance_from_epitope"].value
        probe_target_type = options_dictionary[probes_gui["mock_type"].value]
        probe_target_value = probes_gui["mock_type_options1"].value
        probe_target_value2 = probes_gui["mock_type_options2"].value
        as_linker = probes_gui["as_linker"].value
        if probes_gui["wobble"].value:
            probe_wobble_theta = probes_gui["wobble_theta"].value
        else:
            probe_wobble_theta = None

        # Handle defect parameters
        defect_fraction = probes_gui["defect_fraction"].value
        defect_small_cluster = probes_gui["defect_small_cluster"].value
        defect_large_cluster = probes_gui["defect_large_cluster"].value

        # Set defect parameters in experiment if all are provided and non-zero
        if (
            defect_fraction > 0
            and defect_small_cluster > 0
            and defect_large_cluster > 0
        ):
            experiment.defect_eps["defect"] = float(defect_fraction)
            experiment.defect_eps["eps1"] = float(defect_small_cluster)
            experiment.defect_eps["eps2"] = float(defect_large_cluster)
            experiment.defect_eps["use_defects"] = True
        else:
            experiment.defect_eps["defect"] = 0.0
            experiment.defect_eps["eps1"] = 100.0
            experiment.defect_eps["eps2"] = 200.0
            experiment.defect_eps["use_defects"] = False

        if as_linker:
            options_per_type1["Primary_Probe"] = [
                probe_name,
            ]
        if probe_target_type == "Sequence":
            experiment.add_probe(
                probe_template=probe_template,
                probe_name=probe_name,
                probe_target_type=probe_target_type,
                peptide_motif={
                    "chain_name": probe_target_value,
                    "position": probe_target_value2,
                },
                labelling_efficiency=labelling_efficiency,
                probe_distance_to_epitope=probe_distance_to_epitope,
                as_primary=as_linker,
                probe_wobble_theta=probe_wobble_theta,
            )
        elif probe_target_type == "Atom_residue":
            residue = probes_gui["mock_type_options1"].value
            atom = probes_gui["mock_type_options2"].value
            experiment.add_probe(
                probe_template=probe_template,
                probe_name=probe_name,
                probe_target_type=probe_target_type,
                probe_target_value=dict(atoms=atom, residues=residue),
                labelling_efficiency=labelling_efficiency,
                probe_distance_to_epitope=probe_distance_to_epitope,
                as_primary=as_linker,
                probe_wobble_theta=probe_wobble_theta,
            )
        elif probe_target_type == "Primary":
            experiment.add_probe(
                probe_template=probe_template,
                probe_name=probe_name,
                probe_target_type=probe_target_type,
                probe_target_value=probe_target_value,
                labelling_efficiency=labelling_efficiency,
                probe_distance_to_epitope=probe_distance_to_epitope,
                as_primary=as_linker,
                probe_wobble_theta=probe_wobble_theta,
            )
        probes_gui["create_particle"].disabled = False
        update_probe_list()

    def update_probe_list():
        probes_gui["message1"].value = ""
        for probe in experiment.probe_parameters.keys():
            probes_gui["message1"].value += probe + "<br>"

    def create_particle(b):
        probes_gui["message2"].value = "Creating labelled structure..."
        with io.capture_output() as captured:
            experiment.build(modules=["particle"])
        probes_gui["add_probe"].disabled = True
        probes_gui["create_particle"].disabled = True
        if experiment.generators_status("particle"):
            probes_gui["message2"].value = (
                "Labelled structure created successfully!"
            )
        else:
            probes_gui["message2"].value = (
                "Labelled structure creation failed. Check the logs for details."
            )

    def show_probe_info(change):
        probe_template = probes_gui["select_probe_template"].value
        probes_gui["probe_name"].value = probe_template
        if probe_template in experiment.config_probe_params.keys():
            info_text = "<b>Target: </b>"
            probe_info = experiment.config_probe_params[probe_template]
            if probe_info["target"]["type"] == "Atom_residue":
                target_type = "residue"
                target_value = probe_info["target"]["value"]["residues"]
                info_text += f"This probe targets the {target_type}: "
                info_text += f"{target_value}<br>"
            elif probe_info["target"]["type"] == "Sequence":
                target_type = "protein sequence"
                target_value = probe_info["target"]["value"]
                info_text += f"This probe targets the {target_type}: "
                info_text += f"{target_value}<br>"
            else:
                target_type = probe_info["target"]["type"]
                target_value = probe_info["target"]["value"]
                info_text += f"This probe model does not contain a target.<br> If selected, it will be assigned a random target from the selected structure.<br>"
            info_text += f"<b>Probe Model: </b>{probe_info['model']['ID']}<br>"
            probes_gui["probe_info"].value = info_text
        else:
            probes_gui["probe_info"].value = (
                "No information available for this probe."
            )

    def type_dropdown_change(change):
        probes_gui["mock_type_options1"].options = options_per_type1[
            change.new
        ]
        probes_gui["mock_type_options1"].value = options_per_type1[change.new][
            0
        ]
        probes_gui["mock_type_options2"].options = options_per_type2[
            change.new
        ]
        probes_gui["mock_type_options2"].value = options_per_type2[change.new][
            0
        ]

    def clear_probes(b):
        experiment.remove_probes()
        # Clear defect parameters when clearing probes
        experiment.defect_eps["defect"] = 0.0
        experiment.defect_eps["eps1"] = 20.0
        experiment.defect_eps["eps2"] = 100.0
        experiment.defect_eps["use_defects"] = False
        probes_gui["message1"].value = "No probes selected yet."
        probes_gui["message2"].value = "No labelled structure created yet."
        probes_gui["add_probe"].disabled = False
        probes_gui["create_particle"].disabled = True
        update_probe_list()

    # widgets
    probes_gui.add_HTML(tag="Header_message", value="<b>Selected probes:</b>")
    probes_gui.add_HTML(
        "message1",
        "No probes selected yet.",
        style=dict(font_weight="bold", font_size="15px"),
    )
    probes_gui.add_dropdown(
        "select_probe_template",
        description="Choose a probe:",
        options=probe_options,
    )
    probes_gui.add_HTML("probe_info", "")
    probes_gui["select_probe_template"].observe(show_probe_info, names="value")
    probes_gui.add_button(
        "toggle_advanced_parameters",
        description="Toggle advanced parameters",
        icon=toggle_icon,
    )
    # advanced parameters
    probes_gui.add_HTML(
        "advanced_param_header",
        "<b>Advanced parameters</b>",
        style=dict(font_size="15px"),
    )
    probes_gui.add_text(
        tag="probe_name",
        value=probes_gui["select_probe_template"].value,
        description="Probe name",
    )
    probes_gui.add_float_slider(
        "labelling_efficiency",
        description="Labelling efficiency",
        min=0.0,
        max=1.0,
        value=1,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    probes_gui.add_float_slider(
        "distance_from_epitope",
        description="Distance from epitope (Angstroms)",
        min=0.0,
        max=1000,
        value=1,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    # change target type and value
    options_dictionary = dict(
        Protein="Sequence", Residue="Atom_residue", Primary_Probe="Primary"
    )
    probes_gui.add_dropdown(
        "mock_type",
        options=list(options_dictionary.keys()),
        description="I want this probe to target a: ",
    )
    list_of_proteins = experiment.structure.list_protein_names()
    list_of_residues = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
    options_per_type1 = dict(
        Protein=list_of_proteins,
        Residue=list_of_residues,
        Primary_Probe=[
            None,
        ],
    )
    options_per_type2 = dict(
        Protein=["cterminal", "nterminal"],
        Residue=["CA"],
        Primary_Probe=[
            None,
        ],
    )
    probes_gui.add_dropdown(
        "mock_type_options1",
        options=options_per_type1[probes_gui["mock_type"].value],
        description="Which one: ",
    )
    probes_gui.add_dropdown(
        "mock_type_options2",
        options=options_per_type2[probes_gui["mock_type"].value],
        description="Where: ",
    )
    probes_gui.add_HTML(
        "as_linker_info",
        "Activate this option if you intent to use a secondary that recognises the current probe",
    )
    probes_gui.add_checkbox(
        "as_linker",
        description="Model as primary with epitope for secondary probe",
        value=False,
    )
    probes_gui.add_checkbox("wobble", description="Enable wobble", value=False)
    probes_gui.add_float_slider(
        "wobble_theta",
        value=10,
        min=0,
        max=45,
        step=1,
        description="Wobble cone range (degrees)",
    )

    # Defect parameters section
    probes_gui.add_HTML(
        "defects_section_header",
        "<b>Structural Defect Parameters</b>",
        style=dict(font_size="14px", color="darkblue"),
    )
    probes_gui.add_HTML(
        "defects_info",
        "Model structural defects in the macromolecular complex. All three parameters must be set to enable defects.",
        style=dict(font_size="12px", color="gray"),
    )
    probes_gui.add_float_slider(
        "defect_fraction",
        description="Defect fraction (0-1)",
        min=0.0,
        max=1.0,
        value=0.0,
        step=0.001,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    probes_gui.add_float_text(
        "defect_small_cluster",
        description="Small cluster distance (Å)",
        value=100.0,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    probes_gui.add_float_text(
        "defect_large_cluster",
        description="Large cluster distance (Å)",
        value=200.0,
        continuous_update=False,
        style={"description_width": "initial"},
    )

    probes_gui.add_button(
        "add_custom_probe",
        description="Add probe with custom parameters",
        disabled=False,
        icon=add_icon,
    )
    probes_gui["mock_type"].observe(type_dropdown_change, names="value")

    #
    def toggle_advanced_parameters(b):
        probe_widgets_visibility["advanced_param_header"] = (
            not probe_widgets_visibility["advanced_param_header"]
        )
        probe_widgets_visibility["probe_name"] = not probe_widgets_visibility[
            "probe_name"
        ]
        probe_widgets_visibility["labelling_efficiency"] = (
            not probe_widgets_visibility["labelling_efficiency"]
        )
        probe_widgets_visibility["distance_from_epitope"] = (
            not probe_widgets_visibility["distance_from_epitope"]
        )
        probe_widgets_visibility["mock_type"] = not probe_widgets_visibility[
            "mock_type"
        ]
        probe_widgets_visibility["mock_type_options1"] = (
            not probe_widgets_visibility["mock_type_options1"]
        )
        probe_widgets_visibility["mock_type_options2"] = (
            not probe_widgets_visibility["mock_type_options2"]
        )
        probe_widgets_visibility["as_linker_info"] = (
            not probe_widgets_visibility["as_linker_info"]
        )
        probe_widgets_visibility["as_linker"] = not probe_widgets_visibility[
            "as_linker"
        ]
        probe_widgets_visibility["wobble"] = not probe_widgets_visibility[
            "wobble"
        ]
        probe_widgets_visibility["wobble_theta"] = (
            not probe_widgets_visibility["wobble_theta"]
        )
        # Defect parameters visibility
        probe_widgets_visibility["defects_section_header"] = (
            not probe_widgets_visibility["defects_section_header"]
        )
        probe_widgets_visibility["defects_info"] = (
            not probe_widgets_visibility["defects_info"]
        )
        probe_widgets_visibility["defect_fraction"] = (
            not probe_widgets_visibility["defect_fraction"]
        )
        probe_widgets_visibility["defect_small_cluster"] = (
            not probe_widgets_visibility["defect_small_cluster"]
        )
        probe_widgets_visibility["defect_large_cluster"] = (
            not probe_widgets_visibility["defect_large_cluster"]
        )
        probe_widgets_visibility["add_custom_probe"] = (
            not probe_widgets_visibility["add_custom_probe"]
        )
        update_widgets_visibility(probes_gui, probe_widgets_visibility)

    probes_gui.add_callback(
        "add_probe",
        select_probe,
        probes_gui.elements,
        description="Add probe (with defaults)",
        icon=add_icon,
    )
    probes_gui.add_button(
        "clear_probes",
        description="Clear all probes",
        icon=clear_icon,
        style={"button_color": remove_colour},
    )
    probes_gui.add_button(
        "create_particle",
        description="Select probes and create labelled structure",
        style={"button_color": select_colour},
        icon=select_icon,
        disabled=True,
    )
    probes_gui.add_HTML(
        "message2",
        "No labelled structure created yet.",
        style=dict(font_weight="bold", font_size="15px"),
    )
    probe_widgets_visibility = {}
    _unstyle_widgets(probes_gui, probe_widgets_visibility)
    show_probe_info(True)
    probes_gui["create_particle"].on_click(create_particle)
    probes_gui["clear_probes"].on_click(clear_probes)
    probes_gui["toggle_advanced_parameters"].on_click(
        toggle_advanced_parameters
    )
    probes_gui["add_custom_probe"].on_click(select_custom_probe)
    toggle_advanced_parameters(True)  # Initialize with default visibility
    return probes_gui


def ui_select_sample_parameters(experiment):
    """
    Create a widget for configuring sample parameters such as number of particles and random orientations.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing sample parameters.

    Returns
    -------
    EZInput
        Widget for sample parameter configuration.
    """
    sample_gui = EZInput(title="Sample parameters")
    # Add widgets for sample parameters
    sample_gui.add_label("Current sample parameters selected:")

    sample_gui.add_HTML("message", "")

    def update_message():
        if experiment.virtualsample_params.items() is None:
            sample_gui["message"].value = "No sample parameters selected yet."
        else:
            text = ""
            for param in [
                "number_of_particles",
                "random_orientations",
                "minimal_distance",
            ]:
                value = experiment.virtualsample_params[param]
                text += f"{param}: {value}<br>"
            sample_gui["message"].value = text

    sample_gui.add_int_slider(
        "number_of_particles",
        description="Number of particles",
        min=1,
        max=20,
        value=1,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    sample_gui.add_checkbox(
        "random_orientations", description="Randomise orientations", value=True
    )

    sample_gui.add_button(
        "advanced_parameters",
        description="Toggle advanced parameters",
        icon=toggle_icon,
    )
    ####  advanced parameters ####
    sample_gui.add_checkbox(
        "use_min_from_particle",
        value=True,
        description="Use minimal distance from labelled particle dimensions",
    )
    sample_gui.add_bounded_int_text(
        "minimal_distance_nm",
        description="Set minimal distance between particles (nm)",
        value=100,
        vmin=1,
        vmax=1000,
        step=1,
        style={"description_width": "initial"},
    )
    sample_gui.add_file_upload(
        "File",
        description="Select from file",
        accept="*.tif",
        save_settings=False,
    )
    sample_gui.add_bounded_int_text(
        "pixel_size",
        description="Pixel size (nm)",
        value=100,
        vmin=1,
        vmax=1000,
        step=1,
        style={"description_width": "initial"},
    )
    sample_gui.add_bounded_int_text(
        "background_intensity",
        description="Background intensity",
        value=0,
        vmin=0,
        vmax=100000,
        step=1,
        style={"description_width": "initial"},
    )
    sample_gui.add_bounded_int_text(
        "blur_sigma",
        description="Gaussian blurr sigma (nm)",
        value=0,
        vmin=0,
        vmax=1000,
        step=1,
        style={"description_width": "initial"},
    )
    sample_gui.add_bounded_int_text(
        "intensity_threshold",
        description="Intensity threshold",
        value=0,
        vmin=0,
        vmax=10000,
        step=1,
        style={"description_width": "initial"},
    )
    sample_gui.add_dropdown(
        "detection_method",
        description="Detection method",
        options=["Local Maxima", "Mask"],
        value="Local Maxima",
        style={"description_width": "initial"},
    )
    sample_gui.add_checkbox(
        "random",
        value=True,
        description="Randomise positions (enforced when there is more than one particle)",
        style={"description_width": "initial"},
    )
    sample_gui.add_button(
        "update_sample_parameters",
        description="Update sample parameters",
        icon=update_icon,
        style={"button_color": update_colour},
    )
    sample_gui.add_button(
        "select_sample_parameters",
        description="Select parameters and build virtual sample",
        disabled=False,
        icon=select_icon,
        style={"button_color": select_colour},
    )
    sample_gui.add_button(
        "upload_and_set",
        description="Load image and select parameters",
        disabled=False,
        icon=upload_icon,
        style={"button_color": select_colour},
    )
    sample_gui.add_HTML(
        "advanced_params_feedback", "", style=dict(font_weight="bold")
    )

    def update_parameters(b):
        if sample_gui["use_min_from_particle"].value:
            experiment.set_virtualsample_params(
                number_of_particles=sample_gui["number_of_particles"].value,
                random_orientations=sample_gui["random_orientations"].value,
            )
        else:
            min_distance = sample_gui["minimal_distance_nm"].value
            experiment.set_virtualsample_params(
                number_of_particles=sample_gui["number_of_particles"].value,
                random_orientations=sample_gui["random_orientations"].value,
                minimal_distance=min_distance,
            )
        update_message()

    def select_virtual_sample_parameters(b):
        with io.capture_output() as captured:
            experiment.build(modules=["coordinate_field"])
            if experiment.objects_created["imager"]:
                experiment.build(modules=["imager"])
            update_message()

    def upload_and_set(b):
        filepath = sample_gui["File"].selected
        img = tif.imread(filepath)
        pixelsize = sample_gui["pixel_size"].value
        min_distance = None
        if sample_gui["detection_method"].value == "Local Maxima":
            mode = "localmaxima"
        elif sample_gui["detection_method"].value == "Mask":
            mode = "mask"
        else:
            raise ValueError("Unknown detection method selected.")
        if sample_gui["use_min_from_particle"].value:
            min_distance = experiment.virtualsample_params["minimal_distance"]
        else:
            min_distance = sample_gui["minimal_distance_nm"].value
        sigma = sample_gui["blur_sigma"].value
        background = sample_gui["background_intensity"].value
        threshold = sample_gui["intensity_threshold"].value

        npositions = sample_gui["number_of_particles"].value
        experiment.use_image_for_positioning(
            img=img,
            mode=mode,
            sigma=sigma,
            background=background,
            threshold=threshold,
            pixelsize=pixelsize,
            min_distance=min_distance,
            npositions=npositions,
        )
        sample_gui.save_settings()
        update_message()

    def toggle_advanced_parameters(b):
        widgets_visibility["minimal_distance_nm"] = not widgets_visibility[
            "minimal_distance_nm"
        ]
        widgets_visibility["use_min_from_particle"] = not widgets_visibility[
            "use_min_from_particle"
        ]
        widgets_visibility["select_sample_parameters"] = (
            not widgets_visibility["select_sample_parameters"]
        )
        widgets_visibility["upload_and_set"] = not widgets_visibility[
            "upload_and_set"
        ]
        widgets_visibility["File"] = not widgets_visibility["File"]
        widgets_visibility["pixel_size"] = not widgets_visibility["pixel_size"]
        widgets_visibility["background_intensity"] = not widgets_visibility[
            "background_intensity"
        ]
        widgets_visibility["blur_sigma"] = not widgets_visibility["blur_sigma"]
        widgets_visibility["intensity_threshold"] = not widgets_visibility[
            "intensity_threshold"
        ]
        widgets_visibility["detection_method"] = not widgets_visibility[
            "detection_method"
        ]
        widgets_visibility["random"] = not widgets_visibility["random"]
        update_widgets_visibility(sample_gui, widgets_visibility)

    widgets_visibility = {}
    _unstyle_widgets(sample_gui, widgets_visibility)
    update_widgets_visibility(sample_gui, widgets_visibility)
    sample_gui["select_sample_parameters"].on_click(
        select_virtual_sample_parameters
    )
    sample_gui["advanced_parameters"].on_click(toggle_advanced_parameters)
    sample_gui["upload_and_set"].on_click(upload_and_set)
    sample_gui["update_sample_parameters"].on_click(update_parameters)
    widgets_visibility["select_sample_parameters"] = False
    update_parameters(True)
    select_virtual_sample_parameters(
        True
    )  # Initialize with default parameters
    toggle_advanced_parameters(
        True
    )  # Initialize with advanced parameters hidden
    return sample_gui


def ui_select_modality(experiment):
    """
    Create a widget for selecting and previewing imaging modalities.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing modality configuration.

    Returns
    -------
    EZInput
        Widget for modality selection and preview.
    """
    modalities_default = ["Widefield", "Confocal", "STED", "SMLM", "All"]
    imager, preview_experiment = build_virtual_microscope(
        multimodal=modalities_default
    )
    xy_zoom_in = 0.5
    for mod_names in modalities_default[0 : len(modalities_default) - 1]:
        preview_experiment.add_modality(modality_name=mod_names, save=True)
    preview_experiment.build(modules=["imager"])
    modality_gui = EZInput(title="Modality selection")
    modality_gui.add_label("Current modalities list:")
    modality_gui.add_HTML("message", "No modalities selected yet.")

    def update_message():
        text = ""
        for mod_name, params in experiment.imaging_modalities.items():
            text += f"{mod_name}<br>"
        modality_gui["message"].value = text

    def add_modality(b):
        selected_modality = modality_gui["modality"].value
        if selected_modality == "All":
            for mod_names in modalities_default[
                0 : len(modalities_default) - 1
            ]:
                experiment.add_modality(modality_name=mod_names, save=True)
        else:
            experiment.add_modality(modality_name=selected_modality)
        update_message()

    def update_modality_params(b):
        selected_modality = modality_gui["modality"].value
        if selected_modality != "All":
            experiment.update_modality(
                modality_name=selected_modality,
                depth_of_field_nm=modality_gui["psf_depth"].value,
            )
            update_message()
            update_plot(None)

    def remove_modality(b):
        selected_modality = modality_gui["modality"].value
        if selected_modality == "All":
            for mod_names in modalities_default[
                0 : len(modalities_default) - 1
            ]:
                experiment.update_modality(
                    modality_name=mod_names, remove=True
                )
        if selected_modality in experiment.imaging_modalities:
            experiment.update_modality(
                modality_name=selected_modality, remove=True
            )
        else:
            print(f"Modality {selected_modality} not found.")
        update_message()

    def select_modalities(b):
        with io.capture_output() as captured:
            experiment.build(modules=["imager"])
        b1.disabled = True
        b2.disabled = True
        modality_gui["select_modalities"].disabled = True

    def update_plot(change):
        mod_name = modality_gui["modality"].value

        if mod_name != "All":
            if mod_name not in experiment.imaging_modalities:
                info = experiment.local_modalities_parameters[mod_name]
            else:
                info = experiment.imaging_modalities[mod_name]
            psf_stack = preview_experiment.imager.get_modality_psf_stack(
                mod_name
            )
            psf_shape = psf_stack.shape
            half_xy = int(psf_shape[0] / 2)
            half_z = int(psf_shape[2] / 2)
            psf_stack = psf_stack[
                half_xy
                - int(half_xy * xy_zoom_in) : half_xy
                + int(half_xy * xy_zoom_in),
                half_xy
                - int(half_xy * xy_zoom_in) : half_xy
                + int(half_xy * xy_zoom_in),
                :,
            ]
            dimension_plane = modality_gui["dimension_slice"].value
            if dimension_plane == "YZ":
                dimension = 0
            elif dimension_plane == "XZ":
                dimension = 1
            elif dimension_plane == "XY":
                dimension = 2
            # mod info
            pixelsize = info["detector"]["pixelsize"]
            pixelsize_nm = pixelsize * 1000
            psf_voxel = np.array(info["psf_params"]["voxelsize"])
            psf_sd = np.array(info["psf_params"]["std_devs"])
            psf_depth = info["psf_params"]["depth"] * psf_voxel[0]
            s1 = "Detector pixelsize (nm): " + str(pixelsize_nm)
            psf_sd_metric = np.multiply(psf_voxel, psf_sd)
            s2 = "PSF sd (nm): " + str(psf_sd_metric)
            s3 = "Depth of field (nm): " + str(psf_depth)
            s4 = "PSF preview (on a 1x1 µm field of view)"
            modality_gui["modality_info"].value = (
                "<b>"
                + s1
                + "</b><br>"
                + "<b>"
                + s2
                + "</b><br>"
                + "<b>"
                + s3
                + "</b><br>"
                + "<b>"
                + s4
                + "</b><br>"
            )
            modality_gui["preview_modality"].clear_output()
            with modality_gui["preview_modality"]:
                display(
                    slider_normalised(
                        psf_stack,
                        dimension=dimension,
                        cbar=False,
                    )
                )

    modality_gui.add_dropdown(
        "modality",
        description="Modality",
        options=modalities_default,
        on_change=update_plot,
    )
    b1 = widgets.Button(
        description="Add Modality",
        layout=widgets.Layout(width="50%"),
        icon="fa-plus",
    )
    b2 = widgets.Button(
        description="Remove Modality",
        layout=widgets.Layout(width="50%"),
        icon="fa-minus",
    )
    b3 = widgets.Button(
        description="Update modality parameters",
        style={"button_color": "#4985b7d9"},
        layout=widgets.Layout(width="100%"),
        icon="fa-wrench",
    )
    modality_gui.add_custom_widget(
        "add_remove", widgets.HBox, children=[b1, b2]
    )
    select_modalities_button = widgets.Button(
        description="Select modalities and update virtual microscope",
        style={"button_color": "#4daf4ac7"},
        layout=widgets.Layout(width="100%"),
        icon="fa-check",
    )
    modality_gui.add_custom_widget(
        "select_modalities", widgets.HBox, children=[select_modalities_button]
    )
    button_toggle_advanced_parameters = widgets.Button(
        description="Toggle advanced parameters",
        layout=widgets.Layout(width="100%"),
        icon="eye-slash",
    )
    modality_gui.add_custom_widget(
        "toggle_advanced_parameters",
        widgets.HBox,
        children=[button_toggle_advanced_parameters],
    )
    modality_gui.add_int_slider(
        "psf_depth",
        description="PSF depth (nm)",
        min=10,
        max=1000,
        step=10,
        value=100,
        continuous_update=False,
        style={"description_width": "initial"},
    )
    modality_gui.add_custom_widget(
        "update_modality_params",
        widgets.HBox,
        children=[
            b3,
        ],
    )
    button_toggle_preview = widgets.Button(
        description="Toggle modality info and PSF preview",
        layout=widgets.Layout(width="100%"),
        icon="eye-slash",
    )
    modality_gui.add_custom_widget(
        "toggle_preview", widgets.HBox, children=[button_toggle_preview]
    )
    modality_gui.add_HTML("modality_info", "")
    modality_gui.add_custom_widget(
        "dimension_slice",
        widgets.ToggleButtons,
        options=["XY", "XZ", "YZ"],
        value="XY",
        on_change=update_plot,
        style={"description_width": "initial"},
        description="Plane of view: ",
    )
    modality_gui.add_output(
        "preview_modality",
        description="Preview of selected modality",
        style={"description_width": "initial"},
    )

    def toggle_preview(b):
        widgets_visibility["modality_info"] = not widgets_visibility[
            "modality_info"
        ]
        widgets_visibility["dimension_slice"] = not widgets_visibility[
            "dimension_slice"
        ]
        widgets_visibility["preview_modality"] = not widgets_visibility[
            "preview_modality"
        ]
        update_widgets_visibility(modality_gui, widgets_visibility)

    def toggle_advanced_parameters(b):
        widgets_visibility["psf_depth"] = not widgets_visibility["psf_depth"]
        widgets_visibility["update_modality_params"] = not widgets_visibility[
            "update_modality_params"
        ]
        update_widgets_visibility(modality_gui, widgets_visibility)

    widgets_visibility = {}
    _unstyle_widgets(modality_gui, widgets_visibility)
    modality_gui["dimension_slice"].style = dict(description_width="initial")
    b1.on_click(add_modality)
    b2.on_click(remove_modality)
    select_modalities_button.on_click(select_modalities)
    button_toggle_advanced_parameters.on_click(toggle_advanced_parameters)
    b3.on_click(update_modality_params)
    button_toggle_preview.on_click(toggle_preview)
    toggle_preview(True)  # Initialize with default visibility
    toggle_advanced_parameters(True)  # Initialize with default visibility
    update_message()
    update_plot(True)
    return modality_gui


def ui_run_experiment(experiment):
    """
    Create a widget for running the experiment and saving results.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object to run and save.

    Returns
    -------
    EZInput
        Widget for running the experiment and saving results.
    """
    run_gui = EZInput(title="Run experiment")

    # experiment.build(modules=["imager",])
    def run_simulation(b):

        run_gui["message"].value = "Running simulation..."
        run_gui["Acquire"].disabled = True
        sav_dir = run_gui["saving_directory"].value
        if sav_dir is not None:
            experiment.output_directory = sav_dir
            save = True
        experiment.experiment_id = run_gui["experiment_name"].value
        with io.capture_output() as captured:
            output = experiment.run_simulation(save=save)
        run_gui.save_settings()
        if output is None:
            run_gui["message"].value = (
                "Simulation failed. Make sure all parameters are set correctly."
            )
        else:
            run_gui["message"].value = "Simulation completed successfully."
            run_gui["Acquire"].disabled = False

    experiment_info = experiment.current_settings(
        as_string=True, newline="<br>"
    )
    run_gui.add_HTML("experiment_info", experiment_info)
    run_gui.add_label("Set experiment name")
    run_gui.add_text_area(
        "experiment_name", value="Exp_name", remember_value=True
    )
    run_gui.add_label("Set saving directory")
    run_gui.elements["saving_directory"] = FileChooser(
        experiment.output_directory,
        title="<b>Select output directory</b>",
        show_hidden=False,
        select_default=True,
        show_only_dirs=False,
    )
    run_gui.add_HTML("message", "", style=dict(font_weight="bold"))
    run_gui.add_button(
        "Acquire",
        description="Run Simulation",
        icon=select_icon,
        style={"button_color": select_colour},
    )
    run_gui["Acquire"].on_click(run_simulation)
    return run_gui
