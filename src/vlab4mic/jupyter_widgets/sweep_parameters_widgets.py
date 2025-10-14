"""
sweep_parameters_widgets
------------------------

This module provides widget-based interfaces for configuring parameter sweeps, reference images, and analysis
in the virtual microscopy simulation workflow. It uses EZInput and ipywidgets to allow users to select structures,
probes, modalities, parameter ranges, and to run and save analysis sweeps in Jupyter notebooks.

Functions
---------

- select_structure(sweep_gen):
    Returns a widget for selecting a structure to sweep over.

- select_probes_and_mods(sweep_gen):
    Returns a widget for selecting probes and imaging modalities for the sweep.

- add_parameters_values(sweep_gen):
    Returns a widget for specifying parameter ranges and values to sweep.

- set_reference(sweep_gen):
    Returns a widget for setting and previewing a reference image for analysis.

- analyse_sweep(sweep_gen):
    Returns a widget for running the analysis sweep and saving results.

- create_param_widgets(sweep_gen):
    Helper function to generate widgets for parameter range selection.


Each function returns an EZInput-based widget or ipywidgets element for use in a Jupyter notebook.
"""

import copy
import ipywidgets as widgets
from ezinput import EZInput
from ipyfilechooser import FileChooser
from IPython.utils import io
import matplotlib.pyplot as plt
import numpy as np
from ._widget_generator import widgen
import pprint

select_colour = "#4daf4ac7"
remove_colour = "#ff8000da"
update_colour = "#00bfffda"

select_icon = "fa-check"
add_icon = "fa-plus"
remove_icon = "fa-minus"
loding_icon = "fa-spinner fa-spin"
update_icon = "fa-wrench"  # create
toggle_icon = "fa-eye-slash"
upload_icon = "fa-upload"
view_icon = "fa-eye"
save_icon = "fa-save"


def format_parameter_name(param_name):
    """
    Convert underscore parameter names to natural language format.

    Parameters
    ----------
    param_name : str
        Parameter name with underscores (e.g., 'defect_small_cluster')

    Returns
    -------
    str
        Formatted name in natural language (e.g., 'Defect Small Cluster')
    """
    # Handle specific parameter name mappings
    name_mappings = {
        "labelling_efficiency": "Labeling Efficiency",
        "probe_distance_to_epitope": "Probe Distance to Epitope",
        "defect_fraction": "Defect Fraction",
        "defect_small_cluster": "Defect Small Cluster",
        "defect_large_cluster": "Defect Large Cluster",
        "random_orientations": "Random Orientations",
        "number_of_particles": "Number of Particles",
        "exp_time": "Exposure Time",
    }

    # Return mapped name if available, otherwise convert underscores to spaces
    # and title case
    if param_name in name_mappings:
        return name_mappings[param_name]
    else:
        return param_name.replace("_", " ").title()


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


def select_structure(sweep_gen):
    """
    Create a widget for selecting a structure to sweep over.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing structure information.

    Returns
    -------
    EZInput
        Widget for structure selection.
    """
    ez_sweep_structure = EZInput(title="structure")
    ez_sweep_structure.add_dropdown(
        "structures", options=sweep_gen.structures_info_list.keys()
    )
    ez_sweep_structure.add_HTML(
        "message",
        value="Note: parsing of the structure will be done when running the sweep",
        style={"font_size": "15px"},
    )
    ez_sweep_structure.add_button(
        "Select",
        description="Select",
        icon=select_icon,
        style={"button_color": select_colour},
    )

    def select(b):
        sweep_gen.structures = [
            sweep_gen.structures_info_list[
                ez_sweep_structure["structures"].value
            ],
        ]
        ez_sweep_structure["Select"].disabled = True
        ez_sweep_structure["select_from_database"].disabled = True
        ez_sweep_structure["select_structure_from_file"].disabled = True
    
    def select_structure_from_file(b):
        selected_file = ez_sweep_structure["File"].selected
        filename = ez_sweep_structure["File"].selected_filename
        structure_id = filename.split(".")[0]
        if selected_file:
            sweep_gen.structures = [structure_id,]
            sweep_gen.structure_path = selected_file
            sweep_gen.experiment.select_structure(
                structure_id=structure_id, 
                structure_path=selected_file,
                build=False
            )
            sweep_gen.use_experiment_structure = True
            ez_sweep_structure["message"].value = (
                f"Structure {sweep_gen.structures[0]} selected from file."
            )
            ez_sweep_structure["Select"].disabled = True
            ez_sweep_structure["select_structure_from_file"].disabled = True
        else:
            ez_sweep_structure["message"].value = "No file selected."

    def select_structure_from_database(b):
        structure_id = ez_sweep_structure["structure_id"].value
        database = ez_sweep_structure["database"].value
        if structure_id:
            sweep_gen.structures = [structure_id,]
            sweep_gen.structure_path = None
            sweep_gen.experiment.select_structure(
                structure_id=structure_id, 
                structure_path=None,
                database=database,
                build=True
            )
            sweep_gen.use_experiment_structure = True
            ez_sweep_structure["message"].value = (
                f"Structure {sweep_gen.structures[0]} selected from {database} database."
            )
            ez_sweep_structure["Select"].disabled = True
            ez_sweep_structure["select_from_database"].disabled = True
            ez_sweep_structure["select_structure_from_file"].disabled = True

    
    def toggle_advanced_parameters(b):
        widgets_visibility["advanced_param_header"] = not widgets_visibility[
            "advanced_param_header"
        ]
        widgets_visibility["File"] = not widgets_visibility["File"]
        widgets_visibility["select_structure_from_file"] = not widgets_visibility[
            "select_structure_from_file"
        ]
        #widgets_visibility["database"] = not widgets_visibility["database"]
        #widgets_visibility["structure_id"] = not widgets_visibility["structure_id"]
        #widgets_visibility["select_from_database"] = not widgets_visibility[
        #    "select_from_database"
        #]
        widgets_visibility["divisor1"] = not widgets_visibility["divisor1"]
        #widgets_visibility["divisor2"] = not widgets_visibility["divisor2"]
        widgets_visibility["section1"] = not widgets_visibility["section1"]
        #widgets_visibility["section2"] = not widgets_visibility["section2"]

        update_widgets_visibility(ez_sweep_structure, widgets_visibility)

    ez_sweep_structure.add_button(
        "toggle_advanced_parameters",
        description="Toggle advanced parameters",
        icon=toggle_icon,
    )
    # advanced parameters
    ez_sweep_structure.elements["divisor1"] = widgets.Label("")
    ez_sweep_structure.add_HTML(
        "section1",
        "<b>Select a structure from file</b>",
        style=dict(font_size="15px"),
    )
    ez_sweep_structure.add_HTML(
        "advanced_param_header",
        "<b>Upload a PDB/CIF file</b>",
        style=dict(font_size="15px"),
    )
    ez_sweep_structure.add_file_upload(
        "File",
        description="Select from file",
        accept=["*.pdb", "*.cif"],
        save_settings=False,
    )
    ez_sweep_structure.add_button(
        "select_structure_from_file",
        description="Select structure from file",
        icon=select_icon,
        style={"button_color": select_colour},
    )
    #ez_sweep_structure.elements["divisor2"] = widgets.Label("")
    #ez_sweep_structure.add_HTML(
    #    "section2",
    #    "<b>Or select a structure from the database</b>",
    #    style=dict(font_size="15px"),
    #)
    #ez_sweep_structure.add_dropdown(
    #    "database",
    #    description="Database:",
    #    options=["RCSB", "AlphaFold"]
    #)
    #ez_sweep_structure.add_text(
    #    tag="structure_id",
    #    value="",
    #    description="Structure id",
    #)
    #ez_sweep_structure.add_button(
    #    "select_from_database",
    #    description="Select structure from database",
    #    icon=select_icon,
    #    style={"button_color": select_colour},
    #)
    widgets_visibility = {}
    _unstyle_widgets(ez_sweep_structure, widgets_visibility)
    ez_sweep_structure["toggle_advanced_parameters"].on_click(toggle_advanced_parameters)
    ez_sweep_structure["select_structure_from_file"].on_click(select_structure_from_file)
    #ez_sweep_structure["select_from_database"].on_click(select_structure_from_database)

    ez_sweep_structure["Select"].on_click(select)
    toggle_advanced_parameters(True)
    return ez_sweep_structure


def select_probes_and_mods(sweep_gen):
    """
    Create a widget for selecting probes and imaging modalities for the sweep.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing probe and modality information.

    Returns
    -------
    EZInput
        Widget for probe and modality selection.
    """
    my_exp = sweep_gen.experiment
    probes_per_structure = copy.copy(my_exp.config_probe_per_structure_names)
    vlab_probes = copy.copy(my_exp.config_global_probes_names)
    probe_models = copy.copy(my_exp.config_probe_models_names)
    modalities_default = copy.copy(my_exp.example_modalities)

    ez_sweep = EZInput(title="Sweep")
    probes2show = []
    if sweep_gen.structures[0] in probes_per_structure.keys():
        probe_list = probes_per_structure[sweep_gen.structures[0]]
        probes2show.extend(copy.copy(probe_list))
    probes2show.extend(copy.copy(vlab_probes))
    for probe in probe_models:
        if probe not in probes2show:
            probes2show.append(probe)
    widget_modules = {}
    widget_modules["probes"] = widgets.SelectMultiple(
        description="probes", options=probes2show
    )
    widget_modules["modalities"] = widgets.SelectMultiple(
        description="modalities", options=modalities_default
    )
    tab_name = list(widget_modules.keys())
    children = [widget_modules[name] for name in tab_name]
    ez_sweep.elements["tabs"] = widgets.HBox(children)

    def select_str(b):
        sweep_gen.modalities = list(widget_modules["modalities"].value)
        sweep_gen.probes = list(widget_modules["probes"].value)
        ez_sweep["Select"].disabled = True
        for name in tab_name:
            widget_modules[name].disabled = True

    ez_sweep.add_button(
        "Select",
        description="Select",
        icon=select_icon,
        style={"button_color": select_colour},
    )
    ez_sweep["Select"].on_click(select_str)
    return ez_sweep


def add_parameters_values(sweep_gen):
    """
    Create a widget for specifying parameter ranges and values to sweep.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing parameter settings.

    Returns
    -------
    EZInput
        Widget for parameter range selection and sweep configuration.
    """
    range_widgets = create_param_widgets(sweep_gen)
    sweep_parameter_gui = EZInput(
        title="sweep_parameters",
    )

    for group_name, params in sweep_gen.param_settings.items():
        sweep_parameter_gui.add_HTML(
            tag=group_name,
            value=f"<b>Parameter group: {group_name}</b>",
            style={"font_size": "20px"},
        )
        for param_name, param_info in params.items():
            param_widget = range_widgets[param_name]
            sweep_parameter_gui.elements[param_name] = param_widget
            sweep_parameter_gui.add_label()

    sweep_parameter_gui.add_button(
        "select_parameters",
        "Select parameters for sweep",
        icon=select_icon,
        style={"button_color": select_colour},
    )
    sweep_parameter_gui.add_button(
        "clear_parameters",
        "Clear all parameters",
        icon=remove_icon,
        style={"button_color": remove_colour},
    )
    sweep_parameter_gui.add_HTML(
        tag="message",
        value="No parameters selected",
    )

    def set_param_ranges(b):
        for group_name, params in sweep_gen.param_settings.items():
            for param_name, param_info in params.items():
                use = sweep_parameter_gui[param_name].children[1].value
                if use:
                    print(
                        f"Setting parameter {param_name} in group {group_name}"
                    )
                    if param_info["wtype"] != "logical":
                        start, end = (
                            sweep_parameter_gui[param_name].children[2].value
                        )
                        steps = (
                            sweep_parameter_gui[param_name].children[3].value
                        )
                        param_values = (start, end, steps)
                    else:
                        val = sweep_parameter_gui[param_name].children[2].value
                        if val == "Both":
                            param_values = [True, False]
                        elif val == "True":
                            param_values = [True]
                        else:
                            param_values = [False]
                    sweep_gen.set_parameter_values(
                        param_group=group_name,
                        param_name=param_name,
                        values=param_values,
                    )

        sweep_parameter_gui["message"].value = "Parameters set successfully"

    def clear_parameters(b):
        sweep_gen.clear_sweep_parameters()
        sweep_parameter_gui["message"].value = "All parameters cleared"

    sweep_parameter_gui["select_parameters"].on_click(set_param_ranges)
    sweep_parameter_gui["clear_parameters"].on_click(clear_parameters)
    return sweep_parameter_gui


def set_reference(sweep_gen):
    """
    Create a widget for setting and previewing a reference image for analysis.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing experiment and reference info.

    Returns
    -------
    EZInput
        Widget for reference image selection and preview.
    """
    my_exp = sweep_gen.experiment
    probes_per_structure = copy.copy(my_exp.config_probe_per_structure_names)
    reference = EZInput(title="reference")

    def gen_ref(b):
        reference["set"].disabled = True
        reference["feedback"].value = "Generating Reference..."
        if sweep_gen.use_experiment_structure:
            #reference_structure = reference["structure"].value
            sweep_gen.reference_structure = copy.copy(sweep_gen.structure_path)
            with io.capture_output() as captured:
                sweep_gen.generate_reference_image(override=True)
            reference["feedback"].value = "Reference Set"
            reference["preview"].disabled = False
        else:
            reference_structure = reference["structure"].value
            reference_probe = reference["probe"].value
            sweep_gen.reference_structure = reference_structure
            sweep_gen.set_reference_parameters(
                reference_structure=reference_structure,
                reference_probe=reference_probe,
            )
            with io.capture_output() as captured:
                sweep_gen.generate_reference_image(override=True)
            reference["feedback"].value = "Reference Set"
            reference["preview"].disabled = False

    def show_reference(b):
        reference["output"].clear_output()
        with reference["output"]:
            image = sweep_gen.preview_reference_image(return_image=True)
            if image.shape[0] == 1:
                image = np.squeeze(image)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(image)
            plt.close()
            display(fig)

    reference.add_dropdown(
        "structure",
        options=sweep_gen.structures,
        description="Structure",
        disabled=False,
    )
    options_probes = [
        "NHS_ester",
    ]
    if sweep_gen.structures[0] in probes_per_structure.keys():
        probe_list = probes_per_structure[sweep_gen.structures[0]]
        options_probes.extend(copy.copy(probe_list))
    reference.add_dropdown(
        "probe", options=options_probes, description="Probe", disabled=False
    )
    reference.add_dropdown(
        "modality",
        options=[
            "Reference",
        ],
        description="Modality",
        disabled=True,
    )
    reference.add_button(
        "advanced_parameters",
        description="Toggle advanced parameters",
        icon=toggle_icon,
    )
    # advanced parameters
    reference.add_HTML(
        tag="Upload_ref_message",
        value="<b>Upload image for reference</b>",
        style={"font_size": "15px"},
    )
    reference.add_HTML(
        tag="File_message",
        value="Choose reference image",
        style={"font_size": "15px"},
    )
    reference.add_HTML(
        tag="File_message",
        value="Choose reference image",
        style={"font_size": "15px"},
    )
    reference.add_file_upload(
        "File",
        description="Select from file",
        accept=["*.tif", "*.tiff"],
        save_settings=False,
    )
    reference.add_HTML(
        tag="File_message_mask",
        value="(Optional) Choose reference image mask:",
        style={"font_size": "15px"},
    )
    reference.add_file_upload(
        "File_mask",
        description="(Optional) Select mask from file",
        accept=["*.tif", "*.tiff"],
        save_settings=False,
    )
    reference.add_bounded_float_text(
        "pixel_size",
        description="Pixel size (nm)",
        value=100,
        vmin=1,
        vmax=1000,
        step=0.1,
        style={"description_width": "initial"},
    )
    reference.add_button(
        "upload_and_set",
        description="Upload image reference",
        disabled=False,
        icon=upload_icon,
    )

    def toggle_advanced_parameters(b):
        ref_widgets_visibility["Upload_ref_message"] = (
            not ref_widgets_visibility["Upload_ref_message"]
        )
        ref_widgets_visibility["File"] = not ref_widgets_visibility["File"]
        ref_widgets_visibility["upload_and_set"] = not ref_widgets_visibility[
            "upload_and_set"
        ]
        ref_widgets_visibility["pixel_size"] = not ref_widgets_visibility[
            "pixel_size"
        ]
        ref_widgets_visibility["File_mask"] = not ref_widgets_visibility[
            "File_mask"
        ]
        ref_widgets_visibility["File_message"] = not ref_widgets_visibility[
            "File_message"
        ]
        ref_widgets_visibility["File_message_mask"] = not ref_widgets_visibility[
            "File_message_mask"
        ]
        update_widgets_visibility(reference, ref_widgets_visibility)

    def upload_and_set(b):
        if reference["File_mask"].selected:
            ref_image_mask_path = reference["File_mask"].selected
        else:
            ref_image_mask_path = None
        sweep_gen.load_reference_image(
            ref_image_path=reference["File"].selected,
            ref_pixelsize=reference["pixel_size"].value,
            ref_image_mask_path=ref_image_mask_path
        )
        reference["feedback"].value = "Reference image uploaded and set."
        reference["preview"].disabled = False

    #
    reference.add_button(
        "set",
        description="Generate image reference",
        icon=select_colour,
        style={"button_color": select_colour},
    )
    reference.add_button(
        "preview",
        description="Preview reference",
        disabled=True,
        icon=view_icon,
    )
    reference.elements["feedback"] = widgets.HTML(
        "", style=dict(font_size="15px", font_weight="bold")
    )
    reference.add_output("output", description="Reference output")

    # visibility and layout
    ref_widgets_visibility = {}
    _unstyle_widgets(reference, ref_widgets_visibility)
    reference["set"].on_click(gen_ref)
    reference["preview"].on_click(show_reference)
    reference["advanced_parameters"].on_click(toggle_advanced_parameters)
    reference["upload_and_set"].on_click(upload_and_set)
    toggle_advanced_parameters(True)
    return reference


def analyse_sweep(sweep_gen):
    """
    Create a widget for running the analysis sweep and saving results.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing analysis and output info.

    Returns
    -------
    EZInput
        Widget for running analysis and saving results.
    """
    wgen = widgen()
    ouput_directory = getattr(sweep_gen, "ouput_directory", ".")
    analysis_widget = EZInput(title="analysis")

    def analyse_sweep_action(b):
        analysis_widget["feedback"].value = (
            "Running analysis sweep. This might take some minutes..."
        )
        analysis_widget["preview"].disabled = True
        analysis_widget["analyse"].disabled = True
        update_sliders(from_sweep=False, disable_all=True)
        plots = analysis_widget["plots"].value
        param_names_set = sweep_gen.parameters_with_set_values
        if len(param_names_set) >= 2:
            sweep_gen.set_plot_parameters(
                "heatmaps",
                param1=param_names_set[0],
                param2=param_names_set[1],
            )
        if analysis_widget["metric"].value == "All":
            metric_list = ["ssim", "pearson"]
        elif analysis_widget["metric"].value == "SSIM":
            metric_list = [
                "ssim",
            ]
        elif analysis_widget["metric"].value == "Pearson":
            metric_list = [
                "pearson",
            ]
        sweep_gen.set_number_of_repetitions(analysis_widget["reps"].value)
        sweep_gen.set_analysis_parameters(metrics_list=metric_list)
        with io.capture_output() as captured:
            if sweep_gen.reference_image is None:
                sweep_gen.generate_reference_image()
        with analysis_widget["outputs"]:
            print("Generating Virtual samples.")
            print(
                "Once created, a progress bar will show the image simulation progression"
            )
            sweep_gen.run_analysis(plots=plots, save=False)
        analysis_widget["preview"].disabled = False
        analysis_widget["saving_directory"].disabled = False
        analysis_widget["save"].disabled = False
        analysis_widget["output_name"].disabled = False
        analysis_widget["analyse"].disabled = False
        update_sliders()
        #
    def update_sliders(from_sweep=True, disable_all=True):
        if from_sweep:
            flag = True
        elif disable_all:
            flag = False
        else:
            flag = True
        n_probes = len(sweep_gen.probes)
        n_modalities = len(sweep_gen.modalities)
        if n_modalities > 1 and flag:
            analysis_widget["modality_template"].max = n_modalities - 1
            analysis_widget["modality_template"].disabled = False
        else:
            analysis_widget["modality_template"].disabled = True
        if n_probes > 1 and flag:
            analysis_widget["probe_template"].max = n_probes - 1
            analysis_widget["probe_template"].disabled = False
        else:
            analysis_widget["probe_template"].disabled = True
        if sweep_gen.probe_parameters is not None and flag:
            n_probe_parameters = len(sweep_gen.probe_parameters.keys())
            analysis_widget["probe_parameters"].max = n_probe_parameters - 1
            analysis_widget["probe_parameters"].disabled = False
        else:
            analysis_widget["probe_parameters"].disabled = True
        if sweep_gen.defect_parameters is not None and flag:
            n_defect_parameters = len(sweep_gen.defect_parameters.keys())
            analysis_widget["defect_parameters"].max = n_defect_parameters - 1
            analysis_widget["defect_parameters"].disabled = False
        else:
            analysis_widget["defect_parameters"].disabled = True
        if sweep_gen.vsample_parameters is not None and flag:
            n_vsample_parameters = len(sweep_gen.vsample_parameters.keys())
            analysis_widget["vsample_parameters"].max = n_vsample_parameters - 1
            analysis_widget["vsample_parameters"].disabled = False
        else:
            analysis_widget["vsample_parameters"].disabled = True
        if sweep_gen.acquisition_parameters is not None and flag:
            n_acquisition_parameters = len(
                sweep_gen.acquisition_parameters.keys()
            )
            analysis_widget["acquisition_parameters"].max = (
                n_acquisition_parameters - 1
            )
            analysis_widget["acquisition_parameters"].disabled = False
        else:
            analysis_widget["acquisition_parameters"].disabled = True
        n_replicas = sweep_gen.sweep_repetitions
        if n_replicas >= 2 and flag:
            analysis_widget["replica_number"].max = n_replicas - 1
            analysis_widget["replica_number"].disabled = False
        else:
            analysis_widget["replica_number"].disabled = True
        
    def save_results(b):
        output_directory = analysis_widget["saving_directory"].selected_path
        output_name = analysis_widget["output_name"].value
        save_images = analysis_widget["save_images"].value
        sweep_gen.ouput_directory = output_directory
        sweep_gen.save_analysis(output_name=output_name)
        if save_images:
            sweep_gen.save_images()

    analysis_widget.elements["reps"] = wgen.gen_bound_int(
        value=3,
        description="Repeats per parameter combination",
        style={"description_width": "initial"},
    )
    analysis_widget.add_dropdown(
        "metric",
        options=["SSIM", "Pearson", "All"],
        description="Metric for image comparison",
        disabled=False,
    )
    analysis_widget.add_checkbox(
        "plots", description="Generate plots", value=True
    )
    analysis_widget.add_button(
        "analyse",
        description="Run analysis",
        icon=select_colour,
        style={"button_color": select_colour},
    )
    analysis_widget.elements["feedback"] = widgets.HTML(
        "", style=dict(font_size="15px", font_weight="bold")
    )
    analysis_widget.elements["outputs"] = widgets.Output()
    analysis_widget.elements["saving_directory"] = FileChooser(
        ouput_directory,
        title="<b>Select output directory</b>",
        show_hidden=False,
        select_default=True,
        show_only_dirs=True,
        disabled=True,
    )
    analysis_widget.add_text_area(
        "output_name", value="vlab4mic_analysis", description="Output name"
    )

    # preview
    def update_plot(change):
        modality_template = analysis_widget["modality_template"].value
        probe_template = analysis_widget["probe_template"].value
        probe_parameters = analysis_widget["probe_parameters"].value
        defect_parameters = analysis_widget["defect_parameters"].value
        vsample_parameters = analysis_widget["vsample_parameters"].value
        acquisition_parameters = analysis_widget[
            "acquisition_parameters"
        ].value
        replica_number = analysis_widget["replica_number"].value

        image, parameters = sweep_gen.preview_image_output_by_ID(
            modality_template=modality_template,
            probe_template=probe_template,
            probe_parameters=probe_parameters,
            defect_parameters=defect_parameters,
            virtual_sample_parameters=vsample_parameters,
            acquisition_parameters=acquisition_parameters,
            replica_number=replica_number,
            return_image=True,
        )
        figure, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image, cmap="gray")
        ax.axis("off")
        ax.set_title(f"Preview of parameter sweeps results")
        plt.close()
        analysis_widget["preview_results"].clear_output()
        with analysis_widget["preview_results"]:
            display(figure)
        text = "Structure: " + str(parameters[0]) + "<br>"
        text += "Modality: " + str(parameters[5])+ "<br>"
        text += "Probe: " + str(parameters[1])+ "<br>"
        text += "Probe Parameters: " + str(parameters[2])+ "<br>"
        text += "Defect Parameters: " + str(parameters[3])+ "<br>"
        text += "Virtual Sample Parameters: " + str(parameters[4])+ "<br>"
        text += "Acquisition Parameters: " + str(parameters[6])+ "<br>"
        analysis_widget["params_preview"].value = text

    def toggle_preview(b):
        widgets_visibility["preview_results"] = not widgets_visibility[
            "preview_results"
        ]
        widgets_visibility["params_preview"] = not widgets_visibility[
            "params_preview"
        ]
        widgets_visibility["modality_template"] = not widgets_visibility[
            "modality_template"
        ]
        widgets_visibility["probe_template"] = not widgets_visibility[
            "probe_template"
        ]
        widgets_visibility["probe_parameters"] = not widgets_visibility[
            "probe_parameters"
        ]
        widgets_visibility["defect_parameters"] = not widgets_visibility[
            "defect_parameters"
        ]
        widgets_visibility["vsample_parameters"] = not widgets_visibility[
            "vsample_parameters"
        ]
        widgets_visibility["acquisition_parameters"] = not widgets_visibility[
            "acquisition_parameters"
        ]
        widgets_visibility["replica_number"] = not widgets_visibility[
            "replica_number"
        ]

        update_widgets_visibility(analysis_widget, widgets_visibility)
        update_plot(widgets_visibility["preview_results"])

    # preview widgets
    if sweep_gen.acquisition_outputs is None:
        preview_disabled = True
    else:
        preview_disabled = False
    analysis_widget.add_button(
        "preview", description="Preview results", disabled=preview_disabled, icon=view_icon
    )
    analysis_widget.add_int_slider(
        "modality_template",
        description="Modality",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "probe_template",
        description="Probe",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "probe_parameters",
        description="Probe parameter",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "defect_parameters",
        description="Defect parameter",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "vsample_parameters",
        description="Vsample parameters",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "acquisition_parameters",
        description="acquisition parameters",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    analysis_widget.add_int_slider(
        "replica_number",
        description="Replica number",
        min=0,
        max=1,
        value=0,
        continuous_update=False,
    )
    # connect the preview widgets to the update function
    analysis_widget["modality_template"].observe(update_plot, names="value")
    analysis_widget["probe_template"].observe(update_plot, names="value")
    analysis_widget["probe_parameters"].observe(update_plot, names="value")
    analysis_widget["defect_parameters"].observe(update_plot, names="value")
    analysis_widget["vsample_parameters"].observe(update_plot, names="value")
    analysis_widget["acquisition_parameters"].observe(
        update_plot, names="value"
    )
    analysis_widget["replica_number"].observe(update_plot, names="value")
    # output preview
    analysis_widget.add_output(
        "preview_results", description="Preview results"
    )
    analysis_widget.add_HTML(tag="params_preview", value="")
    # save
    analysis_widget.add_checkbox(
        "save_images", description="Save images", value=False
    )
    analysis_widget.add_button(
        "save", description="save analysis", disabled=True, icon=save_icon
    )
    widgets_visibility = {}
    _unstyle_widgets(analysis_widget, widgets_visibility)
    widgets_visibility["preview_results"] = False
    widgets_visibility["params_preview"] = False
    widgets_visibility["modality_template"] = False
    widgets_visibility["probe_template"] = False
    widgets_visibility["probe_parameters"] = False
    widgets_visibility["defect_parameters"] = False
    widgets_visibility["vsample_parameters"] = False
    widgets_visibility["acquisition_parameters"] = False
    widgets_visibility["replica_number"] = False
    update_widgets_visibility(analysis_widget, widgets_visibility)
    update_sliders()
    analysis_widget["analyse"].on_click(analyse_sweep_action)
    analysis_widget["preview"].on_click(toggle_preview)
    analysis_widget["save"].on_click(save_results)
    return analysis_widget


def create_param_widgets(sweep_gen):
    """
    Helper function to generate widgets for parameter range selection.

    Parameters
    ----------
    sweep_gen : object
        The sweep generator object containing parameter settings.

    Returns
    -------
    dict
        Dictionary of widgets for each parameter.
    """
    wgen = widgen()
    range_widgets = dict()
    for groupname, group_parameters in sweep_gen.param_settings.items():
        for parameter_name, settings in group_parameters.items():
            if (
                settings["wtype"] == "float_slider"
                or settings["wtype"] == "int_slider"
            ):
                if settings["wtype"] == "float_slider":
                    slidertype = "float"
                else:
                    slidertype = "int"
                slider = wgen.gen_range_slider(
                    slidertype=slidertype,
                    minmaxstep=settings["range"],
                    orientation="horizontal",
                    description="Range",
                    style={"description_width": "initial"},
                    layout=widgets.Layout(width="40%"),
                )
                inttext = wgen.gen_bound_int(
                    value=settings["nintervals"], description="Total values"
                )
                # Use formatted parameter name for display
                formatted_name = format_parameter_name(parameter_name)
                name = widgets.HTML(
                    f"<b>{formatted_name}</b>", style={"font_size": "15px"}
                )
                check = widgets.Checkbox(
                    value=False,
                    description="Use parameter",
                    style={"description_width": "initial"},
                )
                items = [name, check, slider, inttext]
                range_widgets[parameter_name] = widgets.VBox(items)
            elif settings["wtype"] == "logical":
                # Use formatted parameter name for display
                formatted_name = format_parameter_name(parameter_name)
                name = widgets.HTML(
                    f"<b>{formatted_name}</b>", style={"font_size": "15px"}
                )
                check = widgets.Checkbox(
                    value=False,
                    description="Use parameter",
                    style={"description_width": "initial"},
                )
                items = [name, check]
                items.append(
                    wgen.gen_logicals(
                        description="Select either True or False or both",
                        layout=widgets.Layout(width="auto", height="auto"),
                    )
                )
                range_widgets[parameter_name] = widgets.VBox(items)
    return range_widgets
