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
    ez_sweep_structure.add_dropdown("structures", options=sweep_gen.structures_info_list.keys())
    ez_sweep_structure.add_button("Select", description="Select")

    def select(b):
        sweep_gen.structures = [sweep_gen.structures_info_list[
            ez_sweep_structure["structures"].value
        ],]
        ez_sweep_structure["structures"].disabled = True

    ez_sweep_structure["Select"].on_click(select)
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
        sweep_gen.selected_modalities = widget_modules["modalities"].value
        sweep_gen.selected_probes = widget_modules["probes"].value
        ez_sweep["Select"].disabled = True
        for name in tab_name:
            widget_modules[name].disabled = True

    ez_sweep.add_button("Select", description="Select")
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
    sweep_parameter_gui = EZInput(title="sweep_parameters",)

    for group_name, params in sweep_gen.param_settings.items():
        sweep_parameter_gui.add_HTML(
            tag=group_name,
            value=f"<b>Parameter group: {group_name}</b>",
            style={'font_size': '20px'}
        )
        for param_name, param_info in params.items():
            param_widget = range_widgets[param_name]
            sweep_parameter_gui.elements[param_name] = param_widget
            sweep_parameter_gui.add_label()

    sweep_parameter_gui.add_button("select_parameters",
                                   "Select parameters for sweep")
    sweep_parameter_gui.add_button("clear_parameters",
                                    "Clear all parameters",
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
                    print(f"Setting parameter {param_name} in group {group_name}")
                    if param_info["wtype"] != "logical":
                        start, end = sweep_parameter_gui[param_name].children[2].value
                        steps = sweep_parameter_gui[param_name].children[3].value
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
        reference_structure = reference["structure"].value
        reference_probe = reference["probe"].value
        sweep_gen.reference_structure = reference_structure
        sweep_gen.set_reference_parameters(
            reference_structure=reference_structure,
            reference_probe=reference_probe)
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
        "structure", options=sweep_gen.structures, 
        description="Structure",
        disabled = False
    )
    options_probes = ["NHS_ester",]
    if sweep_gen.structures[0] in probes_per_structure.keys():
        probe_list = probes_per_structure[sweep_gen.structures[0]]
        options_probes.extend(copy.copy(probe_list))
    reference.add_dropdown(
        "probe", options=options_probes, 
        description="Probe",
        disabled = False
    )
    reference.add_dropdown(
        "modality", options=["Reference",], 
        description="Modality",
        disabled = True
    )
    reference.add_button(
        "advanced_parameters",
        description="Toggle advanced parameters",
    )
    # advanced parameters
    reference.add_file_upload(
        "File", description="Select from file", accept="*.tif", save_settings=False
    )



    def toggle_advanced_parameters(b):
        ref_widgets_visibility["File"] = not ref_widgets_visibility["File"]
        update_widgets_visibility(reference, ref_widgets_visibility)

    #
    reference.add_button(
        "set", description="Set reference"
    )
    reference.add_button(
        "preview", description="Preview reference", disabled = True
    )
    reference.elements["feedback"] = widgets.HTML("", style = dict(font_size= "15px", font_weight='bold'))
    reference.add_output(
        "output", description="Reference output"
    )
    
    # visibility and layout
    ref_widgets_visibility = {}
    for wgt in reference.elements.keys():
        ref_widgets_visibility[wgt] = True
        reference.elements[wgt].layout = widgets.Layout(width="50%", display="inline-flex") 

    reference["set"].on_click(gen_ref)
    reference["preview"].on_click(show_reference)
    reference["advanced_parameters"].on_click(toggle_advanced_parameters)
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
        analysis_widget["feedback"].value = "Running analysis sweep. This might take some minutes..."
        analysis_widget["analyse"].disabled = True
        plots = analysis_widget["plots"].value
        param_names_set = sweep_gen.parameters_with_set_values
        if len(param_names_set) >= 2:
            sweep_gen.set_plot_parameters(
                "heatmaps", 
                param1=param_names_set[0], 
                param2=param_names_set[1])
        if analysis_widget["metric"].value == "All":
            metric_list = ["ssim", "pearson"]
        elif analysis_widget["metric"].value == "SSIM":
            metric_list = ["ssim", ]
        elif analysis_widget["metric"].value == "Pearson":
            metric_list = ["pearson", ]
        sweep_gen.set_number_of_repetitions(analysis_widget["reps"].value)
        sweep_gen.set_analysis_parameters(metrics_list = metric_list)
        with io.capture_output() as captured:
            if sweep_gen.reference_image is None:
                sweep_gen.generate_reference_image()
        with analysis_widget["outputs"]:
            print("Generating Virtual samples.")
            print("Once created, a progress bar will show the image simulation progression")
            sweep_gen.run_analysis(plots=plots, save=False)
        analysis_widget["saving_directory"].disabled = False
        analysis_widget["save"].disabled = False
        analysis_widget["output_name"].disabled = False
    def save_results(b):
        output_directory = analysis_widget["saving_directory"].selected_path
        output_name = analysis_widget["output_name"].value
        save_images = analysis_widget["save_images"].value
        sweep_gen.ouput_directory = output_directory
        sweep_gen.save_analysis(
            output_name=output_name
            )
        if save_images:
            sweep_gen.save_images()
    analysis_widget.elements["reps"] = wgen.gen_bound_int(
            value=3, description="Repeats per parameter combination",
            style={'description_width': 'initial'}
        )
    analysis_widget.add_dropdown(
        "metric", options=["SSIM", "Pearson", "All"], 
        description="Metric for image comparison",
        disabled = False
    )
    analysis_widget.add_checkbox("plots", description="Generate plots", value=True)
    analysis_widget.add_button(
        "analyse", description="Run analysis"
    )
    analysis_widget.elements["feedback"] = widgets.HTML("", style = dict(font_size= "15px", font_weight='bold'))
    analysis_widget.elements["outputs"] = widgets.Output()
    analysis_widget.elements["saving_directory"] = FileChooser(
        ouput_directory,
        title="<b>Select output directory</b>",
        show_hidden=False,
        select_default=True,
        show_only_dirs=True,
        disabled=True
    )
    analysis_widget.add_text_area(
        "output_name", 
        value="vlab4mic_analysis", 
        description="Output name")
    analysis_widget.add_checkbox("save_images", description="Save images", value=False)
    analysis_widget.add_button(
        "save", description="save analysis", disabled=True
    )
    analysis_widget["analyse"].on_click(analyse_sweep_action)
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
                    style={'description_width': 'initial'},
                    layout=widgets.Layout(width='40%')
                )
                inttext = wgen.gen_bound_int(
                    value=settings["nintervals"], description="Total values"
                )
                name = widgets.HTML(f"<b>" + parameter_name + "</b>", style={'font_size': '15px'})
                check = widgets.Checkbox(
                            value=False,
                            description="Use parameter",
                            style={'description_width': 'initial'}
                        )
                items = [name, check, slider, inttext]
                range_widgets[parameter_name] = widgets.VBox(
                    items
                )
            elif settings["wtype"] == "logical":
                name = widgets.HTML(f"<b>" + parameter_name + "</b>", style={'font_size': '15px'})
                check = widgets.Checkbox(
                            value=False,
                            description="Use parameter",
                            style={'description_width': 'initial'}
                        )      
                items = [name, check]
                items.append(wgen.gen_logicals(
                    description="Select either True or False or both",
                    layout=widgets.Layout(width='auto', height='auto')
                ))
                range_widgets[parameter_name] = widgets.VBox(
                    items
                )
    return range_widgets