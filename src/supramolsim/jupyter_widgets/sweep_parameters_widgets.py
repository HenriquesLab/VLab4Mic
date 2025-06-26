import copy
import ipywidgets as widgets
from ezinput import EZInput
from ipyfilechooser import FileChooser
from IPython.utils import io
import matplotlib.pyplot as plt
import numpy as np
from .widget_generator import widgen

def select_structure(sweep_gen):
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
    my_exp = sweep_gen.experiment
    probes_per_structure = copy.copy(my_exp.config_probe_per_structure_names)
    vlab_probes = copy.copy(my_exp.config_global_probes_names)
    modalities_default = copy.copy(my_exp.example_modalities)

    ez_sweep = EZInput(title="Sweep")
    probes2show = []
    if sweep_gen.structures[0] in probes_per_structure.keys():
        probe_list = probes_per_structure[sweep_gen.structures[0]]
        probes2show.extend(copy.copy(probe_list))
    probes2show.extend(copy.copy(vlab_probes))
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

    range_widgets = create_param_widgets(sweep_gen)
    param_ranges = EZInput(title="ranges")
    
    def change_param_list(change):
        new_options = list(sweep_gen.param_settings[change.new].keys())
        param_ranges["parms_per_group"].options = new_options

    def change_param_widget(change):
        param_ranges[change.old].layout.display = "None"
        param_ranges[change.new].layout.display = "inline-flex"

    def set_param_range(b):
        param_group = param_ranges["groups"].value
        param_name = param_ranges["parms_per_group"].value
        if sweep_gen.param_settings[param_group][param_name]["wtype"] != "logical":
            start, end = param_ranges[param_name].children[0].value
            steps = param_ranges[param_name].children[1].value
            param_values = (start, end, steps)
        else:
            param_values = []
            if param_ranges[param_name].value == "Both":
                param_values = [True, False,]
            elif param_ranges[param_name].value ==  "True":
                param_values = [True,]
            if param_ranges[param_name].value == "False":
                param_values = [False,]

        sweep_gen.set_parameter_values(
            param_group=param_group,
            param_name=param_name,
            values=param_values,
        )
    
    def disable_widgets(b):
        param_ranges["groups"].disabled = True
        param_ranges["parms_per_group"].disabled = True
        param_ranges["add_parameter"].disabled = True
        param_ranges["done"].disabled = True

    parameter_group_names = list(sweep_gen.param_settings.keys())
    param_ranges.add_dropdown("groups", options=parameter_group_names, description="Parameter group")
    param_ranges.add_dropdown(
        "parms_per_group",
        options=list(sweep_gen.param_settings[param_ranges["groups"].value].keys()),
        description="Parameter name"
    )
    for wname, wgt in range_widgets.items():
        param_ranges.elements[wname] = wgt
        param_ranges.elements[wname].layout.display = "None"
    param_ranges[param_ranges["parms_per_group"].value].layout.display = (
        "inline-flex"
    )
    param_ranges.add_button(
        "add_parameter", description="Add this parameter for sweep"
    )
    param_ranges.add_button(
        "done", description="Done"
    )
    param_ranges["groups"].observe(change_param_list, names="value")
    param_ranges["parms_per_group"].observe(change_param_widget, names="value")
    param_ranges["add_parameter"].on_click(set_param_range)
    param_ranges["done"].on_click(disable_widgets)
    return param_ranges

def set_reference(sweep_gen):
    my_exp = sweep_gen.my_experiment
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
            sweep_gen.preview_reference_image()

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
        "set", description="Set reference"
    )
    reference.add_button(
        "preview", description="Preview reference", disabled = True
    )
    reference.elements["feedback"] = widgets.HTML("", style = dict(font_size= "15px", font_weight='bold'))
    reference.elements["output"] = widgets.Output()
    reference["set"].on_click(gen_ref)
    reference["preview"].on_click(show_reference)
    return reference

def analyse_sweep(sweep_gen):
    wgen = sweep_gen.wgen
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
                    description=parameter_name,
                    style={'description_width': 'initial'}
                )
                inttext = wgen.gen_bound_int(
                    value=settings["nintervals"], description="Total values"
                )
                range_widgets[parameter_name] = wgen.gen_box(
                    widget1=slider,
                    widget2=inttext,
                    orientation="vertical",
                    layout = widgets.Layout(width="70%")
                )
            elif settings["wtype"] == "logical":
                range_widgets[parameter_name] = wgen.gen_logicals()
    return range_widgets