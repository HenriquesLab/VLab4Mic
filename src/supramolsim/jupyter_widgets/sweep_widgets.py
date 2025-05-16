import supramolsim.utils
from ..sweep_generator import sweep_generator
import ipywidgets as widgets
from .widget_generator import widgen
from .widgets_dataclass import jupyter_gui
import os
import supramolsim
from supramolsim.utils.io import yaml_functions
import copy
from ezinput import EZInput
from ipyfilechooser import FileChooser

class Sweep_gui(jupyter_gui):
    sweep_gen = sweep_generator()
    wgen = widgen()
    selected_structure = None
    selected_probes = None
    selected_modalities = None
    vlab_probes = []
    targetless_probes = []
    probes_per_structure = {}
    probe_parameters = {}
    # widget parameters
    param_settings = None
    range_widgets = {}

    def __init__(self):
        super().__init__()
        self.vlab_probes = []
        param_settings_file = os.path.join(
            self.config_directories["base"], "parameter_settings.yaml"
        )
        self.param_settings = yaml_functions.load_yaml(param_settings_file)
        self._create_param_widgets()
        for file in os.listdir(self.config_directories["probes"]):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                label_config_path = os.path.join(
                    self.config_directories["probes"], file
                )
                label_parmeters = supramolsim.load_yaml(label_config_path)
                lablname = os.path.splitext(file)[0]
                # self.vlab_probes.append(lablname)
                # self.probe_parameters[lablname] = label_parmeters
                if "Mock" in label_parmeters["known_targets"]:
                    self.targetless_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                if "Generic" in label_parmeters["known_targets"]:
                    self.vlab_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                else:
                    #self.vlab_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                    for known_structures in label_parmeters["known_targets"]:
                        if known_structures in self.probes_per_structure.keys():
                            self.probes_per_structure[known_structures].append(lablname)
                        else:
                            self.probes_per_structure[known_structures] = [
                                lablname,
                            ]

    def _create_param_widgets(self):
        for groupname, group_parameters in self.param_settings.items():
            for parameter_name, settings in group_parameters.items():
                if (
                    settings["wtype"] == "float_slider"
                    or settings["wtype"] == "int_slider"
                ):
                    if settings["wtype"] == "float_slider":
                        slidertype = "float"
                    else:
                        slidertype = "int"
                    slider = self.wgen.gen_range_slider(
                        slidertype=slidertype,
                        minmaxstep=settings["range"],
                        orientation="horizontal",
                        description=parameter_name,
                        style={'description_width': 'initial'}
                    )
                    inttext = self.wgen.gen_bound_int(
                        value=settings["nintervals"], description="Total values"
                    )
                    self.range_widgets[parameter_name] = self.wgen.gen_box(
                        widget1=slider,
                        widget2=inttext,
                        orientation="vertical",
                        layout = widgets.Layout(width="70%")
                    )
                elif settings["wtype"] == "logical":
                    self.range_widgets[parameter_name] = self.wgen.gen_logicals()

    def select_structure(self):
        ez_sweep_structure = EZInput(title="structure")
        ez_sweep_structure.add_dropdown("structures", options=self.demo_structures)
        ez_sweep_structure.add_button("Select", description="Select")

        def select(b):
            self.selected_structure = self.structures_info_list[
                ez_sweep_structure["structures"].value
            ]
            ez_sweep_structure["structures"].disabled = True

        ez_sweep_structure["Select"].on_click(select)
        ez_sweep_structure.show()

    def select_probes_and_mods(self, include_probe_models=False):
        ez_sweep = EZInput(title="Sweep")
        probes2show = []
        if self.selected_structure:
            probe_list = self.probes_per_structure[self.selected_structure]
            probes2show.extend(
                copy.copy(probe_list)
            )
        probes2show.extend(
                copy.copy(self.vlab_probes)
            )
        if include_probe_models:
            probes2show.extend(copy.copy(self.targetless_probes))
        # create muliple options widgets
        widget_modules = {}
        widget_modules["probes"] = widgets.SelectMultiple(
            description="", options=probes2show
        )
        widget_modules["modalities"] = widgets.SelectMultiple(
            description="", options=self.modalities_default
        )
        # create tabs
        tab_name = list(widget_modules.keys())
        children = [widget_modules[name] for name in tab_name]
        ez_sweep.elements["tabs"] = widgets.Tab()
        ez_sweep.elements["tabs"].children = children
        ez_sweep.elements["tabs"].titles = tab_name

        # on clicks
        def select_str(b):
            self.selected_modalities = widget_modules["modalities"].value
            self.selected_probes = widget_modules["probes"].value
            ez_sweep["Select"].disabled = True
            for name in tab_name:
                widget_modules[name].disabled = True
            #ez_sweep["modalities"].disabled = True
            #ez_sweep["probes"].disabled = True

        ez_sweep.add_button("Select", description="Select")
        ez_sweep["Select"].on_click(select_str)
        ez_sweep.show()

    def add_parameters_ranges(self):
        param_ranges = EZInput(title="ranges")

        def change_param_list(change):
            new_options = list(self.param_settings[change.new].keys())
            param_ranges["parms_per_group"].options = new_options

        def change_param_widget(change):
            param_ranges[change.old].layout.display = "None"
            param_ranges[change.new].layout.display = "inline-flex"

        def set_param_range(b):
            param_group = param_ranges["groups"].value
            param_name = param_ranges["parms_per_group"].value
            if self.param_settings[param_group][param_name]["wtype"] != "logical":
                param_type = "numeric"
                first, last = param_ranges[param_name].children[0].value
                option = param_ranges[param_name].children[1].value
            else:
                param_type = "logical"
                first = None
                last = None
                option = param_ranges[param_name].value
            self.sweep_gen._set_param_range(
                param_group=param_group,
                param_name=param_name,
                param_type=param_type,
                first=first,
                last=last,
                option=option,
            )
        
        def disable_widgets(b):
            param_ranges["groups"].disabled = True
            param_ranges["parms_per_group"].disabled = True
            param_ranges["add_parameter"].disabled = True
            param_ranges["done"].disabled = True

        parameter_group_names = list(self.param_settings.keys())
        param_ranges.add_dropdown("groups", options=parameter_group_names, description="Parameter group")
        param_ranges.add_dropdown(
            "parms_per_group",
            options=list(self.param_settings[param_ranges["groups"].value].keys()),
            description="Parameter name"
        )
        # add the widgets to list
        for wname, wgt in self.range_widgets.items():
            param_ranges.elements[wname] = wgt
            param_ranges.elements[wname].layout.display = "None"
        # show the first one
        param_ranges[param_ranges["parms_per_group"].value].layout.display = (
            "inline-flex"
        )
        param_ranges.add_button(
            "add_parameter", description="Add this parameter for sweep"
        )
        param_ranges.add_button(
            "done", description="Done"
        )
        # widget actions or updates
        param_ranges["groups"].observe(change_param_list, names="value")
        param_ranges["parms_per_group"].observe(change_param_widget, names="value")
        param_ranges["add_parameter"].on_click(set_param_range)
        param_ranges["done"].on_click(disable_widgets)
        param_ranges.show()


    def generate_simulations(self):
        simulate = EZInput(title="simulate")
        def run_sweeps(b):
            simulate["Run"].disabled = True
            self.sweep_gen.sweep_repetitions = simulate["reps"].value
            self.sweep_gen.create_parameters_iterables()
            with simulate["outputs"]:
                self.sweep_gen.generate_acquisitions()
        simulate.elements["reps"] = self.wgen.gen_bound_int(
                value=3, description="Repeats per parameter combination",
                style={'description_width': 'initial'}
            )
        simulate.add_button(
            "Run", description="Run"
        )
        simulate.elements["outputs"] = widgets.Output()
        simulate["Run"].on_click(run_sweeps)
        simulate.show()

    def set_reference(self):
        reference = EZInput(title="reference")
        def gen_ref(b):
            reference["set"].disabled = True
            self.reference_structure = self.select_structure
            self.sweep_gen.generate_reference_image()
        reference.add_dropdown(
            "structure", options=["same for parameter sweep",], 
            description="Structure",
            disabled = True
        )
        reference.add_dropdown(
            "probe", options=["NHS_ester",], 
            description="Probe",
            disabled = True
        )
        reference.add_dropdown(
            "modality", options=["Reference",], 
            description="Modality",
            disabled = True
        )
        reference.add_button(
            "set", description="Set reference"
        )
        reference["set"].on_click(gen_ref)
        reference.show()
        
    def analyse_sweep(self):
        analysis_widget = EZInput(title="analysis")
        def analyse_sweep(b):
            analysis_widget["analyse"].disabled = True
            
            self.sweep_gen.run_analysis()
            self.sweep_gen.gen_analysis_dataframe()
            analysis_widget["saving_directory"].disabled = False
            analysis_widget["save"].disabled = False
            analysis_widget["output_name"].disabled = False
        def save_results(b):
            output_name = analysis_widget["output_name"].value
            self.sweep_gen.save_analysis(
                output_directory=self.ouput_directory,
                output_name=output_name
                )
            analysis_widget["save"].disabled = True
        analysis_widget.add_dropdown(
            "metric", options=["SSIM",], 
            description="Metric for image comparison",
            disabled = True
        )
        analysis_widget.add_button(
            "analyse", description="Run analysis"
        )
        analysis_widget.elements["saving_directory"] = FileChooser(
            self.ouput_directory,
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
        analysis_widget.add_button(
            "save", description="save analysis", disabled=True
        )
        analysis_widget["analyse"].on_click(analyse_sweep)
        analysis_widget["save"].on_click(save_results)
        analysis_widget.show()