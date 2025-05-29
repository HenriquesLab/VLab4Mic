import supramolsim.utils
from ..sweep_generator import sweep_generator
import ipywidgets as widgets
from .. import experiments
from .widget_generator import widgen
from .widgets_dataclass import jupyter_gui
import os
import supramolsim
from supramolsim.utils.io import yaml_functions
import copy
from ezinput import EZInput
from ipyfilechooser import FileChooser
from IPython.utils import io
from ..generate.labels import construct_label
from ..workflows import probe_model
import matplotlib.pyplot as plt
from ipywidgets import GridspecLayout
from ..utils import data_format
from supramolsim.workflows import create_imaging_system
import numpy as np
import mpl_toolkits.axes_grid1 as axes_grid1


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
        #self.parameters_with_set_values = []
        self.probes_per_structure = copy.copy(self.my_experiment.config_probe_per_structure_names)
        self.probe_parameters = copy.copy(self.my_experiment.config_probe_params)
        self.vlab_probes = copy.copy(self.my_experiment.config_global_probes_names)
        
        #for file in os.listdir(self.config_directories["probes"]):
        #    if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
        #        label_config_path = os.path.join(
        #            self.config_directories["probes"], file
        #        )
        #        label_parmeters = supramolsim.load_yaml(label_config_path)
        #        lablname = os.path.splitext(file)[0]
        #        # self.vlab_probes.append(lablname)
        #        # self.probe_parameters[lablname] = label_parmeters
        #        if "Mock" in label_parmeters["known_targets"]:
        #            self.targetless_probes.append(lablname)
        #            self.probe_parameters[lablname] = label_parmeters
        #        if "Generic" in label_parmeters["known_targets"]:
        #            self.vlab_probes.append(lablname)
        #            self.probe_parameters[lablname] = label_parmeters
        #        else:
        #            #self.vlab_probes.append(lablname)
        #            self.probe_parameters[lablname] = label_parmeters
        #            for known_structures in label_parmeters["known_targets"]:
        #                if known_structures in self.probes_per_structure.keys():
        #                    self.probes_per_structure[known_structures].append(lablname)
        #                else:
        #                    self.probes_per_structure[known_structures] = [
        #                        lablname,
        #                    ]

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
        if self.selected_structure in self.probes_per_structure.keys():
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
            description="probes", options=probes2show
        )
        widget_modules["modalities"] = widgets.SelectMultiple(
            description="modalities", options=self.modalities_default
        )
        # create tabs
        tab_name = list(widget_modules.keys())
        children = [widget_modules[name] for name in tab_name]
        ez_sweep.elements["tabs"] = widgets.HBox(children)

        #ez_sweep.elements["tabs"].children = children
        #ez_sweep.elements["tabs"].titles = tab_name

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

    def add_parameters_values(self):
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
            #self.parameters_with_set_values.append(param_name)
            if self.param_settings[param_group][param_name]["wtype"] != "logical":
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

            self.sweep_gen.set_parameter_values(
                param_group=param_group,
                param_name=param_name,
                values=param_values,
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
            plots = analysis_widget["plots"].value
            param_names_set = self.sweep_gen.parameters_with_set_values
            if len(param_names_set) >= 2:
                self.sweep_gen.set_plot_parameters(
                    "heatmaps", 
                    param1=param_names_set[0], 
                    param2=param_names_set[1])
            if analysis_widget["metric"].value == "All":
                metric_list = ["ssim", "pearson"]
            elif analysis_widget["metric"].value == "SSIM":
                metric_list = ["ssim", ]
            elif analysis_widget["metric"].value == "Pearson":
                metric_list = ["pearson", ]
            self.sweep_gen.set_number_of_repetitions(analysis_widget["reps"].value)
            self.sweep_gen.set_analysis_parameters(metrics_list = metric_list)
            with io.capture_output() as captured:
                if self.sweep_gen.reference_image is None:
                    self.sweep_gen.generate_reference_image()
            with analysis_widget["outputs"]:
                self.sweep_gen.run_analysis(plots=plots, save=False)
            analysis_widget["saving_directory"].disabled = False
            analysis_widget["save"].disabled = False
            analysis_widget["output_name"].disabled = False
        def save_results(b):
            output_name = analysis_widget["output_name"].value
            self.sweep_gen.save_analysis(
                output_directory=self.ouput_directory,
                output_name=output_name
                )
            #analysis_widget["save"].disabled = True
        analysis_widget.elements["reps"] = self.wgen.gen_bound_int(
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
        analysis_widget.elements["outputs"] = widgets.Output()
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


    def structure_probe_ui(self, height = '400px', structures=["9I0K", "1XI5"]):
        main_widget = GridspecLayout(7, 3, height=height)
        params_section = 2
        list_of_experiments = dict()
        structure_target_suggestion = dict()
        with io.capture_output() as captured:
            for struct in structures:
                list_of_experiments[struct] = experiments.ExperimentParametrisation()
                list_of_experiments[struct].structure_id = struct
                list_of_experiments[struct]._build_structure()
                protein_name = None
                sequence = None
                protein_name, _1, site, sequence = (
                    list_of_experiments[struct].structure.get_peptide_motif(position="cterminal")
                )
                structure_target_suggestion[struct] = {}
                structure_target_suggestion[struct]["probe_target_type"] = "Sequence"
                structure_target_suggestion[struct]["probe_target_value"] = sequence
        structure_name = widgets.Dropdown(options=structures)
        n_atoms = widgets.IntSlider(value=1e4, min=0, max=1e5, steps = 100, description="Atoms to display", style = {'description_width': 'initial'}, continuous_update=False)
        h_rotaiton = widgets.IntSlider(value=0, min=-90, max=90, description="Horizontal view",  style = {'description_width': 'initial'}, continuous_update=False)
        v_rotation = widgets.IntSlider(value=0, min=-90, max=90, description="Vertical view",  style = {'description_width': 'initial'}, continuous_update=False)
        structure_params = [structure_name, n_atoms, h_rotaiton, v_rotation]
        structure_output = widgets.Output()
        def on_changes(change):
            structure_output.clear_output()
            with structure_output:
                p1 = structure_name.value
                p2 = n_atoms.value
                p3 = h_rotaiton.value
                p4 = v_rotation.value
                index = structure_params.index(change.owner)
                if index == 0: # structure id
                    p1 = change.new
                elif index == 1:
                    p2 = change.new
                elif index == 2:
                    p3 = change.new
                elif index == 3:
                    p4 = change.new
                total = list_of_experiments[p1].structure.num_assembly_atoms
                if total > p2:
                    fraction = p2/total
                else:
                    fraction = 1.0
                with io.capture_output() as captured:   
                    figure = list_of_experiments[p1].structure.show_assembly_atoms(
                        assembly_fraction=fraction,
                        view_init = [p4,p3,0]
                    )
                plt.close()
                display(figure)
        structure_name.observe(on_changes, names="value")
        n_atoms.observe(on_changes, names="value")
        h_rotaiton.observe(on_changes, names="value")
        v_rotation.observe(on_changes, names="value")
        main_widget[:params_section, 0]  = widgets.VBox(structure_params)
        main_widget[params_section:, 0] = structure_output
        # probes 
        list_of_probe_objects = {}
        with io.capture_output() as captured:
            vsample, experiment = experiments.generate_virtual_sample(
                clear_probes=True,
                )
            for probe_name in self.probe_parameters.keys():
                print(probe_name)
                label_config_path = os.path.join(experiment.configuration_path, "probes", probe_name + ".yaml")
                probe_obj, probe_parameters = construct_label(label_config_path)
                if probe_obj.model:
                    (
                        probe_structure_obj,
                        probe_emitter_sites,
                        anchor_point,
                        direction_point,
                        probe_epitope,
                    ) = probe_model(
                        model=probe_obj.model,
                        binding=probe_obj.binding,
                        conjugation_sites=probe_obj.conjugation,
                        epitope=probe_obj.epitope,
                        config_dir=experiment.configuration_path,
                    )
                    if anchor_point.shape == (3,):
                            print("setting new axis")
                            probe_obj.set_axis(pivot=anchor_point, direction=direction_point)
                    if (
                        probe_epitope["coordinates"] is not None
                        and probe_parameters["as_linker"]
                    ):
                        print("Generating linker from epitope site")
                        # TODO: this decision needs to take into account if there is a 
                        # secondary label for this specific probe
                        probe_obj.set_emitters(probe_epitope["coordinates"])
                    else:
                        probe_obj.set_emitters(probe_emitter_sites)
                    probe_parameters["coordinates"] = probe_obj.gen_labeling_entity()
                    list_of_probe_objects[probe_name] = {}
                    list_of_probe_objects[probe_name]["probe_object"] = probe_obj
                    list_of_probe_objects[probe_name]["probe_structure"] = probe_structure_obj
                else:
                    list_of_probe_objects[probe_name] = {}
                    list_of_probe_objects[probe_name]["probe_structure"] = None
                    list_of_probe_objects[probe_name]["probe_object"] = probe_obj
                    if probe_parameters["target"]["type"] is not None:
                        if probe_parameters["target"]["type"] == "Sequence":
                            text = "This probe targets the sequence: "
                            text = text + probe_parameters["target"]["value"]
                        else:
                            text = "This probe targets a residue: "
                            text = text + probe_parameters["target"]["value"]["residues"]
                        list_of_probe_objects[probe_name]["probe_info_text"] = text

        def show_probe(probe, n_atoms, h_rotation=0, v_rotation=0):
            if probe in list_of_probe_objects.keys():
                if probe == "Linker":
                    list_of_probe_objects[probe]["probe_object"].plot_emitters()
                else:
                    if list_of_probe_objects[probe]["probe_structure"] is not None:
                        total = list_of_probe_objects[probe]["probe_structure"].num_assembly_atoms
                        if total > n_atoms:
                            fraction = n_atoms/total
                        else:
                            fraction = 1.0
                        list_of_probe_objects[probe]["probe_structure"].plotting_params["assemblyatoms"]["plotalpha"] = 0.3
                        list_of_probe_objects[probe]["probe_structure"].show_target_labels(
                            with_assembly_atoms = True,
                            assembly_fraction=fraction,
                            view_init = [v_rotation, h_rotation,0],
                            show_axis = False 
                        )
                    else:
                        fig, ax = plt.subplots(figsize=[4,4])
                        ax.text(0.5, 0.5, list_of_probe_objects[probe]["probe_info_text"], fontsize=14, ha='center')
                        ax.set_axis_off()  # This hides the axes
                        plt.show()
        probes2show = []
        current_structure = structure_name.value
        if current_structure in self.probes_per_structure.keys():
            probes2show.extend(
                copy.copy(self.probes_per_structure[current_structure])
            )
        probes2show.extend(copy.copy(self.vlab_probes))

        probes_widget_2 = self.wgen.gen_interactive_dropdown(
            options=probes2show,
            orientation="vertical", routine=show_probe,
            n_atoms=["int_slider", [100,0,10000,1]],
            h_rotation=["int_slider", [0,-90,90,1]],
            v_rotation=["int_slider", [0,-90,90,1]],
            height=height)
        
        def my_update(change):
            probes2show = []
            if change.new in self.probes_per_structure:
                probe_list = self.probes_per_structure[change.new]
                probes2show.extend(
                    copy.copy(probe_list)
                )
            probes2show.extend(copy.copy(self.vlab_probes))
            probes_widget_2.children[0].children[0].options = probes2show
        
        structure_name.observe(my_update, names="value")
        
        main_widget[:params_section, 1] = probes_widget_2.children[0]
        main_widget[params_section:, 1] = probes_widget_2.children[1]
        ## show particle


        def calculate_labelled_particle(b):
            struct = structure_name.value
            probe_name = probes_widget_2.children[0].children[0].value
            probe_target_type=None
            probe_target_value=None

            if self.probe_parameters[probe_name]["target"]["type"] is None:
                probe_target_type = structure_target_suggestion[struct]["probe_target_type"]
                probe_target_value = structure_target_suggestion[struct]["probe_target_value"]
            particle_output.clear_output()
            with particle_output:
                with io.capture_output() as captured:   
                    #vsample, experiment = .generate_virtual_sample(
                    self.my_experiment.structure_id = struct
                    self.my_experiment.remove_probes()
                    self.my_experiment.add_probe(probe_name,
                        probe_target_type=probe_target_type,
                        probe_target_value=probe_target_value
                        )
                    list_of_experiments[struct].remove_probes()
                    list_of_experiments[struct].add_probe(probe_name,
                        probe_target_type=probe_target_type,
                        probe_target_value=probe_target_value
                        )
                    #list_of_experiments[struct].add_probe(probe_name, **target_probe_params)
                    list_of_experiments[struct].build(modules=["particle",])    
                    figure = show_particle(struct = struct)
                    plt.close()
                display(figure)
            
        def update_plot(change):
            struct = structure_name.value
            particle_output.clear_output()
            with particle_output:
                psize1 = emitter_plotsize.value
                psize2 = epitope_plotsize.value
                hview = particle_h_rotation.value
                vview = particle_v_rotation.value
                figure = show_particle(
                    struct = struct,
                    emitter_plotsize=psize1,
                    source_plotsize=psize2,
                    hview=hview,
                    vview=vview
                    )
                plt.close()
                display(figure)
        
        def show_particle(struct= None,
                          emitter_plotsize = 1, 
                          source_plotsize = 1, 
                          hview=0,
                          vview=0):
            particle_output.clear_output()
            with io.capture_output() as captured:   
                fig = plt.figure()
                ax = fig.add_subplot(111, projection="3d")
                list_of_experiments[struct].particle.gen_axis_plot(
                            axis_object=ax,
                            with_sources=True, 
                            axesoff=True,
                            emitter_plotsize=emitter_plotsize,
                            source_plotsize=source_plotsize,
                            view_init=[vview, hview, 0]
                            )
                plt.close()
                return fig


        def select_model_action(b):
            self.my_experiment.structure = copy.deepcopy(list_of_experiments[self.my_experiment.structure_id].structure)
            self.my_experiment.objects_created["structure"] = True
            self.my_experiment.build(modules=["particle",])
        
        emitter_plotsize = widgets.IntSlider(value=1, min=1, max=24, description="Emitters size",  style = {'description_width': 'initial'}, continuous_update=False)
        epitope_plotsize = widgets.IntSlider(value=1, min=1, max=24, description="Epitope size",  style = {'description_width': 'initial'}, continuous_update=False)
        particle_h_rotation = widgets.IntSlider(value=0, min=-90, max=90, description="Horizontal view",  style = {'description_width': 'initial'}, continuous_update=False)
        particle_v_rotation = widgets.IntSlider(value=0, min=-90, max=90, description="Vertical view",  style = {'description_width': 'initial'}, continuous_update=False)
        particle_output = widgets.Output()
        preview_button = widgets.Button(description = "Preview labelling")
        set_button = widgets.Button(description = "Use this model")
        preview_button.on_click(calculate_labelled_particle)
        set_button.on_click(select_model_action)
        #
        emitter_plotsize.observe(update_plot, names="value")
        epitope_plotsize.observe(update_plot, names="value")
        particle_h_rotation.observe(update_plot, names="value")
        particle_v_rotation.observe(update_plot, names="value")
        buttons_widget = widgets.HBox([preview_button, set_button])
        main_widget[:params_section, 2]  = widgets.VBox([buttons_widget, emitter_plotsize, epitope_plotsize, particle_h_rotation , particle_v_rotation])
        main_widget[params_section:, 2] = particle_output


        #static = self.wgen.gen_action_with_options(
        #    param_widget=left_parameters_linkd, 
        #    routine=calculate_labelled_particle, 
        #    emitter_plotsize=["int_slider", [1,0,30,1]], 
        #    source_plotsize=["int_slider", [1,0,30,1]],
        #    select_model = ["button", ["Use this model", select_model_action]],
        #    options=structure_target_suggestion,
        #    action_name="Preview labelled particle",
        #    height=height)
        
        #main_widget = self.wgen.gen_box(widget1=left_parameters_linkd, widget2=static)





        #main_widget[0,0] = left_parameters_linkd.children[0]
        #main_widget[0,1] = left_parameters_linkd.children[1]
        #main_widget[0,2] = static
        #main_widget.layout = widgets.Layout(width='100%',display='inline-flex')
        structure_name.value = structures[1]
        return main_widget
    


    def vsample_vmicroscope_ui(self, mode = "default", height = "500px"):
        grid = GridspecLayout(7, 3, height=height)
        preview_exp = copy.deepcopy(self.my_experiment)
        with io.capture_output() as captured:
            self.my_experiment._build_coordinate_field(
                        keep=True,
                        nparticles=1
                    )
        nparticles = widgets.IntSlider(value=1, min=1, max=20, description="Number of particles",  style = {'description_width': 'initial'},continuous_update=False)
        angle_view = widgets.IntSlider(value=20, min=-90, max=90, description="Angle view",  style = {'description_width': 'initial'},continuous_update=False)
        vsample_params = [nparticles, angle_view]
        vsample_output = widgets.Output()
        def on_changes(change):
            vsample_output.clear_output()
            with vsample_output:
                index = vsample_params.index(change.owner)
                if index == 1:
                    with io.capture_output() as captured:
                        plot = self.my_experiment.coordinate_field.show_field(
                            return_fig=True,
                            view_init=[change.new,0,0]
                            )
                        plt.close()
                    display(plot)
                elif index == 0:
                    with io.capture_output() as captured:
                        self.my_experiment._build_coordinate_field(
                            keep=True,
                            nparticles=change.new
                        )
                        plot = self.my_experiment.coordinate_field.show_field(
                            return_fig=True,
                            view_init=[angle_view.value,0,0]
                        )
                        plt.close()
                    display(plot)
        nparticles.observe(on_changes, names="value")
        angle_view.observe(on_changes, names="value")
        grid[:2, 0]  = widgets.VBox(vsample_params)
        grid[2:, 0] = vsample_output
        angle_view.value = 22
        # modalities
        modalities_options = []
        if mode == "default":
            modalities_list = self.modalities_default
        else:
            modalities_list = self.my_experiment.local_modalities_names
        modalities_options = copy.copy(modalities_list)
        modalities_options.append("All")
        modality_info = {}
        for mod in modalities_list:
            mod_info = data_format.configuration_format.compile_modality_parameters(
                mod, self.config_directories["base"]
            )
            modality_info[mod] = mod_info
        with io.capture_output() as captured:
            temp_imager, tmp_modality_parameters = create_imaging_system(
                modalities_id_list=modalities_list,
                config_dir=self.config_directories["base"],
            )


        def show_modality(modality_name):
            if modality_name != "All":
                pixelsize = modality_info[modality_name]["detector"]["pixelsize"]
                pixelsize_nm = pixelsize * 1000
                psf_sd = np.array(
                    modality_info[modality_name]["psf_params"]["std_devs"]
                )
                psf_voxel = np.array(
                    modality_info[modality_name]["psf_params"]["voxelsize"]
                )
                psf_sd_metric = np.multiply(psf_voxel, psf_sd)
                fig, axs = plt.subplots()
                modality_preview = temp_imager.modalities[modality_name]["psf"][
                    "psf_stack"
                ]
                psf_shapes = modality_preview.shape
                stack_max = np.max(modality_preview)
                axs.imshow(
                    modality_preview[:, :, int(psf_shapes[2] / 2)],
                    cmap="gray",
                    interpolation="none",
                    vmin=0,
                    vmax=stack_max,
                )
                axs.set_xticks([])
                axs.set_yticks([])
                s1 = "Detector pixelsize (nm): " + str(pixelsize_nm)
                s2 = "PSF sd (nm): " + str(psf_sd_metric)
                s3 = "PSF preview (on a 1x1 µm field of view)"
                axs.text(0.05, 0.1, s1, transform=axs.transAxes, size = 10, color = "w")
                axs.text(0.05, 0.15, s2, transform=axs.transAxes, size = 10, color = "w")
                axs.text(0.05, 0.2, s3, transform=axs.transAxes, size = 10, color = "w")


        wgt2 = self.wgen.gen_interactive_dropdown(
                    options=modalities_options,
                    orientation="vertical",
                    routine=show_modality,
                    height=height
        )
        grid[:2, 1]  = wgt2.children[0]
        grid[2:, 1] = wgt2.children[1]


        current_acq = dict()

        def preview_acquisition(widget, exposure_time, noise):
            field = self.my_experiment.coordinate_field.export_field()
            preview_exp.exported_coordinate_field = field
            preview_exp.objects_created["exported_coordinate_field"] = True
            selected_mod = widget.children[0].children[0].value
            fig = plt.figure()
            ax = fig.add_subplot(111)
            with io.capture_output() as captured:
                preview_exp.update_modality(modality_name=selected_mod,remove=True)
                preview_exp.add_modality(modality_name=selected_mod, save=False)
                preview_exp.set_modality_acq(modality_name=selected_mod, exp_time=exposure_time, noise=noise)
                preview_exp.build(modules=["imager",])
                # consider using run_simulation
                timeseries, calibration_beads = (
                    preview_exp.imager.generate_imaging(
                        modality=selected_mod, exp_time=exposure_time, noise=noise
                    )
                )
                current_acq = preview_exp.selected_mods[selected_mod]
            min_val = np.min(timeseries[0])
            max_val = np.max(timeseries[0])
            preview_image=ax.imshow(
                timeseries[0],
                cmap="gray",
                interpolation="none",
                vmin=min_val,
                vmax=max_val,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            #ax.set_title("preview channel:" + single_channel)
            #ax.cbar_axes.colorbar(preview_image)
            # grid[i].set_visible(False)
            plt.close()
            return fig

        def button_method(b):
            selected_mod = list(preview_exp.imaging_modalities.keys())[0]
            mod_acq =copy.deepcopy(preview_exp.selected_mods[selected_mod])
            self.my_experiment.add_modality(modality_name=selected_mod, save=True)
            self.my_experiment.set_modality_acq(modality_name=selected_mod, **mod_acq)

        def button_method2(b):
            modalities_set = list(self.my_experiment.imaging_modalities.keys())
            for mod in modalities_set:
                self.my_experiment.update_modality(modality_name=mod, remove=True)


        static = self.wgen.gen_action_with_options(
            param_widget=wgt2, 
            routine=preview_acquisition,
            exposure_time = ["float_slider", [0.01,0,0.05,0.001]],
            noise = ["checkbox", True],
            button1 = ["button", ["Set parameters of preview", button_method]],
            button2 = ["button", ["Clear all modalities", button_method2]],
            options=None,
            action_name="Preview acquisition")
        grid[:2, 2]  = static.children[0]
        grid[2:, 2] = static.children[1]
        return grid