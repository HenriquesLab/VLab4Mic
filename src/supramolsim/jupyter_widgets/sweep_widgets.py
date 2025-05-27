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
            self.sweep_gen.set_analysis_parameters("metrics_list", metric_list)
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


    def structure_probe_ui(self, height = '400px'):
        structures = ["9I0K", "1XI5"]
        #structures = ["7R5K", "1XI5"]
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

        def plot_structure2(structure_id, n_atoms=1000, h_rotation=0, v_rotation=0):
            total = list_of_experiments[structure_id].structure.num_assembly_atoms
            if total > n_atoms:
                fraction = n_atoms/total
            else:
                fraction = 1.0
            list_of_experiments[structure_id].structure.show_assembly_atoms(
                assembly_fraction=fraction,
                view_init = [v_rotation,h_rotation,0]

            )
        structure_widget = self.wgen.gen_interactive_dropdown(
            options=structures,
            orientation="vertical",
            routine=plot_structure2,
            n_atoms=["int_slider", [1e4,0,1e5,100]],
            h_rotation=["int_slider", [0,-90,90,1]],
            v_rotation=["int_slider", [0,-90,90,1]],
            height=height)
        # 
        vlabprobes = []
        unspecific_probes = copy.copy(self.probes_per_structure["Mock"])
        list_of_probe_objects = {}
        for structurename, probe_names in  self.probes_per_structure.items():
            if structurename != "Mock":
                vlabprobes = vlabprobes + probe_names
        with io.capture_output() as captured:
            vsample, experiment = experiments.generate_virtual_sample(
                clear_probes=True,
                )
            for probe_name in unspecific_probes:
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
                    list_of_probe_objects[probe_name]["probe_object"] = probe_obj


        def show_probe(probe, n_atoms, h_rotation=0, v_rotation=0):
            if probe == "Linker":
                list_of_probe_objects[probe_name]["probe_object"].plot_emitters()
            else:
                total = list_of_probe_objects[probe]["probe_structure"].num_assembly_atoms
                if total > n_atoms:
                    fraction = n_atoms/total
                else:
                    fraction = 1.0
                list_of_probe_objects[probe]["probe_structure"].plotting_params["assemblyatoms"]["plotalpha"] = 0.3
                #list_of_probe_objects[probe]["probe_structure"].show_assembly_atoms(
                #assembly_fraction=fraction,
                #view_init = [30,degree,0]
                #)
                list_of_probe_objects[probe]["probe_structure"].show_target_labels(
                    with_assembly_atoms = True,
                    assembly_fraction=fraction,
                    view_init = [v_rotation, h_rotation,0],
                    show_axis = False 
                )

        probes_widget_2 = self.wgen.gen_interactive_dropdown(
            options=list(list_of_probe_objects.keys()),
            orientation="vertical", routine=show_probe,
            n_atoms=["int_slider", [100,0,10000,1]],
            h_rotation=["int_slider", [0,-90,90,1]],
            v_rotation=["int_slider", [0,-90,90,1]],
            height=height)
        
        def my_update(new_value, dependant, update_params):
            #print("change")
            #dependant.options = update_params["options"][new_value]
            pass


        left_parameters_linkd = self.wgen.gen_box_linked(
            w1=structure_widget, 
            w2=probes_widget_2, 
            observed=structure_widget.children[0].children[0],
            dependant=probes_widget_2.children[0].children[0],
            update_method = my_update,
            update_params = copy.copy(self.probes_per_structure)
            )
        left_parameters_linkd.layout = widgets.Layout(width='50%',display='inline-flex')

        def calculate_labelled_particle(widget, options, emitter_plotsize, source_plotsize):
            structure = widget.children[0].children[0].children[0].value
            probe_name = widget.children[1].children[0].children[0].value
            probe_params = options[structure]
            #print(structure, probe_name, probe_params)
            with io.capture_output() as captured2:
                #vsample, experiment = .generate_virtual_sample(
                self.my_experiment.remove_probes()
                self.my_experiment.structure_id = structure
                self.my_experiment.add_probe(probe_name, **probe_params)
                self.my_experiment.build(modules=["structure", "particle"])
                fig = plt.figure()
                ax = fig.add_subplot(111, projection="3d")
                self.my_experiment.particle.gen_axis_plot(
                    axis_object=ax,
                    with_sources=True, 
                    axesoff=True,
                    emitter_plotsize=emitter_plotsize,
                    source_plotsize=source_plotsize
                    )
                #print(experiment.particle.emitters)
                plt.close()
                return fig


        static = self.wgen.gen_action_with_options(
            param_widget=left_parameters_linkd, 
            routine=calculate_labelled_particle, 
            emitter_plotsize=["int_slider", [1,0,30,1]], 
            source_plotsize=["int_slider", [1,0,30,1]],
            options=structure_target_suggestion,
            action_name="Generate labelled particle",
            height=height)
        
        #main_widget = self.wgen.gen_box(widget1=left_parameters_linkd, widget2=static)
        main_widget = GridspecLayout(1, 3)
        main_widget[0,0] = left_parameters_linkd.children[0]
        main_widget[0,1] = left_parameters_linkd.children[1]
        main_widget[0,2] = static
        #main_widget.layout = widgets.Layout(width='100%',display='inline-flex')
        return main_widget
    


    def vsample_vmicroscope_ui(self, mode = "default", height = "500px"):
        grid = GridspecLayout(5, 3, height=height)
        preview_exp = copy.deepcopy(self.my_experiment)
        def create_field(field_config = None,
                         nparticles = 1,
                         random_pl = None,
                         min_distance = None,
                         random_orientations = None,
                         **kwargs):
            with io.capture_output() as captured:
                self.my_experiment._build_coordinate_field(
                    keep=True,
                    nparticles=nparticles,
                    random_placing=random_pl,
                    minimal_distance=min_distance,
                    random_orientations=random_orientations,
                )
                plot = self.my_experiment.coordinate_field.show_field(
                    return_fig=True
                    )
            return plot
        
        wgt1 = self.wgen.gen_interactive_dropdown(
                    options=["option1",],
                    orientation="vertical",
                    routine=create_field,
                    nparticles=["int_slider", [1,0,20,1]],
                    height=height
        )
        grid[:2, 0]  = wgt1.children[0]
        grid[2:, 0] = wgt1.children[1]

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
                print(f"Detector pixelsize (nm): {pixelsize_nm}")
                print(f"PSF sd (nm): {psf_sd_metric}")
                print(f"PSF preview (on a 1x1 µm field of view)")
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


        wgt2 = self.wgen.gen_interactive_dropdown(
                    options=modalities_options,
                    orientation="vertical",
                    routine=show_modality,
                    height=height
        )
        grid[:2, 1]  = wgt2.children[0]
        grid[2:, 1] = wgt2.children[1]


        current_acq = dict()

        def preview_acquisition(widget, exp_time, noise):
            field = self.my_experiment.coordinate_field.export_field()
            preview_exp.exported_coordinate_field = field
            preview_exp.objects_created["exported_coordinate_field"] = True
            selected_mod = widget.children[0].children[0].value
            print(f"Preview for {selected_mod}")
            fig = plt.figure()
            ax = fig.add_subplot(111)
            with io.capture_output() as captured:
                preview_exp.update_modality(modality_name=selected_mod,remove=True)
                preview_exp.add_modality(modality_name=selected_mod, save=False)
                preview_exp.set_modality_acq(modality_name=selected_mod, exp_time=exp_time, noise=noise)
                preview_exp.build(modules=["imager",])
                # consider using run_simulation
                timeseries, calibration_beads = (
                    preview_exp.imager.generate_imaging(
                        modality=selected_mod, exp_time=exp_time, noise=noise
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
            exp_time = ["float_slider", [0.01,0,0.05,0.001]],
            noise = ["checkbox", True],
            button1 = ["button", ["Set parameters of preview", button_method]],
            button2 = ["button", ["Clear all modalities", button_method2]],
            options=None,
            action_name="Preview acquisition")
        grid[:2, 2]  = static.children[0]
        grid[2:, 2] = static.children[1]
        return grid