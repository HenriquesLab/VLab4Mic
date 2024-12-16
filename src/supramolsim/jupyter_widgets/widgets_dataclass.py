from dataclasses import dataclass, field, fields
from typing import List, Dict
from ..experiments import ExperimentParametrisation
import easy_gui_jupyter
import os
from ..utils import data_format
from ..utils.io.yaml_functions import load_yaml
import matplotlib.pyplot as plt
from IPython.display import clear_output
from ..generate.molecular_structure import build_structure_cif
import supramolsim
import ipywidgets as widgets
from IPython.utils import io
from supramolsim.workflows import create_imaging_system
import numpy as np


@dataclass
class jupyter_gui:
    my_experiment = ExperimentParametrisation()

    def __post_init__(self):
        self.demo_structures = []
        # get available structure IDs
        self.structures_info_list = dict()
        structure_dir = os.path.join(
            self.my_experiment.configuration_path, "structures"
        )
        fluorophores_dir = os.path.join(
            self.my_experiment.configuration_path, "fluorophores"
        )
        labels_dir = os.path.join(self.my_experiment.configuration_path, "labels")
        modalities_dir = os.path.join(
            self.my_experiment.configuration_path, "modalities"
        )
        for file in os.listdir(structure_dir):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                structure_params = load_yaml(os.path.join(structure_dir, file))
                struct_id = structure_params["id"]
                strict_title = structure_params["title"]
                id_title = struct_id + ": " + strict_title
                self.structures_info_list[id_title] = struct_id
                self.demo_structures.append(id_title)
        self.config_directories = dict(
            structure=structure_dir,
            fluorophores=fluorophores_dir,
            labels=labels_dir,
            modalities=modalities_dir,
            base=self.my_experiment.configuration_path,
        )

    def structure_gui(self):
        structure_gui = easy_gui_jupyter.EasyGUI("Structure")
        active_widgets = dict()
        nostructure = True

        def select_structure(b):
            structure_id = self.structures_info_list[
                structure_gui["struct_dropdown"].value
            ]
            self.my_experiment.structure_id = structure_id
            self.my_experiment._build_structure()
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False

        def view_structure(b):
            plt.clf()
            clear_output()
            active_widgets["Elevation"] = None
            active_widgets["Azimuthal"] = None
            active_widgets["Roll"] = None
            for widgetname in active_widgets:
                display(structure_gui[widgetname])
            fraction = structure_gui["Fraction"].value
            el = structure_gui["Elevation"].value
            az = structure_gui["Azimuthal"].value
            rll = structure_gui["Roll"].value
            axes_off = structure_gui["hidegrid"].value
            self.my_experiment.structure.show_assembly_atoms(
                assembly_fraction=fraction, view_init=[el, az, rll], axesoff=axes_off
            )

        def upload_file(b):
            filepath = structure_gui["File"].selected
            filename = structure_gui["File"].selected_filename
            structure = build_structure_cif(
                cif_file=filepath, struct_title=filename, cif_id=filename
            )
            structure_param = data_format.structural_format.struct_params_format()
            print("Structure Loaded!")
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False
            structure_id = filename.split(".")[0]

        def activate_demos(b):
            plt.clf()
            clear_output()
            active_widgets["struct_dropdown"] = None
            active_widgets["Select"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            active_widgets["hidegrid"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        def activate_upload(b):
            plt.clf()
            clear_output()
            active_widgets["File"] = None
            active_widgets["Upload"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            active_widgets["hidegrid"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        structure_gui.add_button("Demos", description="Load from demo structures")
        structure_gui.add_button("Fileupload", description="Upload CIF file")
        structure_gui.add_label("Select structure model to load")
        structure_gui.add_dropdown("struct_dropdown", options=self.demo_structures)
        structure_gui.add_button("Select", description="Load structure")
        structure_gui.add_file_upload(
            "File", description="Select from file", accept="*.cif", save_settings=False
        )
        structure_gui.add_button("Upload", description="Load from file")
        structure_gui.add_button(
            "View", description="Show structure", disabled=nostructure
        )
        structure_gui.add_float_slider(
            "Fraction",
            value=0.01,
            min=0,
            max=1,
            step=0.001,
            description="Fraction of total atoms to show",
            disabled=nostructure,
            readout_format=".3f",
        )

        structure_gui.add_int_slider(
            "Elevation",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Elevation",
            disabled=nostructure,
        )
        structure_gui.add_int_slider(
            "Azimuthal",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Azimuthal",
            disabled=nostructure,
        )
        structure_gui.add_int_slider(
            "Roll",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Roll",
            disabled=nostructure,
        )
        structure_gui.add_checkbox(
            "hidegrid", description="Hide plot grids", value=True, disabled=nostructure
        )
        structure_gui["Select"].on_click(select_structure)
        structure_gui["View"].on_click(view_structure)
        structure_gui["Upload"].on_click(upload_file)
        structure_gui["Demos"].on_click(activate_demos)
        structure_gui["Fileupload"].on_click(activate_upload)
        display(structure_gui["Demos"], structure_gui["Fileupload"])

    def structural_model_gui(self):
        # ensure that structure has no labels associated
        labels_gui = easy_gui_jupyter.EasyGUI("Labels")
        particle_created = False
        current_labels = dict()
        generic_labels = []
        structure_labels = []
        fluorophores_list = []
        structure_id = self.my_experiment.structure_id
        for fluoid in os.listdir(self.config_directories["fluorophores"]):
            if os.path.splitext(fluoid)[-1] == ".yaml" and "_template" not in fluoid:
                fluorophores_list.append(os.path.splitext(fluoid)[0])
        for file in os.listdir(self.config_directories["labels"]):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                lablname = os.path.splitext(file)[0]
                if lablname.split("_")[0] == "Generic":
                    generic_labels.append(lablname)
                elif lablname.split("_")[0] == structure_id:
                    structure_labels.append(lablname)

        def build_label(b):
            label_id = labels_gui["label_dropdown"].value
            if label_id == "<None>":
                print("Invalid label, no label added")
            else:
                fluorophore_id = labels_gui["fluo_dropdown"].value
                lab_eff = labels_gui["Labelling_efficiency"].value
                tmp_label = data_format.structural_format.label_builder_format(
                    label_id, fluorophore_id, lab_eff
                )
                unique_name = label_id + "_conjugated_" + fluorophore_id
                if unique_name in current_labels.keys():
                    print("label already exist")
                else:
                    current_labels[unique_name] = tmp_label
                    print(f"label added: {unique_name}")

        def build_generic_label(b):
            label_id = labels_gui["generic_label_dropdown"].value
            fluorophore_id = labels_gui["generc_fluo_dropdown"].value
            lab_eff = labels_gui["generic_Labelling_efficiency"].value
            tmp_label = data_format.structural_format.label_builder_format(
                label_id, fluorophore_id, lab_eff
            )
            unique_name = label_id + "_conjugated_" + fluorophore_id
            if unique_name in current_labels.keys():
                print("label already exist")
            else:
                current_labels[unique_name] = tmp_label
                print(f"label added: {unique_name}")

        def clear(b):
            current_labels.clear()

        def show(b):
            for lab in current_labels.keys():
                print(lab)

        def label_struct(b):
            particle_created = True
            labels_list = []
            if len(current_labels.keys()) > 0:
                self.nlabels = len(current_labels)
                for keys, values in current_labels.items():
                    labels_list.append(values)
                # print(labels_list)
                # self.my_experiment._build_particle(keep=True)
                # using this way since Experiment method assumes only one label
                particle = supramolsim.particle_from_structure(
                    self.my_experiment.structure,
                    labels_list,
                    self.my_experiment.configuration_path,
                )
                self.my_experiment.particle = particle
                self.my_experiment.objects_created["particle"] = True
                print("Structure has been labelled")
            else:
                print("No label has been added")

        labels_gui.add_label("Structure specific labels")
        if self.my_experiment.structure is not None:
            labels_gui.add_dropdown("label_dropdown", options=structure_labels)
            labels_gui.add_dropdown("fluo_dropdown", options=fluorophores_list)
            labels_gui.add_float_slider(
                "Labelling_efficiency",
                value=1,
                min=0,
                max=1,
                step=0.01,
                description="Labelling efficiency",
            )
            labels_gui.add_button(
                "Add", description="Add specific label", disabled=(not structure_labels)
            )
            labels_gui["Add"].on_click(build_label)

        labels_gui.add_label("Generic labels")
        labels_gui.add_dropdown("generic_label_dropdown", options=generic_labels)
        labels_gui.add_dropdown("generc_fluo_dropdown", options=fluorophores_list)
        labels_gui.add_float_slider(
            "generic_Labelling_efficiency",
            value=1,
            min=0,
            max=1,
            step=0.01,
            description="Labelling efficiency",
        )
        labels_gui.add_button("Add_generic", description="Add generic label")
        labels_gui["Add_generic"].on_click(build_generic_label)

        labels_gui.add_button("Clear", description="Clear Labels")
        labels_gui.add_button("Show", description="Display current labels")
        labels_gui.add_label(
            "After adding labels, create a labelled model of your structure"
        )
        labels_gui.add_button("Label", description="Label structure")
        labels_gui["Clear"].on_click(clear)
        labels_gui["Show"].on_click(show)
        labels_gui["Label"].on_click(label_struct)
        if self.my_experiment.structure is None:
            print("No structure has been loaded")
        else:
            self.my_experiment.structure._clear_labels()
            labels_gui.show()

    def refine_model_gui(self):
        structural_model_gui = easy_gui_jupyter.EasyGUI("StructuralModel")

        def show_model(b):
            plt.clf()
            clear_output()
            emitter_plotsize = structural_model_gui["emitterplotsize"].value
            source_size = structural_model_gui["sourceplotsize"].value
            structural_model_gui.show()
            if self.my_experiment.particle:
                particle = self.my_experiment.particle
                fig, axs = plt.subplots(1, 3, subplot_kw={"projection": "3d"})
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[0, 0, 0],
                    axis_object=axs[0],
                    emitter_plotsize=emitter_plotsize,
                )
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[30, 0, 0],
                    axis_object=axs[1],
                    emitter_plotsize=emitter_plotsize,
                )
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[90, 0, 0],
                    axis_object=axs[2],
                    emitter_plotsize=emitter_plotsize,
                )
                plt.subplots_adjust(wspace=0.5)
                plt.show()
            else:
                print(
                    "You have not created a labelled structure. "
                    "Make sure you select 'Label structure' button on previous cell"
                )

        def add_defects(b):
            if self.nlabels == 1:
                self.my_experiment.particle.add_defects(
                    eps1=structural_model_gui["eps1"].value,
                    xmer_neigh_distance=structural_model_gui[
                        "xmer_neigh_distance"
                    ].value,
                    deg_dissasembly=structural_model_gui["Defect"].value,
                )
                message = "Defects added"
            else:
                message = (
                    "Defect modelling is currently unsupported for more than one label"
                )
            structural_model_gui.save_settings()
            show_model(b)
            print(message)

        def relabel(b):
            self.my_experiment.particle.generate_instance()
            show_model(b)

        structural_model_gui.add_button("Show", description="Show current model")
        structural_model_gui.add_label("Visualisation parameters")
        structural_model_gui.add_float_slider(
            "emitterplotsize",
            value=24,
            min=0,
            max=50,
            step=1,
            description="Emitter size",
        )
        structural_model_gui.add_float_slider(
            "sourceplotsize", value=1, min=0, max=50, step=1, description="Target size"
        )
        structural_model_gui.add_checkbox(
            "WTarget", description="With target site", value=True
        )
        structural_model_gui.add_checkbox("Axes", description="Hide Axes", value=True)
        structural_model_gui.add_button(
            "Relabel", description="Recalculate labelled particle and show"
        )
        structural_model_gui.add_label("Model defects parameters (optional):")
        structural_model_gui._widgets["eps1"] = widgets.BoundedIntText(
            value=300,
            min=0,
            max=100000,
            description="Short distance cluster",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui._widgets["xmer_neigh_distance"] = widgets.BoundedIntText(
            value=600,
            min=0,
            max=100000,
            description="Long distance cluster",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui._widgets["Defect"] = widgets.BoundedFloatText(
            value=0.5,
            min=0,
            max=1,
            description="percentage of defect",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui.add_button(
            "Defects", description="Model defects and show model"
        )
        structural_model_gui["Show"].on_click(show_model)
        structural_model_gui["Defects"].on_click(add_defects)
        structural_model_gui["Relabel"].on_click(relabel)
        if self.my_experiment.particle:
            structural_model_gui.show()
        else:
            print("No particle has been created")

    def create_field(self):
        field_gui = easy_gui_jupyter.EasyGUI("field")

        def createmin(b):
            random_pl = True
            use = field_gui["use_particle"].value
            self.my_experiment._build_coordinate_field(
                use_self_particle=use, keep=True, random_placing=random_pl
            )
            field_gui["show"].disabled = False

        def upload_custom(b):
            display(field_gui["File"], field_gui["load_custom"])

        def create_custom(b):
            path_to_custom_field = field_gui["File"].selected
            use = field_gui["use_particle"].value
            self.my_experiment._build_coordinate_field(
                use_self_particle=use,
                keep=True,
                coordinate_field_path=path_to_custom_field,
            )
            field_gui["show"].disabled = False

        def showfield(b):
            plt.clf()
            self.my_experiment.coordinate_field.show_field(view_init=[90, 0, 0])

        field_gui.add_button("minimal", description="Create minimal field (random)")
        field_gui.add_button("custom", description="Use custom field params")
        field_gui.add_button("show", description="Show field", disabled=True)
        field_gui.add_button("load_custom", description="Load custom field")
        field_gui.add_file_upload(
            "File",
            description="Upload custom field",
            accept="*.yaml",
            save_settings=False,
        )
        if self.my_experiment.particle is None:
            print("There is no particle initalised")
            use_p_disabled = True
            is_particle = False
        else:
            print("Particle model exists")
            use_p_disabled = False
            is_particle = True
        field_gui.add_checkbox(
            "use_particle",
            description="Use existing particle",
            value=is_particle,
            disabled=use_p_disabled,
        )
        field_gui["minimal"].on_click(createmin)
        field_gui["show"].on_click(showfield)
        field_gui["custom"].on_click(upload_custom)
        field_gui["load_custom"].on_click(create_custom)
        # field_gui.show()
        display(
            field_gui["minimal"],
            field_gui["custom"],
            field_gui["show"],
            field_gui["use_particle"],
        )

    def set_image_modalities(self):
        imager_created = False
        imaging_gui = easy_gui_jupyter.EasyGUI("imaging")
        # if configuration_path is None:
        #    print("No configuration path set")
        # else:
        modalities_dir = self.config_directories["modalities"]
        modalities_list = []
        selected_mods = []
        modality_info = {}

        for mods in os.listdir(modalities_dir):
            if os.path.splitext(mods)[-1] == ".yaml" and "_template" not in mods:
                modalities_list.append(os.path.splitext(mods)[0])

        for mod in modalities_list:
            mod_info = data_format.configuration_format.compile_modality_parameters(
                mod, self.config_directories["base"]
            )
            modality_info[mod] = mod_info

        # create mock imager to show a pre-visualisation of PSF and noise model
        with io.capture_output() as captured:
            temp_imager = create_imaging_system(
                modalities_id_list=modalities_list,
                config_dir=self.config_directories["base"],
            )

        def add_mod(b):
            selected_mods.append(imaging_gui["modalities_dropdown"].value)
            print(f"{selected_mods[-1]} modality added")

        def clear(b):
            selected_mods.clear()

        def create_imager(b):
            if exported_field is not None:
                if len(selected_mods) == 0:
                    print("No modalites had been added")
                else:
                    with io.capture_output() as captured:
                        imaging_system = create_imaging_system(
                            exported_field, selected_mods, configuration_path[0]
                        )
                    imager_created = True
                    print("Imaging system created")
                imaging_gui["Show"].disabled = False
            else:
                print("No field info")

        def preview(b):
            def get_info(imaging_gui):
                def preview_info(Modality):
                    # print(Modality)
                    # print(modality_info[Modality])
                    pixelsize = modality_info[Modality]["detector"]["pixelsize"]
                    pixelsize_nm = pixelsize * 1000
                    psf_sd = np.array(modality_info[Modality]["psf_params"]["std_devs"])
                    psf_voxel = np.array(
                        modality_info[Modality]["psf_params"]["voxelsize"]
                    )
                    psf_sd_metric = np.multiply(psf_voxel, psf_sd)
                    print(f"Detector pixelsize (nm): {pixelsize_nm}")
                    print(f"PSF sd (nm): {psf_sd_metric}")
                    # show PSF
                    fig, axs = plt.subplots()
                    modality_preview = temp_imager.modalities[Modality]["psf"][
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

                widgets.interact(
                    preview_info, Modality=imaging_gui["modalities_dropdown"]
                )

            get_info(imaging_gui)

        def showfov(b):
            imaging_system.show_field()

        imaging_gui.add_label("Select modalities")
        imaging_gui.add_dropdown("modalities_dropdown", options=modalities_list)
        imaging_gui.add_button("Add", description="Add modality")
        imaging_gui.add_button("Clear", description="Clear selection")
        if self.my_experiment.generators_status("exported_coordinate_field") is None:
            disable_button = True
        else:
            disable_button = False
        imaging_gui.add_button(
            "Create", description="Create imaging system", disabled=disable_button
        )
        imaging_gui.add_button("Show", description="Show field of view", disabled=True)
        imaging_gui.add_int_slider(
            "PSF_nslice", min=0, max=400, continuous_update=False
        )
        imaging_gui["Add"].on_click(add_mod)
        imaging_gui["Clear"].on_click(clear)
        imaging_gui["Create"].on_click(create_imager)
        imaging_gui["Show"].on_click(showfov)
        preview(True)
        display(
            imaging_gui["Add"],
            imaging_gui["Clear"],
            imaging_gui["Create"],
            imaging_gui["Show"],
        )
