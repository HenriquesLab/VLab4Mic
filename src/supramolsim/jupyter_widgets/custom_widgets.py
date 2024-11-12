import easy_gui_jupyter
import ipywidgets as widgets
import supramolsim
import supramolsim.generate.coordinates_field as field
from ..utils import data_format
from ..workflows import *
import matplotlib.pyplot as plt
from IPython.display import clear_output
import os
from IPython.utils import io
import numpy as np
import mpl_toolkits.axes_grid1 as axes_grid1

global structure, particle, configuration_path, coordinates_field, exported_field
particle = None
structure = None
configuration_path = None
coordinates_field = None
exported_field = None


def set_directory():
    config_gui = easy_gui_jupyter.EasyGUI("Config")

    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    local_dir = os.path.join(pck_dir, "configuration")

    def clear(b):
        global configuration_path
        configuration_path = []
        configuration_path.clear()
        configuration_path.append(local_dir)
        print(configuration_path)

    def select_custom(b):
        global configuration_path
        configuration_path = []
        configuration_path.clear()
        configuration_path.append(config_gui["custom_path"].value)
        config_gui.save_settings()
        print(configuration_path)

    config_gui.add_textarea("custom_path", remember_value=True)
    config_gui.add_button("Set", description="Set configuration path")
    config_gui["Set"].on_click(select_custom)
    config_gui.add_button("demos", description="Use demos")
    config_gui["demos"].on_click(clear)
    config_gui.show()


def select_structure():
    global structure, structure_params, configuration_path, nostructure, active_widgets
    active_widgets = dict()
    nostructure = True
    demo_structures = []
    structure_gui = easy_gui_jupyter.EasyGUI("Structure")
    if configuration_path is None:
        print("No configuration path")
    else:
        print(configuration_path[0])
        generic_labels = []
        # get available structure IDs
        structures_info_list = dict()
        # configuration_path = configuration_path[0]
        if os.path.isdir(configuration_path[0]):
            structure_dir = os.path.join(configuration_path[0], "structures")
            for file in os.listdir(structure_dir):
                if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                    structure_params = supramolsim.load_yaml(
                        os.path.join(structure_dir, file)
                    )
                    struct_id = structure_params["id"]
                    strict_title = structure_params["title"]
                    id_title = struct_id + ": " + strict_title
                    structures_info_list[id_title] = struct_id
                    demo_structures.append(id_title)

        def select_structure(b):
            global structure, structure_param, struct_ready, structure_id
            structure_id = structures_info_list[structure_gui["struct_dropdown"].value]
            # structure_id = structure_gui["struct_dropdown"].value
            structure, structure_param = supramolsim.load_structure(
                structure_id, configuration_path[0]
            )
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False

        def view_structure(b):
            plt.clf()
            clear_output()
            # structure_gui.show()
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
            structure.show_assembly_atoms(
                assembly_fraction=fraction, view_init=[el, az, rll], axesoff=axes_off
            )

        def upload_file(b):
            global structure, structure_param, structure_id
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
            global active_widgets
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
            global active_widgets
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
        structure_gui.add_dropdown("struct_dropdown", options=demo_structures)
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
        # structure_gui.show()


def create_structural_model():
    global \
        labels_list, \
        current_labels, \
        particle, \
        particle_created, \
        structure, \
        structure_param, \
        structure_id
    # ensure that structure has no labels associated
    labels_gui = easy_gui_jupyter.EasyGUI("Labels")
    if configuration_path is None:
        print("No configuration_path has been loaded")
    else:
        particle_created = False
        current_labels = dict()
        generic_labels = []
        structure_labels = []
        fluorophores_list = []
        fluorophores_dir = os.path.join(configuration_path[0], "fluorophores")
        labels_dir = os.path.join(configuration_path[0], "labels")
        for fluoid in os.listdir(fluorophores_dir):
            if os.path.splitext(fluoid)[-1] == ".yaml" and "_template" not in fluoid:
                fluorophores_list.append(os.path.splitext(fluoid)[0])
        for file in os.listdir(labels_dir):
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
            global current_labels
            for lab in current_labels.keys():
                print(lab)

        def label_struct(b):
            global particle, configuration_path, particle_created, nlabels
            particle_created = True
            labels_list = []
            if len(current_labels.keys()) > 0:
                nlabels = len(current_labels)
                for keys, values in current_labels.items():
                    labels_list.append(values)
                # print(labels_list)
                particle = supramolsim.particle_from_structure(
                    structure, labels_list, configuration_path[0]
                )
                print("Structure has been labelled")
            else:
                print("No label has been added")

        labels_gui.add_label("Structure specific labels")
        if structure is not None:
            labels_gui.add_dropdown("label_dropdown", options=structure_param["labels"])
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
        if structure is None:
            print("No structure has been loaded")
        else:
            structure._clear_labels()
            labels_gui.show()


def refine_structural_model():
    global particle, particle_created, plot_exists
    structural_model_gui = easy_gui_jupyter.EasyGUI("StructuralModel")

    def show_model(b):
        global particle, particle_created, plot_exists
        plt.clf()
        clear_output()
        emitter_plotsize = structural_model_gui["emitterplotsize"].value
        source_size = structural_model_gui["sourceplotsize"].value
        structural_model_gui.show()
        if particle_created:
            # fig = plt.figure()
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
        global particle, nlabels
        if nlabels == 1:
            particle.add_defects(
                eps1=structural_model_gui["eps1"].value,
                xmer_neigh_distance=structural_model_gui["xmer_neigh_distance"].value,
                deg_dissasembly=structural_model_gui["Defect"].value,
            )
            # particle.generate_instance()
            # print("Defects added")
            message = "Defects added"
        else:
            # print("Defect modelling is currently unsupported for more than one label")
            message = (
                "Defect modelling is currently unsupported for more than one label"
            )
        structural_model_gui.save_settings()
        show_model(b)
        print(message)

    def relabel(b):
        particle.generate_instance()
        show_model(b)

    structural_model_gui.add_button("Show", description="Show current model")
    structural_model_gui.add_label("Visualisation parameters")
    structural_model_gui.add_float_slider(
        "emitterplotsize", value=24, min=0, max=50, step=1, description="Emitter size"
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
    if particle:
        structural_model_gui.show()
    else:
        print("No particle has been created")


def create_field():
    global particle_field, exported_field, particle, field_gui, coordinates_field
    field_gui = easy_gui_jupyter.EasyGUI("field")

    def createmin(b):
        global coordinates_field, exported_field
        random_pl = True
        coordinates_field = field.create_min_field(random_placing=random_pl)
        exported_field = coordinates_field.export_field()
        field_gui["show"].disabled = False

    def showfield(b):
        plt.clf()
        global coordinates_field
        coordinates_field.show_field(view_init=[90, 0, 0])

    field_gui.add_button("minimal", description="Create minimal field (random)")
    field_gui.add_button("show", description="Show field", disabled=True)
    field_gui["minimal"].on_click(createmin)
    field_gui["show"].on_click(showfield)

    if particle is None:
        print("There is no particle initalised")
        field_gui.show()
    else:
        print("Using existing structural model")
        field_gui["show"].disabled = False
        display(field_gui["show"])
        exported_field, coordinates_field = field_from_particle(particle)


def set_image_modalities():
    global \
        selected_mods, \
        imager_created, \
        modality_info, \
        temp_imager, \
        configuration_path, \
        exported_field
    imager_created = False
    imaging_gui = easy_gui_jupyter.EasyGUI("imaging")
    if configuration_path is None:
        print("No configuration path set")
    else:
        modalities_dir = os.path.join(configuration_path[0], "modalities")
        modalities_list = []
        selected_mods = []
        modality_info = {}

        for mods in os.listdir(modalities_dir):
            if os.path.splitext(mods)[-1] == ".yaml" and "_template" not in mods:
                modalities_list.append(os.path.splitext(mods)[0])

        for mod in modalities_list:
            mod_info = data_format.configuration_format.compile_modality_parameters(
                mod, configuration_path[0]
            )
            modality_info[mod] = mod_info

        # create mock imager to show a pre-visualisation of PSF and noise model
        with io.capture_output() as captured:
            temp_imager = create_imaging_system(
                modalities_id_list=modalities_list, config_dir=configuration_path[0]
            )

        def add_mod(b):
            global selected_mods
            selected_mods.append(imaging_gui["modalities_dropdown"].value)
            print(f"{selected_mods[-1]} modality added")

        def clear(b):
            global selected_mods
            selected_mods.clear()

        def create_imager(b):
            global imager_created, selected_mods, exported_field, imaging_system
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
                    global modality_info
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

            global imaging_system
            get_info(imaging_gui)

        def showfov(b):
            global imaging_system
            imaging_system.show_field()

        imaging_gui.add_label("Select modalities")
        imaging_gui.add_dropdown("modalities_dropdown", options=modalities_list)
        imaging_gui.add_button("Add", description="Add modality")
        imaging_gui.add_button("Clear", description="Clear selection")
        if exported_field is None:
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


def set_acq_params():
    global selected_mods, acq_params_per_mod, nchannels, imager_channels
    acquisition_gui = easy_gui_jupyter.EasyGUI("acquisition_params")
    acq_params_per_mod = dict()
    imager_channels = []
    anymod = list(imaging_system.modalities.keys())[0]
    for chann in imaging_system.modalities[anymod]["filters"].keys():
        print(chann)
        imager_channels.append(chann)
    nchannels = len(imager_channels)

    def set_params(b):
        mod_id = acquisition_gui["modalities_dropdown"].value
        exp_time = acquisition_gui["Exposure"].value
        noise = acquisition_gui["Noise"].value
        save = True
        nframes = acquisition_gui["Frames"].value
        if acquisition_gui["Channels"].value:
            global imaging_system
            channels = []
            for chann in imaging_system.modalities[mod_id]["filters"].keys():
                channels.append(chann)
            print(f"using all channels: {channels}")
        else:
            channels = [
                "ch0",
            ]
        acq_params_per_mod[mod_id] = format_modality_acquisition_params(
            exp_time=exp_time,
            noise=noise,
            save=save,
            nframes=nframes,
            channels=channels,
        )
        print(f"Acquisition parameters added for {mod_id}")
        print(acq_params_per_mod[mod_id])
        acquisition_gui.save_settings()

    def preview_mod(b):
        preview_image = None

        def get_preview(imaging_system, acq_gui):
            global nchannels, channels

            def preview_exposure(Modality, Exposure, Noise):
                global preview_image
                fig = plt.figure()
                grid = axes_grid1.AxesGrid(
                    fig,
                    111,
                    nrows_ncols=(1, nchannels),
                    axes_pad=1,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="10%",
                    cbar_pad="20%",
                )
                i = 0
                for single_channel in imager_channels:
                    single_mod_acq_params = dict(
                        exp_time=Exposure,
                        noise=Noise,
                        save=False,
                        nframes=1,
                        channel=single_channel,
                    )
                    with io.capture_output() as captured:
                        timeseries, calibration_beads = imaging_system.generate_imaging(
                            modality=Modality, **single_mod_acq_params
                        )
                        min_val = np.min(timeseries[0])
                        max_val = np.max(timeseries[0])
                    preview_image = grid[i].imshow(
                        timeseries[0],
                        cmap="gray",
                        interpolation="none",
                        vmin=min_val,
                        vmax=max_val,
                    )
                    grid[i].set_title("preview: " + single_channel)
                    grid.cbar_axes[i].colorbar(preview_image)
                    i = i + 1

            widgets.interact(
                preview_exposure,
                Modality=acq_gui["modalities_dropdown"],
                Exposure=acq_gui["Exposure"],
                Noise=acq_gui["Noise"],
            )

        global imaging_system
        # global selected_mods
        get_preview(imaging_system, acquisition_gui)

    def clear(b):
        acq_params_per_mod.clear()
        print("Acquisition parameters cleared")
        acquisition_gui.save_settings()

    acquisition_gui.add_label("Set acquisition parameters")
    acquisition_gui.add_dropdown("modalities_dropdown", options=selected_mods)
    acquisition_gui.add_checkbox("Noise", description="Use Noise", value=True)
    acquisition_gui.add_checkbox("Channels", description="Use all channels", value=True)
    ## bounded int Text
    acquisition_gui._widgets["Frames"] = widgets.BoundedIntText(
        value=1,
        min=1,
        max=100000,
        description="Frames (not used for preview)",
        layout=acquisition_gui._layout,
        style=acquisition_gui._style,
        remember_value=True,
    )
    acquisition_gui._widgets["Exposure"] = widgets.BoundedFloatText(
        value=0.005,
        min=0.000000,
        step=0.0001,
        description="exposure (sec)",
        layout=acquisition_gui._layout,
        style=acquisition_gui._style,
        remember_value=True,
    )
    acquisition_gui.add_button("Set", description="Add params")
    acquisition_gui.add_button("Preview", description="Preview (Expermiental feature)")
    acquisition_gui.add_button("Clear", description="Clear params")
    acquisition_gui["Preview"].on_click(preview_mod)
    acquisition_gui["Set"].on_click(set_params)
    acquisition_gui["Clear"].on_click(clear)
    display(acquisition_gui["Set"], acquisition_gui["Frames"])
    preview_mod(True)


def run_simulation():
    experiment_gui = easy_gui_jupyter.EasyGUI("experiment")

    def run_simulation(b):
        sav_dir = experiment_gui["saving_directory"].value
        exp_name = experiment_gui["experiment_name"].value
        print(len(acq_params_per_mod.keys()))
        if len(acq_params_per_mod.keys()) == 0:
            generate_multi_imaging_modalities(
                image_generator=imaging_system,
                experiment_name=exp_name,
                savingdir=sav_dir,
                acquisition_param=None,
            )
        else:
            generate_multi_imaging_modalities(
                image_generator=imaging_system,
                experiment_name=exp_name,
                savingdir=sav_dir,
                acquisition_param=acq_params_per_mod,
            )
        experiment_gui.save_settings()

    experiment_gui.add_label("Set experiment name")
    experiment_gui.add_textarea(
        "experiment_name", value="Exp_name", remember_value=True
    )
    experiment_gui.add_label("Set saving directory")
    experiment_gui.add_textarea("saving_directory", remember_value=True)
    experiment_gui.add_button("Acquire", description="Run Simulation")
    experiment_gui["Acquire"].on_click(run_simulation)
    experiment_gui.show()
