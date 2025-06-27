from ezinput import EZInput
import matplotlib.pyplot as plt
import copy
from IPython.display import display, clear_output
import ipywidgets as widgets
import io
from ipyfilechooser import FileChooser
import copy
from supramolsim.utils.visualisation.matplotlib_plots import slider_normalised

def ui_select_structure(experiment):
    gui = EZInput("Select_structure")
    def select_structure(elements):
        elements["label_1"].value = "Current structure selected: Loading..."
        elements["select_structure"].disabled = True
        experiment.structure_id = experiment.structures_info_list[elements["structures"].value]
        experiment.build(modules="structure")
        update_structure_list()
        elements["select_structure"].disabled = False

    def update_structure_list():
        if experiment.structure_id is not None:
            gui["label_1"].value = "Current structure selected: " + experiment.structure_id
        else:
            gui["label_1"].value = "Current structure selected: " + "No structure selected yet." 

    if experiment.structure_id is not None:
        gui.add_label(value="Current structure selected: " + experiment.structure_id)
    else:
        gui.add_label(value="Current structure selected: No structure selected yet")
    gui.add_dropdown("structures", description="Select Structure:", options=experiment.structures_info_list.keys())
    gui.add_label("Note: Time for structure loading varyies depending on the size of the structure")
    gui.add_callback(
        "select_structure",
        select_structure,
        gui.elements,
        description="Select structure",
    )
    
    return gui

def update_widgets_visibility(ezwidget, visibility_dictionary):
    for widgetname in visibility_dictionary.keys():
        if visibility_dictionary[widgetname]:
            ezwidget[widgetname].layout.display = "inline-flex"
        else:
            visibility_dictionary[widgetname].layout.display = "None"   


def ui_select_probe(experiment, **kwargs):
    experiment.remove_probes()
    probes_gui = EZInput(title="Labels")
    visibility_widgets = dict()
    probe_options = []
    if experiment.structure_id in experiment.config_probe_per_structure_names.keys():
        probe_list = experiment.config_probe_per_structure_names[experiment.structure_id]
        probe_options.extend(
            copy.copy(probe_list)
        )
    # add probes with no targets
    probe_options.extend(
            copy.copy(experiment.config_global_probes_names)
        )
    # methods
    def select_probe(values):
        experiment.add_probe(
            probe_name=values["select_probe"].value
        )
        probes_gui["create_particle"].disabled = False
        update_probe_list()

    def update_probe_list():
        probes_gui["message1"].value = ""
        for probe in experiment.probe_parameters.keys():
            probes_gui["message1"].value += probe + "<br>"

    def create_particle(b):
        experiment.build(modules=["particle"])
        probes_gui["add_probe"].disabled = True
        probes_gui["create_particle"].disabled = True
        if experiment.generators_status("particle"):
            probes_gui["message2"].value = "Particle created successfully!"
        else:
            probes_gui["message2"].value = "Particle creation failed. Check the logs for details."

    # widgets
    ## Feedback labels
    probes_gui.add_label("Seleced probes:")
    probes_gui.add_HTML("message1", "No probes selected yet.", style = dict(font_weight='bold'))
    # pre-built probes
    probes_gui.add_dropdown("select_probe",
                            description="Choose a probe:",
                            options=probe_options)
    probes_gui.add_callback(
        "add_probe",
        select_probe,
        probes_gui.elements,
        description="Select probe",
    )
    probes_gui.add_button("create_particle", 
                          description="Create labelled structure",
                          disabled=True)
    probes_gui.add_HTML("message2", "No particle created yet.", style = dict(font_weight='bold'))
    probes_gui["create_particle"].on_click(create_particle)

    return probes_gui   

def ui_select_sample_parameters(experiment):
    sample_gui = EZInput(title="Sample parameters")
    # Add widgets for sample parameters
    sample_gui.add_label(
        "Current sample parameters selected:"
    )

    sample_gui.add_HTML(
        "message", ""
    )
    def update_message():
        text = ""
        for key, value in experiment.virtualsample_params.items():
            if key in sample_gui.elements.keys():
                text += f"{key}: {value}<br>"
        sample_gui["message"].value = text
    

    sample_gui.add_int_slider(
        "number_of_particles",
        description="Number of particles",
        min=1,
        max=20,
        value=1,
        continuous_update=False,
        style={'description_width': 'initial'}
    )
    sample_gui.add_checkbox(
        "random_orientations",
        description="Randomise orientations",
        value=True
    )
    sample_gui.add_button(
        "select_sample_parameters",
        description="Select sample parameters",
        disabled=False
    )
    def select_virtual_sample_parameters(b):
        experiment.set_virtualsample_params(
            number_of_particles=sample_gui["number_of_particles"].value,
            random_orientations=sample_gui["random_orientations"].value
        )
        experiment.build(modules=["coordinate_field"])
        if experiment.objects_created["imager"]:
            experiment.build(modules=["imager"])
        update_message()

    sample_gui["select_sample_parameters"].on_click(select_virtual_sample_parameters)
    return sample_gui

def ui_select_modality(experiment):
    modalities_default = ["Widefield", "Confocal", "STED", "SMLM", "All"]
    preview_experiment = copy.deepcopy(experiment)
    xy_zoom_in = 0.5
    for mod_names in modalities_default[0:len(modalities_default)-1]:
        preview_experiment.add_modality(modality_name=mod_names, save=True)
    preview_experiment.build(modules=["imager"])
    modality_gui = EZInput(title="Modality selection")
    modality_gui.add_label("Current modalities selected:")
    modality_gui.add_HTML("message", "No modalities selected yet.")
    def update_message():
        text = ""
        for mod_name, params in experiment.imaging_modalities.items():
            text += f"{mod_name}<br>"
        modality_gui["message"].value = text
    def add_modality(b):
        selected_modality = modality_gui["modality"].value
        if selected_modality == "All":
            for mod_names in modalities_default[0:len(modalities_default)-1]:
                experiment.add_modality(modality_name=mod_names, save=True)
        else:   
            experiment.add_modality(
                modality_name=selected_modality
            )
        update_message()
    
    def remove_modality(b):
        selected_modality = modality_gui["modality"].value
        if selected_modality == "All":
            for mod_names in modalities_default[0:len(modalities_default)-1]:
                experiment.update_modality(
                    modality_name=mod_names,
                    remove=True
                )
        if selected_modality in experiment.imaging_modalities:
            experiment.update_modality(
                modality_name=selected_modality,
                remove=True
            )
        else:
            print(f"Modality {selected_modality} not found.")
        update_message()
    
    def select_modalities(b):
        experiment.build(modules=["imager"])
        modality_gui["add_modality"].disabled = True
        modality_gui["remove_modality"].disabled = True
        modality_gui["select_modalities"].disabled = True 

    def update_plot(change):
        mod_name = modality_gui["modality"].value
        if mod_name != "All":
            psf_stack = preview_experiment.imager.get_modality_psf_stack(mod_name)
            psf_shape = psf_stack.shape
            half_xy = int(psf_shape[0] / 2)
            half_z = int(psf_shape[2] / 2)
            psf_stack = psf_stack[
                half_xy - int(half_xy * xy_zoom_in) : half_xy + int(half_xy * xy_zoom_in),
                half_xy - int(half_xy * xy_zoom_in) : half_xy + int(half_xy * xy_zoom_in),
                :]
            dimension_plane = modality_gui["dimension_slice"].value
            if dimension_plane == "YZ plane":
                dimension = 0
            elif dimension_plane == "XZ plane":
                dimension = 1
            elif dimension_plane == "XY plane":
                dimension = 2
            modality_gui["preview_modality"].clear_output()
            with modality_gui["preview_modality"]:
                display(slider_normalised(
                    psf_stack,
                    dimension=dimension,
                    cbar=False,))

    modality_gui.add_dropdown(
        "modality",
        description="Modality",
        options=modalities_default,
        on_change=update_plot,
    )
    modality_gui.add_button(
        "add_modality",
        description="Add modality",
        disabled=False
    )
    modality_gui.add_button(
        "remove_modality",
        description="Remove modality",
        disabled=False
    )
    modality_gui.add_button(
        "select_modalities",
        description="Select and update virtual modalities",
    )
    modality_gui.add_custom_widget(
        "dimension_slice",
        widgets.ToggleButtons,
        options=["XY plane", "XZ plane", "YZ plane"],
        value="XY plane",
        on_change=update_plot,
        style={"description_width": "initial"},
        description="Plane of view: ",
    )
    modality_gui.add_output(
        "preview_modality",
        description="Preview of selected modality",
        style={"description_width": "initial"},
    )
    modality_gui["add_modality"].on_click(add_modality)
    modality_gui["remove_modality"].on_click(remove_modality)
    modality_gui["select_modalities"].on_click(select_modalities)
    update_message()
    update_plot(True)
    return modality_gui

def ui_run_experiment(experiment):
    run_gui = EZInput(title="Run experiment")
    #experiment.build(modules=["imager",])
    def run_simulation(b):
        run_gui["Acquire"].disabled = True
        sav_dir = run_gui["saving_directory"].value
        if sav_dir is not None:
            experiment.output_directory = sav_dir
            save = True
        experiment.experiment_id = run_gui["experiment_name"].value
        experiment.run_simulation(save=save)
        run_gui.save_settings()

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
    run_gui.add_button("Acquire", description="Run Simulation")
    run_gui["Acquire"].on_click(run_simulation)
    return run_gui