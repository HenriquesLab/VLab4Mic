from ezinput import EZInput
import matplotlib.pyplot as plt
import copy
from IPython.display import display, clear_output
import ipywidgets as widgets
import io

def ui_select_structure(experiment):
    gui = EZInput("Select_structure")
    def select_structure(elements):
        elements["button"].disabled = True
        elements["label_2"].value = "Loading structure, this can take a few seconds..."
        experiment.structure_id = elements["structures"].value
        experiment.build(modules="structure")
        elements["label_2"].value = experiment.structure_id

    gui.add_label(value="Current structure selected:")
    gui.add_label(value=experiment.structure_id)
    gui.add_dropdown("structures", description="Select Structure:", options=experiment.config_probe_per_structure_names.keys())
    gui.add_callback(
        "button",
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

    # widgets
    ## Feedback labels
    probes_gui.add_label("Seleced probes:")
    probes_gui.add_HTML("message1", "")
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
        update_message()

    sample_gui["select_sample_parameters"].on_click(select_virtual_sample_parameters)
    return sample_gui