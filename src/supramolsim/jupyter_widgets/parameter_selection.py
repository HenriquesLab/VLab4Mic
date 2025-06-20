from ezinput import EZInput
import matplotlib.pyplot as plt
import copy
from IPython.display import display, clear_output
import ipywidgets as widgets


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
        update_probe_list()

    def update_probe_list():
        probes_gui["message1"].value = ""
        for probe in experiment.probe_parameters.keys():
            probes_gui["message1"].value += probe + "<br>"

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
    
    return probes_gui   


