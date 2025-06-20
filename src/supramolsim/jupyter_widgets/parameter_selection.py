from ezinput import EZInput
import matplotlib.pyplot as plt


def ui_select_structure(experiment):
    gui = EZInput("Select_structure")
    def select_structure(elements):
        elements["button"].disabled = True
        elements["label_2"].value = "Loading structure, this can take a few seconds..."
        experiment.structure_id = elements["structures"].value
        experiment.build(modules="structure")
        elements["label_2"].value = experiment.structure_id

    gui.add_label(value="Selected structure:")
    gui.add_label(value="")
    gui.add_dropdown("structures", description="Select Structure:", options=experiment.config_probe_per_structure_names.keys())
    gui.add_callback(
        "button",
        select_structure,
        gui.elements,
        description="Select structure",
    )
    return gui
    