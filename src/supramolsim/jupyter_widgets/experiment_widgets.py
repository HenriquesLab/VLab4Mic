from ..experiments import ExperimentParametrisation
from ..jupyter_widgets import parameter_selection, visualisation
from ezinput import EZInput
import os
import sys
import matplotlib.pyplot as plt
import copy
from IPython.display import display, clear_output
import ipywidgets as widgets
import io
from ipyfilechooser import FileChooser
from pathlib import Path

IN_COLAB = 'google.colab' in sys.modules
if IN_COLAB:
    output_path = "/content/vlab4mic_outputs"
else:
    output_path = Path.home() / "vlab4mic_outputs"


if not os.path.exists(output_path):
    os.makedirs(output_path)


def _bind_widgets(parameters=None, visualisation=None):
    """
    Bind the widgets to the parameters and visualisation objects.
    """
    gui = EZInput("Main_widget")
    for tag, element in parameters.elements.items():
        gui.elements[tag] = element
    if visualisation is not None:
        for tag, element in visualisation.elements.items():
            gui.elements[tag] = element
    return gui

def select_structure_widget(experiment):
    structure_params = parameter_selection.ui_select_structure(experiment)
    view_structure = visualisation.ui_show_structure(experiment)
    select_structure = _bind_widgets(structure_params, view_structure)
    return select_structure

def select_probe_widget(experiment):
    probes_params = parameter_selection.ui_select_probe(experiment)
    view_probes = visualisation.ui_show_labelled_structure(experiment)
    select_probe = _bind_widgets(probes_params, view_probes)
    return select_probe

def select_sample_parameters_widget(experiment):
    sample_params = parameter_selection.ui_select_sample_parameters(experiment)
    view_sample = visualisation.ui_show_virtual_sample(experiment)
    select_sample = _bind_widgets(sample_params, view_sample)
    return select_sample