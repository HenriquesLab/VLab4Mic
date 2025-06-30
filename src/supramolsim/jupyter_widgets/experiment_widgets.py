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
    if parameters is not None:
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

def select_modalities_widget(experiment):
    modalities_params = parameter_selection.ui_select_modality(experiment)
    select_modalities = _bind_widgets(modalities_params)
    return select_modalities


def select_acquisition_parameters_widget(experiment):
    """
    Create a widget to select acquisition parameters.
    """
    view_acq_params = visualisation.ui_set_acq_params(experiment)
    select_acq_params = _bind_widgets(visualisation = view_acq_params)
    return select_acq_params

def run_experiment_widget(experiment):
    """
    Create a widget to run the experiment.
    """
    run_experiment = parameter_selection.ui_run_experiment(experiment)
    preview_results = visualisation.ui_preview_results(experiment)
    run_experiment = _bind_widgets(run_experiment, preview_results)
    return run_experiment
