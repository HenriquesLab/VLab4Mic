"""
experiment_widgets
------------------

This module provides high-level widget composition functions for interactive experiment setup and visualisation
in Jupyter notebooks. It combines parameter selection and visualisation widgets for each step of the workflow,
using EZInput-based GUIs and ipywidgets.

Functions
---------

- _bind_widgets(parameters=None, visualisation=None, add_border=False):
    Binds parameter and visualisation widgets into a single GUI layout.

- select_structure_widget(experiment):
    Returns a widget for selecting and visualising the experiment structure.

- select_probe_widget(experiment):
    Returns a widget for selecting probes and visualising the labelled structure.

- select_sample_parameters_widget(experiment):
    Returns a widget for selecting sample parameters and visualising the virtual sample.

- select_modalities_widget(experiment):
    Returns a widget for selecting imaging modalities.

- select_acquisition_parameters_widget(experiment):
    Returns a widget for selecting acquisition parameters and previewing acquisition settings.

- run_experiment_widget(experiment):
    Returns a widget for running the experiment and previewing results.

Typical usage
-------------

    from supramolsim.jupyter_widgets import experiment_widgets
    my_experiment = experiments.ExperimentParametrisation()
    structure_widget = experiment_widgets.select_structure_widget(my_experiment)
    structure_widget.show()

Each function returns an EZInput-based widget that can be displayed in a Jupyter notebook.
"""
from . import _experiment_parameters
from ..experiments import ExperimentParametrisation
from . import _experiment_visualisation
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


def _bind_widgets(parameters=None, visualisation=None, add_border=False, merge_into_parameters=False):
    """
    Bind the widgets to the parameters and visualisation objects.

    Parameters
    ----------
    parameters : EZInput or None
        The parameter selection widget (EZInput instance).
    visualisation : EZInput or None
        The visualisation widget (EZInput instance).
    add_border : bool
        If True, adds section headers and a border between parameter and visualisation sections.

    Returns
    -------
    EZInput
        A combined EZInput widget containing both parameter and visualisation elements.
    """
    if merge_into_parameters:
        if visualisation is not None:
            if add_border:
                parameters.elements["divisor"] = widgets.Label("")
                parameters["divisor"].layout=widgets.Layout(border='1px solid black', height='0px', width='50%')
                parameters.add_HTML("visualisation_section","Visualisation section", style = dict(font_weight='bold',  font_size='16px') )

            for tag, element in visualisation.elements.items():
                parameters.elements[tag] = element
        return parameters
    else:
        gui = EZInput("Main_widget")
        if add_border:
            gui.add_HTML("parameters_section", "Parameters section", style=dict(font_weight='bold', font_size='16px'))
    
        if parameters is not None:
            for tag, element in parameters.elements.items():
                gui.elements[tag] = element
        if visualisation is not None:
            if add_border:
                gui.elements["divisor"] = widgets.Label("")
                gui["divisor"].layout=widgets.Layout(border='1px solid black', height='0px', width='50%')
                gui.add_HTML("visualisation_section","Visualisation section", style = dict(font_weight='bold',  font_size='16px') )

            for tag, element in visualisation.elements.items():
                gui.elements[tag] = element
        return gui

def select_structure_widget(experiment):
    """
    Create a widget for selecting and visualising the experiment structure.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for structure selection and visualisation.
    """
    structure_params = _experiment_parameters.ui_select_structure(experiment)
    view_structure = _experiment_visualisation.ui_show_structure(experiment)
    select_structure = _bind_widgets(structure_params, view_structure, merge_into_parameters=True)
    return select_structure

def select_probe_widget(experiment):
    """
    Create a widget for selecting probes and visualising the labelled structure.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for probe selection and labelled structure visualisation.
    """
    probes_params = _experiment_parameters.ui_select_probe(experiment)
    view_probes = _experiment_visualisation.ui_show_labelled_structure(experiment)
    select_probe = _bind_widgets(probes_params, view_probes)
    return select_probe

def select_sample_parameters_widget(experiment):
    """
    Create a widget for selecting sample parameters and visualising the virtual sample.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for sample parameter selection and virtual sample visualisation.
    """
    sample_params = _experiment_parameters.ui_select_sample_parameters(experiment)
    view_sample = _experiment_visualisation.ui_show_virtual_sample(experiment)
    select_sample = _bind_widgets(sample_params, view_sample)
    return select_sample

def select_modalities_widget(experiment):
    """
    Create a widget for selecting imaging modalities.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for modality selection.
    """
    modalities_params = _experiment_parameters.ui_select_modality(experiment)
    select_modalities = _bind_widgets(modalities_params)
    return select_modalities


def select_acquisition_parameters_widget(experiment):
    """
    Create a widget to select acquisition parameters and preview acquisition settings.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for acquisition parameter selection and preview.
    """
    view_acq_params = _experiment_visualisation.ui_set_acq_params(experiment)
    select_acq_params = _bind_widgets(visualisation = view_acq_params)
    return select_acq_params

def run_experiment_widget(experiment):
    """
    Create a widget to run the experiment and preview results.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object.

    Returns
    -------
    EZInput
        Widget for running the experiment and previewing results.
    """
    run_experiment = _experiment_parameters.ui_run_experiment(experiment)
    preview_results = _experiment_visualisation.ui_preview_results(experiment)
    run_experiment = _bind_widgets(run_experiment, preview_results)
    return run_experiment
