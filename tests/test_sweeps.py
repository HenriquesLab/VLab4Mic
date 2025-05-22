from supramolsim import experiments, sweep_generator
import numpy as np
import pytest
import copy


def test_run_default_analysis(sweep_gen):
    test_sweep = copy.deepcopy(sweep_gen)
    test_sweep.run_analysis(save=True, plots=True)


def test_set_parameter_values(sweep_gen):
    test_sweep = copy.deepcopy(sweep_gen)
    test_sweep.set_parameter_values("probe", "labelling_efficiency", values=(0.5,1,2))
    test_sweep.set_parameter_values("virtual_sample", "random_orientations", values=[True, False])
    test_sweep.set_analysis_parameters("metrics_list", ["ssim", "pearson"])
    test_sweep.set_plot_parameters("heatmaps", category="modality_name", param1="labelling_efficiency", param2="random_orientations")
    test_sweep.run_analysis(save=True, plots=True)