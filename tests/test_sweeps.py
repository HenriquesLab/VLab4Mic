from supramolsim import experiments, sweep_generator
import numpy as np
import pytest
import copy


def test_run_default_analysis(sweep_gen):
    test_sweep = copy.deepcopy(sweep_gen)
    test_sweep.run_analysis(save=False, plots=False)
    test_sweep.set_parameter_values("probe", "labelling_efficiency", values=(0.4,1,3))
    test_sweep.set_analysis_parameters(metrics_list = ["ssim", "pearson"])
    test_sweep.run_analysis(save=True, plots=True)