from supramolsim import experiments, sweep_generator
import numpy as np
import pytest


def test_run_analysis(sweep_gen):
    sweep_gen._set_param_range("probe", "labelling_efficiency", "numeric", first=0.5, last=1, option=2)
    sweep_gen._set_param_range("virtual_sample", "random_orientations", "logical", option="True")
    sweep_gen._set_param_range("acquisition", "exp_time","numeric", first= 0.001, last = 0.1, option= 2)
    sweep_gen.generate_acquisitions()
    sweep_gen.generate_reference_image()
    sweep_gen.run_analysis()
    sweep_gen.gen_analysis_dataframe()
