from supramolsim import experiments, sweep_generator
import numpy as np
import pytest


def test_run_analysis(sweep_gen):
    sweep_gen.generate_acquisitions()
    sweep_gen.generate_reference_image()
    sweep_gen.run_analysis()
    sweep_gen.gen_analysis_dataframe()
