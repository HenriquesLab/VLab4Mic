from supramolsim import experiments, sweep_generator
import numpy as np
import pytest


def test_run_analysis(sweep_gen):
    sweep_gen.run_analysis(save=True, plots=True)
