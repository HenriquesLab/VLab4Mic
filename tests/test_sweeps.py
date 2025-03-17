from supramolsim.analysis import sweep
from supramolsim import experiments
import numpy as np
import pytest

def test_sweep_vasmples_empty():
    test_experiment = experiments.ExperimentParametrisation()
    test_experiment, outputs, params = sweep.sweep_vasmples(test_experiment)

def test_sweep_vasmples_directprobes():
    test_experiment = experiments.ExperimentParametrisation()
    structures = ["1XI5",]
    directprobes = ["NHS_ester",]
    repetitions = 3
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=directprobes,
        repetitions=repetitions
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0

indirectprobes = ["Mock_antibody", "Mock_linker", ]
labelling_efficiency = np.linspace(0.5, 1, 3)
distance_to_epitope = np.linspace(10,200,3)
@pytest.mark.parametrize("indirectprobe", indirectprobes)
@pytest.mark.parametrize("efficiency", labelling_efficiency)
@pytest.mark.parametrize("distance", distance_to_epitope)
def test_sweep_vasmples_indirectprobes(indirectprobe, efficiency, distance):
    test_experiment = experiments.ExperimentParametrisation()
    structures = ["1XI5",]
    probes = [indirectprobe, ]
    # prepare probe parameters to sweep
    probe_parameters_vectors = dict(
        target_info = [dict(type="Sequence", value="EQATETQ"),],
        labelling_efficiency = [efficiency, ],
        distance_to_epitope = [distance, ]
    )
    # generate probe combinations
    probe_parameters = sweep.create_probe_param_combinations(**probe_parameters_vectors)
    repetitions = 3
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=probes,
        probe_parameters=probe_parameters,
        repetitions=repetitions
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0
