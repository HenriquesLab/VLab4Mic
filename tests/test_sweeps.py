from supramolsim.analysis import sweep
from supramolsim import experiments, sweep_generator
import numpy as np
import pytest


def test_sweep_generator():
    sweep_gen = sweep_generator.sweep_generator()
    sweep_gen.generate_acquisitions()
    sweep_gen.generate_reference_image()
    sweep_gen.run_analysis()
    sweep_gen.gen_analysis_dataframe()

def test_sweep_vasmples_empty():
    test_experiment = experiments.ExperimentParametrisation()
    test_experiment, outputs, params = sweep.sweep_vasmples(test_experiment)


def test_sweep_vasmples_directprobes():
    test_experiment = experiments.ExperimentParametrisation()
    structures = [
        "1XI5",
    ]
    directprobes = [
        "NHS_ester",
    ]
    repetitions = 50
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=directprobes,
        repetitions=repetitions,
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0


indirectprobes = [
    "Antibody",
    "Linker",
]
labelling_efficiency = np.linspace(0.5, 1, 3)
distance_to_epitope = np.linspace(10, 200, 3)


@pytest.mark.parametrize("indirectprobe", indirectprobes)
@pytest.mark.parametrize("efficiency", labelling_efficiency)
@pytest.mark.parametrize("distance", distance_to_epitope)
def test_sweep_vasmples_indirectprobes(indirectprobe, efficiency, distance):
    test_experiment = experiments.ExperimentParametrisation()
    structures = [
        "1XI5",
    ]
    probes = [
        indirectprobe,
    ]
    # prepare probe parameters to sweep
    probe_parameters_vectors = dict(
        target_info=[
            dict(type="Sequence", value="EQATETQ"),
        ],
        labelling_efficiency=[
            efficiency,
        ],
        distance_to_epitope=[
            distance,
        ],
    )
    # generate probe combinations
    probe_parameters = sweep.create_param_combinations(**probe_parameters_vectors)
    repetitions = 3
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=probes,
        probe_parameters=probe_parameters,
        repetitions=repetitions,
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0

def test_sweep_modalities_updatemod():
    myexperiment, vsmpl_output, vsampl_pars = sweep.sweep_vasmples()
    myexperiment, mod_outputs, mod_params, mod_pixelsizes  = sweep.sweep_modalities_updatemod(
        experiment = myexperiment,
        vsample_outputs = vsmpl_output,
        vsampl_pars = vsampl_pars,
        )


def test_parameter_generators():
    structures = ["1XI5", ]
    probes=["NHS_ester", ]
    modalities = ["STED",]
    probe_parameters = sweep.probe_parameters_sweep(labelling_efficiency=[0.5,1,1])
    vsample_parameters = sweep.virtual_sample_parameters_sweep(random_orientations=[False,])
    modality_parameters = sweep.modality_parameters_sweep(lateral_resolution_nm=[10,100,1])
    acquisition_params = sweep.acquisition_parameters_sweep(exp_time=[0.001, 0.1, 1])
    ref_probe_parameters = dict(
        target_info = dict(type="Sequence", value="EQATETQ"),
        labelling_efficiency = 1,
        distance_to_epitope = 0
    )
    ref_vsample, ref_params = sweep.generate_global_reference_sample(structure=structures[0], probe="Linker", probe_parameters=ref_probe_parameters)
    ref_image, ref_image_pars = sweep.generate_global_reference_modality(reference_vsample=ref_vsample, reference_vsample_params=ref_params)
    myexperiment, vsmpl_output, vsampl_pars = sweep.sweep_vasmples(
        structures=structures,
        probes=probes, 
        probe_parameters=probe_parameters,
        virtual_samples=vsample_parameters,
        repetitions=3
    )
    myexperiment, mod_outputs, mod_params, mod_pixelsizes  = sweep.sweep_modalities_updatemod(
        experiment = myexperiment,
        vsample_outputs = vsmpl_output,
        vsampl_pars = vsampl_pars,
        modalities=modalities,
        modality_params=modality_parameters,
        modality_acq_prams=acquisition_params
    )
    measurements, inputs = sweep.analyse_sweep_single_reference(mod_outputs, mod_params, ref_image[0], ref_image_pars)