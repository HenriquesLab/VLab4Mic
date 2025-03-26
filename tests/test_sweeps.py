from supramolsim.analysis import sweep
from supramolsim import experiments
import numpy as np
import pytest


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
    "Mock_antibody",
    "Mock_linker",
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


def test_sweep_vsample_modality_analysis():
    structures = [
        "1XI5",
    ]
    test_experiment = experiments.ExperimentParametrisation()
    # prepare probe parameters to sweep
    probe_parameters_vectors = dict(
        target_info=[
            dict(type="Sequence", value="EQATETQ"),
        ],
        labelling_efficiency=np.linspace(0.5, 1, 2),
        distance_to_epitope=np.linspace(10, 200, 2),
    )
    indirectprobes = [
        "Mock_antibody",
        "Mock_linker",
    ]
    # generate probe combinations
    probe_parameters = sweep.create_param_combinations(**probe_parameters_vectors)
    repetitions = 3
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=indirectprobes,
        probe_parameters=probe_parameters,
        repetitions=repetitions,
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0
    modalities = dict(
        STED=dict(nframes=2),  # change nframes
        Confocal=dict(nframes=2),
    )
    test_experiment, mod_outputs, mod_params, mod_pixelsizes = sweep.sweep_modalities(
        test_experiment, outputs, params, modalities=modalities
    )

    ref_vsample, ref_params = sweep.generate_global_reference_sample(
        structure="1XI5", probe="NHS_ester"
    )
    ref_image, ref_image_pars = sweep.generate_global_reference_modality(
        reference_vsample=ref_vsample, reference_vsample_params=ref_params
    )
    #
    sweep_analyse_parameters = dict()
    ref_pixelsize = ref_image_pars["ref_pixelsize"]

    for mod_name, pixelsize in mod_pixelsizes.items():
        sweep_analyse_parameters[mod_name] = {}
        sweep_analyse_parameters[mod_name]["metric"] = "ssim"
        sweep_analyse_parameters[mod_name]["modality_pixelsize"] = pixelsize
        sweep_analyse_parameters[mod_name]["ref_pixelsize"] = ref_pixelsize
        sweep_analyse_parameters[mod_name]["force_match"] = True

    measurements, inputs = sweep.analyse_image_sweep(
        mod_outputs, mod_params, ref_image, sweep_analyse_parameters
    )
    dframe, combined = sweep.measurements_dataframe(measurements, probe_parameters)

    df_categories, titles = sweep.pivot_dataframes_byCategory(
        combined, "modality", "labelling_efficiency", "distance_to_epitope"
    )
