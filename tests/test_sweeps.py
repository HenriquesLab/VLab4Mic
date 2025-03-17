#from supramolsim.analysis.scripts import (
#    parameter_sweep_reps,
#    analyse_sweep,
#)
from supramolsim.analysis import sweep, scripts
from supramolsim import experiments
import numpy as np


def test_simple_param_sweep(experiment_7r5k_base):
    assert experiment_7r5k_base.generators_status("structure") == True
    assert experiment_7r5k_base.generators_status("imager") == True
    sweep_pars = dict(
        labelling_efficiency=dict(start=0.3, end=1, nintervals=2, ideal=1),
        defects=dict(start=0, end=0.5, nintervals=2, ideal=0),
    )
    total_combinations = (
        sweep_pars["defects"]["nintervals"]
        * sweep_pars["labelling_efficiency"]["nintervals"]
    )
    sweep_out, sweep_out_pars, ref_out, ref_params = scripts.parameter_sweep_reps(
        Experiment=experiment_7r5k_base, sweep_parameters=sweep_pars, repetitions=3
    )
    assert len(sweep_out_pars) == total_combinations


def test_sweep_analysis(experiment_7r5k_base):
    assert experiment_7r5k_base.generators_status("structure") == True
    assert experiment_7r5k_base.generators_status("imager") == True
    sweep_pars = dict(
        labelling_efficiency=dict(start=0.3, end=1, nintervals=2, ideal=1),
        defects=dict(start=0, end=0.5, nintervals=2, ideal=0),
    )
    total_combinations = (
        sweep_pars["defects"]["nintervals"]
        * sweep_pars["labelling_efficiency"]["nintervals"]
    )
    replicas = 3
    sweep_out, sweep_out_pars, ref_out, ref_params = scripts.parameter_sweep_reps(
        Experiment=experiment_7r5k_base,
        sweep_parameters=sweep_pars,
        repetitions=replicas,
    )
    # get pixelsizes
    imager_scale = experiment_7r5k_base.imager.roi_params["scale"]
    scalefactor = np.ceil(imager_scale / 1e-9)
    sweep_analyse_pars = dict(
        STED=dict(
            metric="ssim",
            subregion=False,
            force_match=True,
            modality_pixelsize=experiment_7r5k_base.imager.modalities["STED"][
                "detector"
            ]["pixelsize"]
            * scalefactor,
            ref_pixelsize=ref_params["STED"]["ref_pixelsize"],
        ),
        Confocal=dict(
            metric="ssim",
            subregion=False,
            force_match=True,
            modality_pixelsize=experiment_7r5k_base.imager.modalities["Confocal"][
                "detector"
            ]["pixelsize"]
            * scalefactor,
            ref_pixelsize=ref_params["Confocal"]["ref_pixelsize"],
        ),
    )
    conditions = len(list(sweep_analyse_pars.keys()))
    data_frame, qries, references = scripts.analyse_sweep(
        sweep_out, sweep_out_pars, ref_out, sweep_analyse_pars
    )
    assert len(data_frame["Metric"]) == total_combinations * replicas * conditions


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


def test_sweep_vasmples_indirectprobes():
    test_experiment = experiments.ExperimentParametrisation()
    structures = ["1XI5",]
    indirectprobes = ["Mock_antibody", "Mock_linker", ]
    # prepare probe parameters to sweep
    probe_parameters_vectors = dict(
        target_info = [dict(type="Sequence", value="EQATETQ"),],
        labelling_efficiency = np.linspace(0.5,1,2),
        distance_to_epitope = np.linspace(10,200,2)
    )
    # generate probe combinations
    probe_parameters = sweep.create_probe_param_combinations(**probe_parameters_vectors)
    repetitions = 3
    experiment, outputs, params = sweep.sweep_vasmples(
        experiment=test_experiment,
        structures=structures,
        probes=indirectprobes,
        probe_parameters=probe_parameters,
        repetitions=repetitions
    )
    assert len(outputs.keys()) > 0
    assert len(params.keys()) > 0
