from vlab4mic import experiments, sweep_generator
import numpy as np
import pytest
import copy


def test_run_default_analysis(sweep_gen):
    test_sweep = copy.deepcopy(sweep_gen)
    test_sweep.set_parameter_values(
        "probe", "labelling_efficiency", values=(0.4, 1, 3)
    )
    test_sweep.set_analysis_parameters(metrics_list=["ssim", "pearson"])
    test_sweep.run_analysis(save=True, plots=True)


def test_parameters_combinations():
    test_sweep = sweep_generator.sweep_generator()
    test_sweep.probes = [
        "CCP_heavy_chain_Cterminal",
    ]
    test_sweep.set_parameter_values(
        "probe", "labelling_efficiency", values=(0, 1, 5)
    )
    test_sweep.set_parameter_values(
        "particle_defect", "defect", values=(0.1, 0.5, 4)
    )
    test_sweep.set_parameter_values(
        "particle_defect",
        "defect_small_cluster",
        values=[
            300,
        ],
    )
    test_sweep.set_parameter_values(
        "particle_defect",
        "defect_large_cluster",
        values=[
            600,
        ],
    )
    test_sweep.set_parameter_values(
        "virtual_sample", "random_orientations", values=[True, False]
    )
    test_sweep.set_parameter_values(
        "virtual_sample", "number_of_particles", values=[1, 2, 10]
    )

    param_groups = list(test_sweep.params_by_group.keys())
    total_combinations = 1
    for group_name in param_groups:
        if len(test_sweep.params_by_group[group_name]) > 0:
            for param_name in test_sweep.params_by_group[group_name]:
                total_combinations *= len(
                    test_sweep.params_by_group[group_name][param_name]
                )
    test_sweep.generate_virtual_samples()
    vsamples_unique_ids = len(test_sweep.virtual_samples_parameters.keys())
    assert total_combinations == vsamples_unique_ids
