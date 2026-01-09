from vlab4mic import experiments, sweep_generator
import numpy as np
import pytest
import copy


def test_run_parameter_sweep():
    sweep_gen_test = sweep_generator.run_parameter_sweep(
        structures=["7R5K",],
        probe_templates=["NPC_Nup96_Cterminal_direct",],
        sweep_repetitions=3,
        # parameters for sweep
        labelling_efficiency=(0,1,2),
        defect=(0,1,2),
        defect_small_cluster=[300,],
        defect_large_cluster=[600,],
        #exp_time=[0.001, 0.01,],
        # output and analysis
        output_name="vlab_script",
        return_generator=True,
        save_sweep_images=False,
        save_analysis_results=True,
        run_analysis=True
        )
    
    assert sweep_gen_test.analysis["dataframes"] is not None
    param_groups = list(sweep_gen_test.params_by_group.keys())
    total_combinations = 1
    for group_name in param_groups:
        if len(sweep_gen_test.params_by_group[group_name]) > 0:
            for param_name in sweep_gen_test.params_by_group[group_name]:
                total_combinations *= len(
                    sweep_gen_test.params_by_group[group_name][param_name]
                )
    vsamples_unique_ids = len(sweep_gen_test.virtual_samples_parameters.keys())
    assert total_combinations == vsamples_unique_ids