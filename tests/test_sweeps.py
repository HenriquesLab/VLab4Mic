from vlab4mic import experiments, sweep_generator
import numpy as np
import pytest
import copy


@pytest.mark.network
def test_run_parameter_sweep():
    sweep_gen_test = sweep_generator.run_parameter_sweep(
        structures=[
            "7R5K",
        ],
        probe_templates=[
            "NPC_Nup96_Cterminal_direct",
        ],
        sweep_repetitions=3,
        # parameters for sweep
        labelling_efficiency=(0, 1, 1),
        structural_integrity=[0, 0.5, 1],
        structural_integrity_small_cluster=[
            300,
        ],
        structural_integrity_large_cluster=[
            600,
        ],
        # exp_time=[0.001, 0.01,],
        # output and analysis
        output_name="vlab_script",
        return_generator=True,
        save_sweep_images=False,
        save_analysis_results=False,
        run_analysis=True,
        analysis_plots=False,
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


@pytest.mark.network
def test_custom_metric():

    def mean_value(
        reference_image=None,
        reference_image_pixelsize_nm=None,
        simulated_image=None,
        simulated_image_pixelsize_nm=None,
        image_mask=None,
        resized_reference_image=None,
        resized_simulated_image=None,
        *args,
        **kwargs,
    ):
        return np.mean(simulated_image)

    sweep_gen = sweep_generator.run_parameter_sweep(
        sweep_repetitions=3,
        # parameters for sweep
        labelling_efficiency=(
            0,
            1,
            0.5,
        ),  # values between 0 and 1 with step of 0.5
        return_generator=True,
        analysis_plots=False,
        save_sweep_images=False,  # By default, the saving directory is set to the home path of the user
        save_analysis_results=False,
        run_analysis=True,
        custom_metrics=[
            mean_value,
        ],
    )

    # the custom metric must actually be registered and computed, not ignored
    assert "mean_value" in sweep_gen.metrics
    results = sweep_gen.analysis["dataframes"]
    assert "mean_value" in results.columns
    assert results["mean_value"].notna().any()


def test_probe_secondary_epitope_sweep_parameter_reaches_add_probe():
    sweep_gen = sweep_generator.sweep_generator()
    secondary_epitope = {
        "target": {
            "type": "Sequence",
            "value": "SENTINEL_EPITOPE",
        }
    }

    sweep_gen.set_sweep_parameters(probe_secondary_epitope=[secondary_epitope])
    sweep_gen.create_parameters_iterables()

    probe_parameters = sweep_gen.probe_parameters[0]
    assert probe_parameters == {"probe_secondary_epitope": secondary_epitope}

    sweep_gen.experiment.add_probe(
        probe_template="NHS_ester", **probe_parameters
    )

    configured_probe = sweep_gen.experiment.probe_parameters["NHS_ester"]
    assert configured_probe.get("epitope_target_info") == secondary_epitope


def test_sweep_based_on_image():
    # small deterministic mask with three positive pixels
    img_mask = np.zeros((8, 8))
    img_mask[2, 2] = img_mask[5, 5] = img_mask[2, 6] = 1
    pixelsize = 100  # in nm
    image_parameters = dict(
        pixelsize=pixelsize,
        mode="mask",
        npositions=3,
    )

    sweep_gen = sweep_generator.run_parameter_sweep(
        sweep_repetitions=1,
        # parameters for sweep
        labelling_efficiency=(0.5, 1, 0.5),
        modalities=["Widefield"],
        # Image for virtual sample positioning
        image4vsample=img_mask,
        image4vsample_parameters=image_parameters,
        # reference image
        reference_image=img_mask,
        reference_image_parameters={"ref_pixelsize": pixelsize},
        # output and analysis
        output_name="vlab_example_sweep",
        return_generator=True,
        analysis_plots=False,  # Do not create plots on tests
        save_sweep_images=False,  # By default, the saving directory is set to the home path of the user
        save_analysis_results=False,
        run_analysis=True,
        capture_outputs=True,
        random_seed=1,
    )

    # one particle per positive pixel, sample sized to the image
    positions = sweep_gen.experiment.virtualsample_params["relative_positions"]
    assert np.shape(positions) == (3, 3)
    assert sweep_gen.experiment.virtualsample_params["sample_dimensions"][:2] == [
        img_mask.shape[0] * pixelsize,
        img_mask.shape[1] * pixelsize,
    ]
    # the provided image is used as reference
    np.testing.assert_array_equal(sweep_gen.reference_image, img_mask)
    # one row per labelling efficiency value, with metrics computed
    results = sweep_gen.analysis["dataframes"]
    assert sorted(results["labelling_efficiency"].unique()) == [0.5, 1.0]
    assert results[["ssim", "pearson"]].notna().all().all()
