from vlab4mic import workflows, experiments
from vlab4mic.utils import data_format
import pytest
import copy


def test_simple_imaging_system():
    imager, _ = experiments.build_virtual_microscope()
    assert imager.get_absolute_reference_point().shape == (1, 3)


@pytest.mark.network
def test_get_raw_volume():
    images, noisless, testexperiment = experiments.image_vsample(
        multimodal=[
            "SMLM",
        ],
        run_simulation=False,
    )
    images_volumes, beads, img_noiseless, b_noiseless = (
        testexperiment.imager.generate_imaging(
            modality="SMLM", convolution_type="raw_volume", exp_time=0.001
        )
    )
    assert images_volumes["raw_volume"][0].shape == (500, 500, 150)


def test_multi_imaging_system():
    imager, _ = experiments.build_virtual_microscope(
        multimodal=["STED", "Confocal"]
    )
    assert imager.get_absolute_reference_point().shape == (1, 3)


@pytest.mark.network
def test_image_from_field(configuration_directory, gt_structural_model_field):
    configuration_path = configuration_directory
    selected_mods = [
        "STED",
    ]
    imgs, imgs_noiseless, experiment_test = experiments.image_vsample(
        vsample=gt_structural_model_field,
        run_simulation=True,
        save=True,
        multimodal=selected_mods,
    )
    assert experiment_test.imager.get_absolute_reference_point().shape == (
        1,
        3,
    )


@pytest.mark.network
def test_imager_optional_methods():
    images, noisless, testexperiment = experiments.image_vsample(
        multimodal=[
            "SMLM",
        ],
        run_simulation=False,
    )
    testexperiment.imager.set_experiment_name("test")
    testexperiment.imager.set_focus(0.23)
    testexperiment.imager.set_roi_position(x=0.1, y=0.2)
    testexperiment.imager.set_roi_sizes(x=2, y=2)
    testexperiment.imager.recenter_roi()
    # verify the setters actually took effect on the imager's ROI params
    assert testexperiment.imager.roi_params["dimension_sizes"] == [2, 2, 1]
    # recenter_roi resets the reference point back to the field centre
    assert testexperiment.imager.roi_params["reference_point"] == [0.5, 0.5]
    # set_focus must persist after recenter_roi (recenter_roi is lateral-only)
    assert testexperiment.imager.roi_params["focus_plane"] == 0.23


@pytest.mark.network
def test_smlm_with_locs():
    images, noisless, testexperiment = experiments.image_vsample(
        multimodal=[
            "SMLM",
        ],
        run_simulation=True,
    )
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"][
            "lateral_precision"
        ]
        is not None
    )
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"]["axial_precision"]
        is not None
    )
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"]["nlocalisations"]
        is not None
    )
    testexperiment.update_modality(
        modality_name="SMLM", simulate_localistations=False
    )
    testexperiment.build(modules=["imager"])
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"][
            "lateral_precision"
        ]
        is None
    )
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"]["axial_precision"]
        is None
    )
    assert (
        testexperiment.imager.modalities["SMLM"]["emitters"]["nlocalisations"]
        is None
    )
