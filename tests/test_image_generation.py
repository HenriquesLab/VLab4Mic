from vlab4mic import workflows, experiments
from vlab4mic.utils import data_format
import pytest
import copy


def test_simple_imaging_system():
    imager, _ = experiments.build_virtual_microscope()
    assert imager.get_absolute_reference_point().shape == (1, 3)


def test_get_raw_volume():
    images, noisless, testexperiment = experiments.image_vsample(
        multimodal=["SMLM",],
        run_simulation=False)
    images_volumes, beads, img_noiseless, b_noiseless = testexperiment.imager.generate_imaging(
        modality="SMLM", convolution_type="raw_volume", exp_time=0.001
    )
    assert images_volumes["raw_volume"][0].shape == (500, 500, 150)


def test_multi_imaging_system():
    imager, _ = experiments.build_virtual_microscope(
        multimodal=["STED", "Confocal"]
    )
    assert imager.get_absolute_reference_point().shape == (1, 3)


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


def test_imager_optional_methods():
    images, noisless, testexperiment = experiments.image_vsample(
        multimodal=["SMLM",],
        run_simulation=False)
    testexperiment.imager.set_experiment_name("test")
    testexperiment.imager.set_focus(0.23)
    testexperiment.imager.set_roi_position(x=0.1, y=0.2)
    testexperiment.imager.set_roi_sizes(x=2, y=2)
    testexperiment.imager.recenter_roi()
