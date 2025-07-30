from supramolsim import workflows, experiments
from supramolsim.utils import data_format
import pytest
import copy


def test_simple_imaging_system():
    imager, _ = experiments.build_virtual_microscope()
    assert imager.get_absoulte_reference_point().shape == (1, 3)


def test_get_raw_volume(experiment_7r5k_base):
    experiment_7r5k_base
    images_volumes, beads = experiment_7r5k_base.imager.generate_imaging(
        modality="SMLM", convolution_type="raw_volume", exp_time=0.001
    )
    assert images_volumes[0].shape == (200, 200, 150)


def test_multi_imaging_system():
    imager, _ = experiments.build_virtual_microscope(multimodal=["STED", "Confocal"])
    assert imager.get_absoulte_reference_point().shape == (1, 3)


def test_image_from_field(configuration_directory, gt_structural_model_field):
    configuration_path = configuration_directory
    selected_mods = [
        "STED",
    ]
    imgs, experiment_test = experiments.image_vsample(vsample=gt_structural_model_field, run_simulation=False, modalities=selected_mods, configuration_path=configuration_path)
    assert experiment_test.imager.get_absoulte_reference_point().shape == (1, 3)


def test_imager_optional_methods(experiment_7r5k_base):
    imager = copy.deepcopy(experiment_7r5k_base.imager)
    imager.set_experiment_name("test")
    imager.set_focus(0.23)
    imager.set_roi_position(x=0.1, y=0.2)
    imager.set_roi_sizes(x=2, y=2)
    imager.recenter_roi()
    imager.show_field()
