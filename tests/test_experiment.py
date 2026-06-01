from vlab4mic import experiments
import pytest
import numpy as np
import math


def test_gen_virtual_sample():
    vsample, testexperiment = experiments.generate_virtual_sample()
    assert len(vsample["reference_point"]) == 3


def test_image_output_shape():
    images,images_noiseless, experiment = experiments.image_vsample()
    for modality in images.keys():
        assert type(images[modality]) is dict
        assert len(images[modality]["ch0"].shape) == 3


def test_build_virtual_microscope():
    vmicroscope, experiment = experiments.build_virtual_microscope()
    assert vmicroscope.get_absolute_reference_point().shape == (1, 3)
    assert len(experiment.imager.modalities) > 0


def test_build_vmicroscope_multimodal():
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images, experiment = experiments.build_virtual_microscope(
        multimodal=modalities,
    )
    for modality in modalities:
        assert modality in experiment.imager.modalities.keys()


# structure_list = ["2RCJ", "7R5K", "3J3Y", "1XI5", "1HZH"]
structure_list = [
    "2RCJ",
    "1XI5",
]


@pytest.mark.parametrize("structure", structure_list)
def test_image_sample_structures(structure):
    images,images_noiseless, experiment = experiments.image_vsample(
        structure=structure,
        labelling_efficiency=0.01,
    )
    modalityname = list(images.keys())[0]
    assert len(images[modalityname]["ch0"].shape) == 3


def test_multimodal_imaging():
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images,images_noiseless, experiment = experiments.image_vsample(multimodal=modalities)
    modalityname = list(images.keys())[0]
    assert len(images[modalityname]["ch0"].shape) == 3
    assert len(list(images.keys())) == len(modalities)


def test_download_structure():
    structure9I0K = "9I0K"
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images,images_noiseless, experiment = experiments.image_vsample(
        structure=structure9I0K, multimodal=modalities
    )
    modalityname = list(images.keys())[0]
    assert len(images[modalityname]["ch0"].shape) == 3
    assert len(list(images.keys())) == len(modalities)


def test_vsample_function():
    axial_offset = [0,100]
    yz_orientations = [20,40,60,]
    sample, experiment = experiments.generate_virtual_sample(
        random_orientations=True,
        random_placing=True,
        random_rotations=True,
        rotation_angles=[0,20,30],
        yz_orientations=yz_orientations,
        axial_offset=axial_offset,
        probe_DoL=5,
        number_of_particles=10,
        clear_experiment=True
    )
    z_pos = experiment.coordinate_field.molecules[0].params["ref_point"][2]
    axis = experiment.coordinate_field.molecules[0].axis["direction"]
    original_axis = experiment.coordinate_field.molecules[0].axis_reset["direction"]
    assert z_pos in axial_offset
    original_norm = original_axis / np.linalg.norm(original_axis)
    new_norm = axis / np.linalg.norm(axis)
    dot_product = np.dot(original_norm, new_norm)
    radians=np.arccos(np.clip(dot_product, -1.0, 1.0))
    degree = math.degrees(radians)
    found = False
    for i in yz_orientations:
        if (np.absolute(i - degree)) < 0.00001:
            found = True
    assert found
    experiment.coordinate_field.molecules[0].reset_axis_orientation()
    axis2 = experiment.coordinate_field.molecules[0].axis["direction"]
    axis_diff = axis2 - original_axis
    assert np.linalg.norm(axis_diff) < 0.00001


def test_generate_virtual_sample_probe_secondary_epitope_applied():
    """Regression test: probe_seconday_epitope passed to generate_virtual_sample
    must reach add_probe and be stored as epitope_target_info (issue #3)."""
    secondary_epitope = {"target": {"type": "Sequence", "value": "TESTEPITOPE"}}
    _vsample, experiment = experiments.generate_virtual_sample(
        probe_seconday_epitope=secondary_epitope,
    )
    probe_name = list(experiment.probe_parameters.keys())[0]
    configured_probe = experiment.probe_parameters[probe_name]
    assert configured_probe.get("epitope_target_info") == secondary_epitope


def test_add_probe_secondary_epitope_correct_spelling():
    """Regression test: add_probe must accept the correctly-spelled
    probe_secondary_epitope parameter (issue #4)."""
    from vlab4mic.experiments import Experiment

    experiment = Experiment()
    secondary_epitope = {"target": {"type": "Sequence", "value": "TESTEPITOPE2"}}
    experiment.add_probe(
        probe_template="NHS_ester",
        probe_secondary_epitope=secondary_epitope,
    )
    configured_probe = experiment.probe_parameters["NHS_ester"]
    assert configured_probe.get("epitope_target_info") == secondary_epitope