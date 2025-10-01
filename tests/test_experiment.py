from vlab4mic import experiments
import pytest
import numpy as np
import math


def test_gen_virtual_sample():
    vsample, testexperiment = experiments.generate_virtual_sample()
    assert len(vsample["reference_point"]) == 3


def test_image_sample():
    images,images_noiseless, experiment = experiments.image_vsample()


def test_build_virtual_microscope():
    vmicroscope, experiment = experiments.build_virtual_microscope()


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
    assert len(images[modalityname].shape) == 3


def test_multimodal_imaging():
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images,images_noiseless, experiment = experiments.image_vsample(multimodal=modalities)
    modalityname = list(images.keys())[0]
    assert len(images[modalityname].shape) == 3
    assert len(list(images.keys())) == len(modalities)


def test_download_structure():
    structure9I0K = "9I0K"
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images,images_noiseless, experiment = experiments.image_vsample(
        structure=structure9I0K, multimodal=modalities
    )
    modalityname = list(images.keys())[0]
    assert len(images[modalityname].shape) == 3
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