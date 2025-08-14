from vlab4mic import experiments
import pytest


def test_gen_virtual_sample():
    vsample, testexperiment = experiments.generate_virtual_sample()
    assert len(vsample["reference_point"]) == 3


def test_image_sample():
    images, experiment = experiments.image_vsample()


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
    images, experiment = experiments.image_vsample(
        structure=structure,
        labelling_efficiency=0.01,
    )
    modalityname = list(images.keys())[0]
    assert len(images[modalityname].shape) == 3


def test_multimodal_imaging():
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images, experiment = experiments.image_vsample(multimodal=modalities)
    modalityname = list(images.keys())[0]
    assert len(images[modalityname].shape) == 3
    assert len(list(images.keys())) == len(modalities)


def test_download_structure():
    structure9I0K = "9I0K"
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    images, experiment = experiments.image_vsample(
        structure=structure9I0K, multimodal=modalities
    )
    modalityname = list(images.keys())[0]
    assert len(images[modalityname].shape) == 3
    assert len(list(images.keys())) == len(modalities)
