import numpy as np
import pytest

from vlab4mic import experiments


@pytest.fixture(scope="module")
def imager():
    imager, _ = experiments.build_virtual_microscope(multimodal=["Confocal", "Widefield"])
    return imager


def test_camera_exposure_is_unchanged(imager):
    assert imager.get_photon_exposure("Widefield", 0.01) == 0.01


def test_scanning_uses_dwell_time_times_psf_footprint(imager):
    imager.modalities["Confocal"]["detector"]["scanning"] = True
    try:
        # Confocal: sigma 100 nm, pixel 70 nm
        expected = 1e-5 * 2 * np.pi * 100 * 100 / 70**2
        assert imager.get_photon_exposure("Confocal", 1e-5) == pytest.approx(expected)
    finally:
        imager.modalities["Confocal"]["detector"]["scanning"] = False


def test_shipped_modalities_are_not_scanning_by_default(imager):
    assert imager.modalities["Confocal"]["detector"]["scanning"] is False
