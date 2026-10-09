import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.utils.transform.points_transforms import (
    photon_limited_lateral_precision,
)


def test_precision_without_background_or_pixelation():
    # Mortensen et al. 2010: sigma^2 = (16/9) s^2 / N in this limit
    precision = photon_limited_lateral_precision(
        photons=1000, psf_sigma_nm=100, pixelsize_nm=1e-6, background_photons=0
    )
    assert precision == pytest.approx(100 * 4 / 3 / np.sqrt(1000), rel=1e-6)


def test_precision_scales_with_inverse_sqrt_photons_without_background():
    p1 = photon_limited_lateral_precision(100, 100, 100, background_photons=0)
    p2 = photon_limited_lateral_precision(10000, 100, 100, background_photons=0)
    assert p1 / p2 == pytest.approx(10)


def test_precision_worsens_with_background_and_emccd():
    base = photon_limited_lateral_precision(1000, 100, 100, background_photons=0)
    with_bg = photon_limited_lateral_precision(1000, 100, 100, background_photons=20)
    emccd = photon_limited_lateral_precision(
        1000, 100, 100, background_photons=0, excess_noise_factor=2
    )
    assert with_bg > base
    assert emccd == pytest.approx(base * np.sqrt(2))


def test_precision_needs_photons():
    with pytest.raises(ValueError):
        photon_limited_lateral_precision(0, 100, 100)


@pytest.fixture(scope="module")
def smlm_imager():
    imager, _ = experiments.build_virtual_microscope(multimodal=["SMLM"])
    imager.fluorophore_params["test_fluo"] = {"photons_per_second": 100000}
    return imager


def test_fixed_precision_is_default(smlm_imager):
    emitters = smlm_imager.modalities["SMLM"]["emitters"]
    assert smlm_imager.get_localisation_precision(
        "SMLM", "test_fluo", exp_time=0.001
    ) == (emitters["lateral_precision"], emitters["axial_precision"])


def test_photon_limited_precision_follows_exposure(smlm_imager):
    emitters = smlm_imager.modalities["SMLM"]["emitters"]
    emitters["precision_model"] = "photon_limited"
    try:
        short = smlm_imager.get_localisation_precision("SMLM", "test_fluo", 0.001)
        long = smlm_imager.get_localisation_precision("SMLM", "test_fluo", 0.01)
        expected = photon_limited_lateral_precision(
            100,
            emitters["detection_psf_sigma_nm"],
            emitters["camera_pixelsize_nm"],
            emitters["background_photons"],
            emitters["excess_noise_factor"],
        )
        assert short[0] == pytest.approx(expected)
        assert short[1] == pytest.approx(expected * emitters["axial_precision_ratio"])
        assert long[0] < short[0]
    finally:
        emitters["precision_model"] = "fixed"


def test_unknown_precision_model_raises(smlm_imager):
    emitters = smlm_imager.modalities["SMLM"]["emitters"]
    emitters["precision_model"] = "unknown"
    try:
        with pytest.raises(ValueError):
            smlm_imager.get_localisation_precision("SMLM", "test_fluo", 0.001)
    finally:
        emitters["precision_model"] = "fixed"


def test_update_modality_sets_precision_model():
    _, _, experiment = experiments.image_vsample(
        multimodal=["SMLM"], run_simulation=False, random_seed=1
    )
    experiment.update_modality(
        modality_name="SMLM",
        precision_model="photon_limited",
        background_photons=5,
    )
    emitters = experiment.imager.modalities["SMLM"]["emitters"]
    assert emitters["precision_model"] == "photon_limited"
    assert emitters["background_photons"] == 5
