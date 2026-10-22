import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.psfs import elliptical_gaussian_3sigmas
from vlab4mic.utils.transform.image_convolution import frame_by_volume_convolution


def _frame(render_sigma_xy=None, photons=100, render_record=None):
    coordinates = np.array([[20.0, 20.0, 20.0]])
    photons = np.array([photons])
    ranges = [(0, 40), (0, 40), (0, 40)]
    kernel = elliptical_gaussian_3sigmas((21, 21, 21), (3, 3, 3))
    return frame_by_volume_convolution(
        coordinates, photons, ranges, 1, kernel, 20,
        projection_depth=10, render_sigma_xy=render_sigma_xy,
        render_record=render_record,
    )


@pytest.mark.parametrize(
    "sigma, photons, expected",
    [(2, 100, [True]), (0, 100, [False]), (2, 0, [False]), (None, 100, [])],
)
def test_render_record_tells_whether_the_gaussian_was_applied(sigma, photons, expected):
    render_record = []
    _frame(render_sigma_xy=sigma, photons=photons, render_record=render_record)
    assert render_record == expected


def test_histogram_rendering_keeps_each_localisation_in_one_voxel():
    frame = _frame(render_sigma_xy=0)
    assert frame.max() == pytest.approx(100)
    assert np.count_nonzero(frame) == 1


def test_rendering_kernel_blurs_laterally_and_keeps_photons():
    frame = _frame(render_sigma_xy=2)
    assert frame.sum() == pytest.approx(100, rel=1e-6)
    assert frame.max() < 100


def test_psf_convolution_is_used_without_rendering():
    frame = _frame(render_sigma_xy=None)
    # the PSF spreads the emitter laterally and axially
    assert frame.max() < 100


@pytest.mark.network
def test_smlm_localisations_are_not_convolved_with_psf():
    # regression: localisations were also convolved with the 8 nm SMLM PSF,
    # so the localisation precision did not enter once as described
    _, _, experiment = experiments.image_vsample(
        multimodal=["SMLM"], run_simulation=False, random_seed=2
    )
    imager = experiment.imager
    imager.modalities["SMLM"]["emitters"]["lateral_precision"] = 1e-6
    imager.modalities["SMLM"]["emitters"]["axial_precision"] = 1e-6
    # ROI coordinates are in micrometres; keep the emitter off voxel edges
    centre = np.array(imager.get_roi_params("ranges")).mean(axis=1) + [7e-4, 7e-4, 0]
    imager.emitters_by_fluorophore = {"AF647": [centre.tolist()]}
    np.random.seed(0)
    _, _, noiseless, _ = imager.generate_imaging(modality="SMLM", exp_time=0.001)
    image = np.asarray(noiseless["ch0"])[0]
    # every localisation of the single emitter lands in one pixel
    assert np.count_nonzero(image) == 1


def _smlm_single_emitter():
    _, _, experiment = experiments.image_vsample(
        multimodal=["SMLM"], run_simulation=False, random_seed=2
    )
    imager = experiment.imager
    centre = np.array(imager.get_roi_params("ranges")).mean(axis=1) + [7e-4, 7e-4, 0]
    imager.emitters_by_fluorophore = {"AF647": [centre.tolist()]}
    return imager


@pytest.mark.network
def test_photon_limited_render_is_filtered_and_not_noised():
    # regression: at low photon counts the rendered localisations are spread
    # below one photon per pixel, and shot noise turned them into isolated
    # single-count pixels that looked unfiltered
    imager = _smlm_single_emitter()
    imager.modalities["SMLM"]["emitters"]["precision_model"] = "photon_limited"
    np.random.seed(0)
    noisy, _, noiseless, _ = imager.generate_imaging(
        modality="SMLM", exp_time=0.001, noise=True
    )
    np.testing.assert_allclose(
        np.asarray(noisy["ch0"]), np.clip(np.asarray(noiseless["ch0"]), 0, None)
    )
    record = imager.get_positions()["SMLM"]["ch0"]["AF647"]
    assert record["parameters"]["rendering_filter_applied"] is True
    assert record["parameters"]["rendering_filtered_frames"] == 1
    assert record["parameters"]["detector_noise_applied"] is False
    assert any("detector noise not applied" in n for n in record["notes"])


@pytest.mark.network
def test_zero_kernel_is_recorded_as_unfiltered():
    imager = _smlm_single_emitter()
    imager.modalities["SMLM"]["emitters"]["rendering_kernel_nm"] = 0
    np.random.seed(0)
    imager.generate_imaging(modality="SMLM", exp_time=0.001)
    record = imager.get_positions()["SMLM"]["ch0"]["AF647"]
    assert record["parameters"]["rendering_filter_applied"] is False
    assert record["parameters"]["detector_noise_applied"] is False
    assert any("histogram" in n for n in record["notes"])


def _imager_with_kernel(kernel):
    imager, _ = experiments.build_virtual_microscope(multimodal=["SMLM"])
    imager.modalities["SMLM"]["emitters"]["rendering_kernel_nm"] = kernel
    return imager


def test_default_kernel_is_the_precision():
    imager, _ = experiments.build_virtual_microscope(multimodal=["SMLM"])
    assert imager.get_rendering_kernel("SMLM", lateral_precision=6.5) == pytest.approx(6.5)


@pytest.mark.parametrize("kernel, expected", [(0, 0.0), (2, 2.0), ("precision", 4.0)])
def test_rendering_kernel_values(kernel, expected):
    imager = _imager_with_kernel(kernel)
    assert imager.get_rendering_kernel("SMLM", lateral_precision=4.0) == expected


def test_invalid_rendering_kernel_raises():
    with pytest.raises(ValueError):
        _imager_with_kernel("blurry").get_rendering_kernel("SMLM", 4.0)
