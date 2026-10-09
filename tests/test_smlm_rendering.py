import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.psfs import elliptical_gaussian_3sigmas
from vlab4mic.utils.transform.image_convolution import frame_by_volume_convolution


def _frame(render_sigma_xy=None):
    coordinates = np.array([[20.0, 20.0, 20.0]])
    photons = np.array([100])
    ranges = [(0, 40), (0, 40), (0, 40)]
    kernel = elliptical_gaussian_3sigmas((21, 21, 21), (3, 3, 3))
    return frame_by_volume_convolution(
        coordinates, photons, ranges, 1, kernel, 20,
        projection_depth=10, render_sigma_xy=render_sigma_xy,
    )


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
