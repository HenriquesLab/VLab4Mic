import numpy as np
import pytest
import skimage

from vlab4mic import experiments
from vlab4mic.utils.transform.image_convolution import lateral_binning_stack


@pytest.mark.parametrize("pixel, n_out", [(25, 40), (12.5, 80), (45, 22), (100, 10)])
def test_binning_any_pixel_size_conserves_intensity(pixel, n_out):
    # regression for issues #96 and #97: pixel sizes that are not a whole
    # multiple of the voxel, or not integers, crashed
    rng = np.random.default_rng(0)
    stack = rng.random((2, 100, 100))
    out = lateral_binning_stack(stack, (10, 10), (pixel, pixel), (n_out, n_out))
    assert out.shape == (2, n_out, n_out)
    # the field covers n_out * pixel nm, a whole number of 10 nm voxels here
    covered = int(round(n_out * pixel / 10))
    np.testing.assert_allclose(
        out.sum(axis=(1, 2)), stack[:, :covered, :covered].sum(axis=(1, 2))
    )


def test_binning_matches_block_sum_for_integer_ratios():
    rng = np.random.default_rng(1)
    stack = rng.random((1, 100, 100))
    out = lateral_binning_stack(stack, (10, 10), (100, 100), (10, 10))
    expected = skimage.measure.block_reduce(stack[0], (10, 10), np.sum)
    np.testing.assert_allclose(out[0], expected)


def test_binning_non_square_field():
    stack = np.ones((1, 100, 50))
    out = lateral_binning_stack(stack, (10, 10), (100, 100), (10, 5))
    assert out.shape == (1, 10, 5)
    np.testing.assert_allclose(out, 100)


@pytest.fixture(scope="module")
def confocal_experiment():
    _, _, experiment = experiments.image_vsample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        multimodal=["Confocal"], run_simulation=False, random_seed=1,
    )
    return experiment


@pytest.mark.network
@pytest.mark.parametrize("pixel", [45, 45.5, 40])
def test_confocal_runs_with_any_pixel_size(confocal_experiment, pixel):
    confocal_experiment.update_modality("Confocal", pixel, 170, 340)
    images, _ = confocal_experiment.run_simulation()
    image = np.asarray(images["Confocal"]["ch0"])
    assert image.shape[1] == int(1000 / pixel)


@pytest.mark.network
def test_voxel_change_keeps_physical_psf(confocal_experiment):
    params = confocal_experiment.imaging_modalities["Confocal"]["psf_params"]
    sigma_nm = [sd * v for sd, v in zip(params["std_devs"], params["voxelsize"])]
    depth_nm = params["depth"] * params["voxelsize"][2]
    confocal_experiment.update_modality("Confocal", psf_voxel_nm=5)
    params = confocal_experiment.imaging_modalities["Confocal"]["psf_params"]
    np.testing.assert_allclose(
        [sd * v for sd, v in zip(params["std_devs"], params["voxelsize"])], sigma_nm
    )
    assert params["depth"] * params["voxelsize"][2] == pytest.approx(depth_nm, abs=5)
