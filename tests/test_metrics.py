import numpy as np
import pytest

from vlab4mic.analysis.metrics import pearson_correlation, structural_similarity


def _images():
    rng = np.random.default_rng(0)
    ref = np.zeros((32, 32))
    ref[8:24, 15:17] = 10
    noisy = ref + rng.normal(0, 1, ref.shape)
    shuffled = rng.permutation(ref.flatten()).reshape(ref.shape)
    return ref, noisy, shuffled


def _ssim(a, b, mask=None):
    return structural_similarity(
        reference_image=a, reference_image_pixelsize_nm=10, reference_image_mask=mask,
        simulated_image=b, simulated_image_pixelsize_nm=10, simulated_image_mask=mask,
    )


def test_ssim_of_identical_images_is_one():
    ref, _, _ = _images()
    assert _ssim(ref, ref) == pytest.approx(1)


def test_ssim_is_spatial():
    # regression: SSIM ran on a flattened 1-D vector of masked pixels; the
    # 2-D SSIM must rank a noisy copy above a pixel-shuffled one
    ref, noisy, shuffled = _images()
    assert _ssim(ref, noisy) > _ssim(ref, shuffled)


def test_metrics_without_masks():
    # regression: no mask raised "win_size exceeds image extent"
    ref, noisy, _ = _images()
    assert np.isfinite(_ssim(ref, noisy))
    r = pearson_correlation(
        reference_image=ref, reference_image_pixelsize_nm=10,
        simulated_image=noisy, simulated_image_pixelsize_nm=10,
    )
    assert r == pytest.approx(np.corrcoef(ref.ravel(), noisy.ravel())[0, 1])


def test_ssim_with_mask():
    ref, noisy, _ = _images()
    mask = ref > 0
    assert 0 < _ssim(ref, noisy, mask) <= 1
