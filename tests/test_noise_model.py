import numpy as np
import pytest

from vlab4mic.utils.transform import noise


def test_stochastic_rounding_keeps_mean_of_dim_pixels():
    # regression: flooring turned 0.7 expected photons into 0
    np.random.seed(0)
    values = np.full(200000, 0.7)
    assert noise.add_binomial_noise(values, p=1.0).mean() == pytest.approx(0.7, abs=0.01)


def test_fixed_gain_has_no_excess_noise():
    np.random.seed(1)
    electrons = np.full(100000, 100)
    out = noise.add_gamma_noise(electrons, g=1.0)
    assert out.var() == pytest.approx(0)


def test_em_gain_doubles_variance():
    np.random.seed(2)
    electrons = np.random.poisson(100, 200000)
    out = noise.add_gamma_noise(electrons, g=1.0, em_gain=True)
    assert out.var() / electrons.var() == pytest.approx(2, rel=0.05)


def test_negative_values_are_clipped_not_shifted():
    from vlab4mic.generate.imaging import Imager

    stack = np.array([[-3.0, 5.0], [2.0, 0.0]])
    clipped = Imager._crop_negative(None, stack)
    np.testing.assert_array_equal(clipped, [[0, 5], [2, 0]])


def test_saturation_is_applied():
    from vlab4mic.generate.imaging import Imager

    class Stub:
        modalities = {"m": {"detector": {"bits_pixel": 8}}}

    out = Imager._adjust_to_pixel_depth(Stub(), "m", np.array([10, 300]))
    np.testing.assert_array_equal(out, [10, 255])
