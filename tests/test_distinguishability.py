import numpy as np
import pytest

from vlab4mic.analysis.distinguishability import (
    distinguishability,
    image_features,
    radial_profile_features,
)


def ring_image(radius, rng, size=41, n=60, sigma=1.5, photons=50, noise=1.0):
    angles = rng.uniform(0, 2 * np.pi, n)
    centre = size / 2 + rng.normal(0, 1, 2)
    x = centre[0] + radius * np.cos(angles) + rng.normal(0, sigma, n)
    y = centre[1] + radius * np.sin(angles) + rng.normal(0, sigma, n)
    image, _, _ = np.histogram2d(x, y, bins=size, range=[[0, size], [0, size]])
    image = image * photons + rng.normal(0, noise, image.shape)
    return image


def test_features_are_rotation_invariant():
    rng = np.random.default_rng(0)
    image = ring_image(8, rng)
    a = radial_profile_features(image)
    b = radial_profile_features(np.rot90(image))
    np.testing.assert_allclose(a, b, atol=1e-9)


def test_feature_matrix_shape_and_frames():
    rng = np.random.default_rng(1)
    images = [ring_image(8, rng)[None].repeat(3, axis=0) for _ in range(4)]
    X = image_features(images)
    assert X.shape == (4, 20)


def test_different_structures_are_distinguished():
    rng = np.random.default_rng(2)
    a = [ring_image(8, rng) for _ in range(40)]
    b = [ring_image(12, rng) for _ in range(40)]
    result = distinguishability(a, b, n_bootstrap=200, random_state=0)
    assert result["auc"] > 0.95
    assert result["accuracy"] > 0.9
    low, high = result["auc_interval"]
    assert low <= result["auc"] <= high


def test_identical_structures_are_not_distinguished():
    rng = np.random.default_rng(3)
    a = [ring_image(10, rng) for _ in range(40)]
    b = [ring_image(10, rng) for _ in range(40)]
    result = distinguishability(a, b, n_bootstrap=200, random_state=0)
    low, high = result["auc_interval"]
    assert low < 0.5 < high or abs(result["auc"] - 0.5) < 0.15


def test_needs_two_images_per_class():
    rng = np.random.default_rng(4)
    with pytest.raises(ValueError):
        distinguishability([ring_image(8, rng)], [ring_image(9, rng)] * 3)


def ellipse_image(a, b, rng, size=61, n=400, angle=None, scale=1.0):
    angle = rng.uniform(0, np.pi) if angle is None else angle
    u = rng.uniform(-1, 1, (n * 3, 2))
    u = u[(u ** 2).sum(axis=1) <= 1][:n] * [a * scale, b * scale]
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    pts = u @ rot.T + size / 2
    image, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], bins=size, range=[[0, size], [0, size]])
    return image + rng.normal(0, 0.2, image.shape)


def test_aspect_ratio_is_size_and_rotation_invariant():
    from vlab4mic.analysis.distinguishability import aspect_ratio_features

    rng = np.random.default_rng(5)
    small = aspect_ratio_features(ellipse_image(12, 6, rng, angle=0.3))[0]
    large = aspect_ratio_features(ellipse_image(12, 6, rng, angle=1.2, scale=1.8))[0]
    assert small == pytest.approx(2, rel=0.15)
    assert large == pytest.approx(small, rel=0.1)


def test_aspect_ratio_separates_elongated_from_round_regardless_of_size():
    rng = np.random.default_rng(6)
    # round particles larger than the elongated ones: size alone would mislead
    round_ = [ellipse_image(14, 14, rng, scale=rng.uniform(1.0, 1.4)) for _ in range(30)]
    elongated = [ellipse_image(13, 8, rng, scale=rng.uniform(0.7, 1.0)) for _ in range(30)]
    result = distinguishability(round_, elongated, features="aspect_ratio",
                                n_bootstrap=100, random_state=0)
    assert result["auc"] > 0.95
