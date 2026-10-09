import numpy as np
import pytest

from vlab4mic.analysis.particle_measures import (
    angular_gaps,
    count_occupied_sectors,
    count_resolved_sites,
    fit_circle,
    has_apparent_break,
    ring_measures,
)


def ring(radius=50, n=8, per_corner=5, jitter=0.5, missing=(), seed=0):
    rng = np.random.default_rng(seed)
    pts = []
    for k in range(n):
        if k in missing:
            continue
        a = 2 * np.pi * k / n
        for _ in range(per_corner):
            pts.append([radius * np.cos(a) + 100, radius * np.sin(a) - 30, 0])
    return np.array(pts) + rng.normal(0, jitter, (len(pts), 3))


def test_circle_fit():
    centre, radius = fit_circle(ring())
    np.testing.assert_allclose(centre, [100, -30], atol=0.5)
    assert radius == pytest.approx(50, abs=0.5)


def test_ring_width():
    assert ring_measures(ring(jitter=3))["width"] == pytest.approx(3, rel=0.4)


def test_breaks_and_sectors():
    assert not has_apparent_break(ring(), max_gap_deg=60)
    assert has_apparent_break(ring(missing=(2,)), max_gap_deg=60)
    assert count_occupied_sectors(ring(), 8) == 8
    assert count_occupied_sectors(ring(missing=(1, 5)), 8) == 6
    np.testing.assert_allclose(angular_gaps(ring(jitter=0)).sum(), 360)


def test_resolved_sites():
    image = np.zeros((40, 40))
    for x, y in [(10, 10), (10, 25), (25, 18)]:
        image[x, y] = 10
    from scipy.ndimage import gaussian_filter
    sharp = gaussian_filter(image, 1)
    blurred = gaussian_filter(image, 8)
    assert count_resolved_sites(sharp, 1, 5) == 3
    assert count_resolved_sites(blurred, 1, 5) < 3


def _triangle_localisations(precision, n=200, side=6.0, seed=0):
    rng = np.random.default_rng(seed)
    corners = side / np.sqrt(3) * np.array(
        [[np.cos(a), np.sin(a)] for a in np.radians([90, 210, 330])]
    )
    return np.concatenate([c + rng.normal(0, precision, (n, 2)) for c in corners])


@pytest.mark.parametrize("precision, expected", [(1, True), (1.5, True), (6.5, False), (10, False)])
def test_sites_resolved_follows_precision(precision, expected):
    from vlab4mic.analysis.particle_measures import sites_resolved

    resolved, n_components = sites_resolved(_triangle_localisations(precision), 3)
    assert resolved is expected
