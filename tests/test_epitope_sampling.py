import numpy as np
import pytest

from vlab4mic.utils.sample.arrays import binomial_epitope_sampling


def flat_grid(n_side=20, spacing=2.0):
    g = np.stack(
        np.meshgrid(np.arange(n_side) * spacing, np.arange(n_side) * spacing), -1
    ).reshape(-1, 2)
    return np.c_[g, np.zeros(len(g))]


def min_pairwise_distance(points):
    d = np.linalg.norm(points[:, None] - points[None], axis=-1)
    np.fill_diagonal(d, np.inf)
    return d.min()


@pytest.mark.parametrize("min_distance", [0.0, 5.0])
def test_zero_efficiency_binds_nothing(min_distance):
    np.random.seed(0)
    _, n, _ = binomial_epitope_sampling(flat_grid(), p=0, min_distance=min_distance)
    assert n == 0


def test_full_efficiency_without_hindrance_binds_everything():
    np.random.seed(0)
    epitopes = flat_grid()
    _, n, _ = binomial_epitope_sampling(epitopes, p=1, min_distance=0.0)
    assert n == len(epitopes)


@pytest.mark.parametrize("p", [0.1, 0.5, 0.9])
def test_bound_count_is_binomial_without_hindrance(p):
    np.random.seed(1)
    epitopes = flat_grid()
    n_sites = len(epitopes)
    counts = np.array(
        [binomial_epitope_sampling(epitopes, p=p)[1] for _ in range(500)]
    )
    expected_sd = np.sqrt(n_sites * p * (1 - p))
    # mean of 500 draws is within 4 standard errors of N*p
    assert abs(counts.mean() - n_sites * p) < 4 * expected_sd / np.sqrt(500)
    assert counts.std() == pytest.approx(expected_sd, rel=0.2)


@pytest.mark.parametrize("p", [0.1, 1.0])
def test_bound_epitopes_respect_min_distance(p):
    np.random.seed(2)
    min_distance = 5.0
    for _ in range(20):
        bound, n, _ = binomial_epitope_sampling(
            flat_grid(), p=p, min_distance=min_distance
        )
        if n > 1:
            assert min_pairwise_distance(bound) >= min_distance


def test_unbound_epitopes_do_not_block_neighbours():
    # with low efficiency, hindrance should rarely matter: the bound count
    # stays close to N*p rather than being capped at p * (steric maximum)
    np.random.seed(3)
    epitopes = flat_grid()
    counts = [
        binomial_epitope_sampling(epitopes, p=0.1, min_distance=5.0)[1]
        for _ in range(200)
    ]
    # the steric maximum on this grid is ~40, so applying hindrance first
    # would give ~4 bound epitopes; binding first gives ~20
    assert np.mean(counts) > 15


def test_normals_follow_selected_epitopes():
    np.random.seed(4)
    epitopes = flat_grid()
    normals = epitopes * 10
    bound, n, bound_normals = binomial_epitope_sampling(
        epitopes, p=0.5, normals=normals, min_distance=5.0
    )
    assert bound_normals.shape == (n, 3)
    np.testing.assert_array_equal(bound_normals, bound * 10)


def test_sampling_is_reproducible_with_seed():
    epitopes = flat_grid()
    np.random.seed(5)
    first = binomial_epitope_sampling(epitopes, p=0.5, min_distance=5.0)[0]
    np.random.seed(5)
    second = binomial_epitope_sampling(epitopes, p=0.5, min_distance=5.0)[0]
    np.testing.assert_array_equal(first, second)
