import warnings

import numpy as np
import pytest

from vlab4mic.utils.transform.structural_integrity import xmersubset_byclustering

PARAMS = dict(eps1=2, minsamples1=1, eps2=15, minsamples2=1)


def line_of_subunits(n=40, spacing=10):
    # n subunits along x, two epitopes each, 1 unit apart
    points = []
    for k in range(n):
        points += [[spacing * k, 0, 0], [spacing * k + 1, 0, 0]]
    return np.array(points, dtype=float)


def kept_subunits(points, removed_fraction, neigh_distance=25, seed=0):
    np.random.seed(seed)
    keep = np.array(
        xmersubset_byclustering(
            points, PARAMS, deg_dissasembly=removed_fraction,
            xmer_neigh_distance=neigh_distance, return_ids=True,
        )
    )
    return np.flatnonzero(keep[::2])


@pytest.mark.parametrize("seed", range(10))
def test_kept_region_is_one_connected_piece(seed):
    # regression: when the remainder split into several clusters, the
    # removed patch was kept instead of the remainder
    kept = kept_subunits(line_of_subunits(), 0.5, seed=seed)
    assert len(kept) > 0
    assert np.all(np.diff(kept) == 1), kept


@pytest.mark.parametrize("seed", range(5))
def test_kept_fraction_matches_request(seed):
    kept = kept_subunits(line_of_subunits(), 0.5, seed=seed)
    assert abs(len(kept) / 40 - 0.5) <= 0.075


def test_no_crash_when_subunits_have_no_neighbour():
    # regression: neighbour distance below the subunit spacing raised
    # ValueError ('a' cannot be empty)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kept = kept_subunits(line_of_subunits(), 0.5, neigh_distance=5)
    assert len(kept) > 0


def test_full_and_no_removal():
    points = line_of_subunits()
    assert len(kept_subunits(points, 0)) == 40
    assert len(kept_subunits(points, 1)) == 0


def test_verification_removes_floating_fragment_and_keeps_proposal_removed():
    # regression: with a split remainder, the function returned the IDs to
    # keep while the caller removes the returned IDs, so the proposed
    # removal was kept and the remainder dropped
    from vlab4mic.utils.transform.structural_integrity import (
        singlecluster_verification,
    )
    centres = np.array([[10.0 * k, 0, 0] for k in range(40)])
    ids = np.arange(40)
    proposal = list(range(10, 20))
    removed = set(singlecluster_verification(centres, ids, proposal, 15, 1).tolist())
    assert set(proposal) <= removed
    # the smaller piece (0-9) floats and is removed; 20-39 is kept
    assert removed == set(range(0, 20))


def test_walk_that_cannot_reach_target_stops_early(capsys):
    # regression: when fewer subunits were reachable than had to be
    # removed, the walk looped until max_iter (100000 iterations)
    import time

    import scipy.spatial

    from vlab4mic.utils.transform.structural_integrity import (
        reachable_subunits,
        xmer_ids_remove,
    )

    # two groups of 5 subunits, 100 apart: from subunit 0 only 5 are reachable
    centres = np.array(
        [[10.0 * k, 0, 0] for k in range(5)] + [[100 + 10.0 * k, 0, 0] for k in range(5)]
    )
    tree = scipy.spatial.cKDTree(centres)
    assert reachable_subunits(0, tree, 3, 15) == {0, 1, 2, 3, 4}
    np.random.seed(0)
    start = time.perf_counter()
    removed = xmer_ids_remove(0, tree, 0.8, 3, 15)  # asks for 8
    assert time.perf_counter() - start < 1
    assert removed == [0, 1, 2, 3, 4]
    assert "overiteration" not in capsys.readouterr().out


def test_walk_reaches_target_when_enough_subunits_are_reachable():
    import scipy.spatial

    from vlab4mic.utils.transform.structural_integrity import xmer_ids_remove

    centres = np.array([[10.0 * k, 0, 0] for k in range(40)])
    tree = scipy.spatial.cKDTree(centres)
    np.random.seed(1)
    removed = xmer_ids_remove(0, tree, 0.5, 3, 15)
    assert len(removed) >= 20
