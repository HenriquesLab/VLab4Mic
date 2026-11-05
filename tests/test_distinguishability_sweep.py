import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from vlab4mic import experiments
from vlab4mic.analysis import distinguishability as dist
from vlab4mic.analysis.distinguishability import (
    distinguishability_sweep,
    plot_distinguishability_sweep,
    scores_vs_parameter,
    swept_parameters,
)


def ring_image(radius, rng, size=41, n=60, sigma=1.5):
    angles = rng.uniform(0, 2 * np.pi, n)
    x = size / 2 + radius * np.cos(angles) + rng.normal(0, sigma, n)
    y = size / 2 + radius * np.sin(angles) + rng.normal(0, sigma, n)
    image, _, _ = np.histogram2d(x, y, bins=size, range=[[0, size], [0, size]])
    return image * 50 + rng.normal(0, 1, image.shape)


@pytest.fixture
def fake_replicates(monkeypatch):
    """run_replicates stand-in: ring radius 8 for "a", 8 + gap for "b"."""
    calls = []

    def run_replicates(n_replicates=10, modality="All", export_directory=None, **kwargs):
        calls.append(kwargs)
        rng = np.random.default_rng(kwargs["random_seed"])
        gap = kwargs.get("gap", 0)
        radius = 8 + (gap if kwargs["structure"] == "b" else 0)
        results = []
        for k in range(n_replicates):
            images = {m: {"ch0": ring_image(radius, rng)} for m in kwargs["multimodal"]}
            results.append(dict(replicate=k, images=images, images_noiseless=images))
        return results, None

    monkeypatch.setattr(experiments, "run_replicates", run_replicates)
    return calls


def test_sweep_table_and_score_trend(fake_replicates):
    table = distinguishability_sweep(
        "a", "b", sweep={"gap": [0, 4]}, modalities=["STED", "SMLM"],
        n_replicates=30, random_seed=1, n_bootstrap=100, random_state=0, progress=False,
    )
    assert len(table) == 4
    assert swept_parameters(table) == ["gap"]
    assert set(table.modality) == {"STED", "SMLM"}
    same = table[table.gap == 0].auc.to_numpy()
    apart = table[table.gap == 4].auc.to_numpy()
    assert (apart > 0.95).all()
    assert (same < apart).all()
    assert (table.auc_low <= table.auc).all() and (table.auc <= table.auc_high).all()
    # independent seeds per structure and per point, recorded
    assert len(set(table.seed_a) | set(table.seed_b)) == 4
    assert all(c["multimodal"] == ["STED", "SMLM"] for c in fake_replicates)


def test_sweep_is_reproducible(fake_replicates):
    run = lambda: distinguishability_sweep(
        "a", "b", sweep={"gap": [2]}, modalities="STED", n_replicates=10,
        random_seed=3, n_bootstrap=50, random_state=0, progress=False,
    )
    pd.testing.assert_frame_equal(run(), run())


def test_arguments_merge_and_modality_parameters(fake_replicates):
    distinguishability_sweep(
        dict(structure="a", probe_DoL=2),
        dict(structure="b"),
        sweep={"SMLM.lateral_precision": [5, 10], "gap": [1]},
        modalities="SMLM",
        common=dict(probe_DoL=4, number_of_particles=1, SMLM={"nlocalisations": 3}),
        n_replicates=4, n_bootstrap=10, progress=False,
    )
    first_a, first_b = fake_replicates[0], fake_replicates[1]
    assert first_a["probe_DoL"] == 2 and first_b["probe_DoL"] == 4
    assert first_a["number_of_particles"] == 1
    assert first_a["SMLM"] == {"nlocalisations": 3, "lateral_precision": 5}
    assert fake_replicates[2]["SMLM"]["lateral_precision"] == 10


def test_scores_vs_parameter_series():
    rows = [
        dict(gap=g, le=le, modality=m, auc=0.5 + g / 10, auc_low=0.4, auc_high=1.0)
        for g in [2, 0] for le in [0.5, 1.0] for m in ["STED"]
    ]
    table = pd.DataFrame(rows)
    with pytest.raises(ValueError):
        scores_vs_parameter(table)
    data = scores_vs_parameter(table, "gap")
    assert list(data.series.unique()) == ["le=0.5", "le=1.0"]
    assert list(data.gap) == [0, 2, 0, 2]


def test_plot_returns_figure(fake_replicates):
    table, figure = distinguishability_sweep(
        "a", "b", sweep={"gap": [0, 3], "random_orientations": [False, True]},
        modalities="STED", n_replicates=6, n_bootstrap=10, progress=False,
        plot=True, plot_parameter="gap",
    )
    ax = figure.axes[0]
    assert ax.get_xlabel() == "gap"
    assert len(ax.get_lines()) == 2 + 1  # two series and the chance line
    figure = plot_distinguishability_sweep(table, "random_orientations", metric="accuracy")
    assert [t.get_text() for t in figure.axes[0].get_xticklabels()] == ["False", "True"]


def test_invalid_sweeps(fake_replicates):
    with pytest.raises(ValueError):
        distinguishability_sweep("a", "b", sweep={}, modalities="STED")
    with pytest.raises(ValueError):
        distinguishability_sweep("a", "b", sweep={"auc": [1]}, modalities="STED")


@pytest.mark.network
def test_sweep_runs_real_simulations():
    table = distinguishability_sweep(
        "7R5K", "7R5K",
        sweep={"labelling_efficiency": [0.5, 1.0]},
        modalities="Widefield",
        common=dict(probe_template="NPC_Nup96_Cterminal_direct", number_of_particles=1),
        n_replicates=3, random_seed=0, n_bootstrap=20, progress=False,
    )
    assert list(table.labelling_efficiency) == [0.5, 1.0]
    assert (table.n_a == 3).all() and (table.n_b == 3).all()
