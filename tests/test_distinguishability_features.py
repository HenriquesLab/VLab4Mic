"""Feature sets, sweep metrics and feature-table scoring."""
import numpy as np
import pandas as pd
import pytest
from scipy import ndimage

from vlab4mic.analysis.distinguishability import (
    feature_ablation,
    feature_metrics,
    feature_names,
    distinguishability_from_features,
    image_features,
    plot_distinguishability_heatmaps,
    same_field_of_view,
    whitened_features,
)


def particle(kind, rng, size=61, scale=1.0, stretch=1.0, angle=0.0, n=600, sigma=1.0):
    """Points in projection: "cap" crowds points towards the rim (a curved
    surface seen from above), "disc" spreads them uniformly."""
    r = np.sqrt(rng.uniform(0, 1, n)) if kind == "disc" else np.sqrt(1 - rng.uniform(0, 1, n) ** 2)
    phi = rng.uniform(0, 2 * np.pi, n)
    x, y = 12 * scale * stretch * r * np.cos(phi), 12 * scale * r * np.sin(phi)
    c, s = np.cos(angle), np.sin(angle)
    x, y = c * x - s * y + size / 2, s * x + c * y + size / 2
    image, _, _ = np.histogram2d(x, y, bins=size, range=[[0, size], [0, size]])
    return ndimage.gaussian_filter(image, sigma) + rng.normal(0, 0.01, image.shape)


def test_same_field_of_view_pads_and_crops():
    image = np.arange(25.0).reshape(5, 5)
    padded = same_field_of_view(image, pixelsize_nm=100, fov_nm=700)
    assert padded.shape == (7, 7)
    np.testing.assert_array_equal(padded[1:6, 1:6], image)
    assert padded[0, 0] == np.median(image)
    assert same_field_of_view(image, 100, 300).shape == (3, 3)


def test_whitened_features_ignore_size_elongation_and_orientation():
    rng = np.random.default_rng(0)
    base = np.mean([whitened_features(particle("cap", rng)) for _ in range(30)], axis=0)
    other = np.mean([whitened_features(particle("cap", rng, scale=1.4, stretch=1.6, angle=0.7))
                     for _ in range(30)], axis=0)
    disc = np.mean([whitened_features(particle("disc", rng)) for _ in range(30)], axis=0)
    # same structure at another size, elongation and orientation: close;
    # a different internal distribution: further apart
    assert np.abs(base - other)[:-1].sum() < 0.5 * np.abs(base - disc)[:-1].sum()


def test_whitened_features_flag_empty_images():
    features = whitened_features(np.zeros((9, 9)))
    assert features[-1] == 1 and np.all(features[:-1] == 0)
    assert len(feature_names("whitened")) == len(features)


def test_feature_metrics_match_image_features():
    rng = np.random.default_rng(1)
    image = particle("cap", rng)
    metrics = feature_metrics(["whitened", "brightness"], fov_nm=60)
    assert [m.__name__ for m in metrics] == feature_names("whitened") + ["brightness"]
    values = [m(simulated_image=image, simulated_image_pixelsize_nm=1.0) for m in metrics]
    expected = np.concatenate([
        image_features([image], "whitened", pixelsize_nm=1.0, fov_nm=60)[0],
        image_features([image], "brightness", pixelsize_nm=1.0, fov_nm=60)[0],
    ])
    np.testing.assert_allclose(values, expected)


def feature_table(rng, n=40):
    rows = []
    for condition, kinds in (("different", ("cap", "disc")), ("same", ("cap", "cap"))):
        for model, kind in zip(("A", "B"), kinds):
            for rep in range(n):
                image = particle(kind, rng, scale=rng.uniform(0.8, 1.2),
                                 stretch=rng.uniform(1, 1.5), angle=rng.uniform(0, np.pi))
                values = list(whitened_features(image)) + list(image_features([image], "brightness")[0])
                rows.append([condition, model, rep] + values)
    return pd.DataFrame(rows, columns=["condition", "model", "Replica"]
                        + feature_names("whitened") + ["brightness"])


@pytest.fixture(scope="module")
def table():
    return feature_table(np.random.default_rng(2))


def test_distinguishability_from_features_per_condition(table):
    scores = distinguishability_from_features(
        table, "model", ("A", "B"), "whitened", by="condition", controls={"brightness": "brightness"},
        replicate_column="Replica", n_bootstrap=100, random_state=0, progress=False,
    ).set_index("condition")
    assert scores.loc["different", "auc"] > 0.9
    assert abs(scores.loc["same", "auc"] - 0.5) < 0.2
    assert {"null_auc", "p_value", "brightness_auc", "n_a", "n_b"} <= set(scores.columns)


def test_parallel_scoring_matches_serial(table):
    kwargs = dict(by="condition", n_bootstrap=50, random_state=0, progress=False)
    serial = distinguishability_from_features(table, "model", ("A", "B"), "whitened", **kwargs)
    parallel = distinguishability_from_features(table, "model", ("A", "B"), "whitened",
                                                n_workers=2, **kwargs)
    pd.testing.assert_frame_equal(serial, parallel)


def test_feature_ablation_has_one_row_per_group_and_condition(table):
    ablation = feature_ablation(
        table, "model", ("A", "B"), by="condition", n_bootstrap=50, progress=False,
        groups={"profile": [f"whitened_{k:02d}" for k in range(8)],
                "moments": ["whitened_m4", "whitened_m6"]},
    )
    assert len(ablation) == 4
    assert set(ablation.feature_group) == {"profile", "moments"}


def test_heatmaps_draw_one_panel_per_facet_with_highlights():
    import matplotlib

    matplotlib.use("Agg")
    scores = pd.DataFrame([dict(method=m, le=le, si=si, auc=0.5 + 0.1 * k)
                           for m in ("a", "b") for k, le in enumerate([0.1, 0.5, 1.0])
                           for si in (0.5, 1.0)])
    figure = plot_distinguishability_heatmaps(
        scores, x="le", y="si", facet="method",
        highlight=[dict(facet="a", x=1.0, y=1.0, colour="red")],
    )
    panels = [ax for ax in figure.axes if ax.get_title()]
    assert [ax.get_title() for ax in panels] == ["a", "b"]
    assert len(panels[0].patches) == 1 and len(panels[1].patches) == 0
