"""Structural distinguishability: how well two candidate structures can be told apart.

Structural distinguishability is the accuracy with which a classifier assigns
simulated images of two candidate structures to the right class, under
stated simulation parameters. Each image is reduced to a feature vector, a
classifier is evaluated by stratified k-fold cross-validation, and the
out-of-fold scores give the classification accuracy and the area under the
receiver operating characteristic curve (ROC AUC), with bootstrap intervals
over realisations. An AUC of 0.5 means the structures cannot be told apart;
an AUC of 1 means every image is assigned correctly.

Typical use::

    results_a, _ = experiments.run_replicates(100, structure=dome, ...)
    results_b, _ = experiments.run_replicates(100, structure=flat, ...)
    score = distinguishability_from_replicates(results_a, results_b, "SMLM")

Across a parameter sweep::

    table = distinguishability_sweep(
        dict(structure=dome_path, structure_is_path=True),
        dict(structure=flat_path, structure_is_path=True),
        sweep={"labelling_efficiency": [0.2, 0.5, 1.0]},
        modalities=["STED", "SMLM"],
        common=dict(probe_list=[probe], number_of_particles=1),
        n_replicates=100,
        plot=True,
    )
"""

import copy
import itertools

import numpy as np
import pandas as pd
from scipy import ndimage
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from tqdm.auto import tqdm


def _as_2d(image):
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        # frames, x, y: integrate the frames
        image = image.sum(axis=0)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2-D image or a stack of frames, got {image.shape}")
    return image


def radial_profile_features(image, n_bins=16, max_radius=None):
    """
    Rotation-invariant features of a single-particle image.

    The image is centred on its intensity centroid. Features are the
    radial intensity profile (n_bins rings, normalised to total intensity),
    the radius of gyration, the two eigenvalues of the second-moment
    matrix (shape and elongation) and the total intensity.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames (frames, x, y) that is summed.
    n_bins : int, optional
        Number of radial bins. Default is 16.
    max_radius : float, optional
        Outer radius of the profile in pixels. Default: half the smaller
        image dimension.

    Returns
    -------
    numpy.ndarray
        Feature vector of length n_bins + 4.
    """
    image = np.clip(_as_2d(image), 0, None)
    total = image.sum()
    if total <= 0:
        return np.zeros(n_bins + 4)
    cx, cy = ndimage.center_of_mass(image)
    x, y = np.indices(image.shape)
    dx, dy = x - cx, y - cy
    r = np.hypot(dx, dy)
    if max_radius is None:
        max_radius = min(image.shape) / 2
    edges = np.linspace(0, max_radius, n_bins + 1)
    profile = np.histogram(r, bins=edges, weights=image)[0] / total
    weights = image / total
    gyration = np.sqrt((weights * r**2).sum())
    cov = np.array(
        [
            [(weights * dx * dx).sum(), (weights * dx * dy).sum()],
            [(weights * dx * dy).sum(), (weights * dy * dy).sum()],
        ]
    )
    eigenvalues = np.sort(np.linalg.eigvalsh(cov))
    return np.concatenate([profile, [gyration], eigenvalues, [np.log1p(total)]])


def pixel_features(image):
    """
    Normalised pixel intensities of an image (not rotation invariant).

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.

    Returns
    -------
    numpy.ndarray
        Flattened image divided by its total intensity.
    """
    image = _as_2d(image)
    total = image.sum()
    return (image / total if total > 0 else image).ravel()


def aspect_ratio_features(image, threshold_rel=0.2, smoothing_px=1.0):
    """
    Aspect ratio of a single-particle image (size and rotation invariant).

    The image is smoothed with a Gaussian of smoothing_px pixels (so that
    isolated noise pixels do not count), the background (median) is
    subtracted and pixels below threshold_rel of the maximum are ignored;
    the aspect ratio is the square root of the
    ratio of the largest to the smallest eigenvalue of the intensity
    second-moment matrix, i.e. the ratio of the long to the short principal
    axis of the particle seen in projection (1 for a circular particle).

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    threshold_rel : float, optional
        Fraction of the maximum below which pixels are ignored. Default 0.2.
    smoothing_px : float, optional
        Standard deviation of the smoothing Gaussian, in pixels. Default 1;
        0 disables smoothing.

    Returns
    -------
    numpy.ndarray
        One-element feature vector [aspect ratio].
    """
    image = _as_2d(image)
    if smoothing_px:
        image = ndimage.gaussian_filter(image, smoothing_px)
    image = np.clip(image - np.median(image), 0, None)
    if image.max() <= 0:
        return np.array([1.0])
    image = np.where(image >= threshold_rel * image.max(), image, 0)
    weights = image / image.sum()
    x, y = np.indices(image.shape)
    cx, cy = (weights * x).sum(), (weights * y).sum()
    dx, dy = x - cx, y - cy
    cov = np.array(
        [
            [(weights * dx * dx).sum(), (weights * dx * dy).sum()],
            [(weights * dx * dy).sum(), (weights * dy * dy).sum()],
        ]
    )
    low, high = np.sort(np.linalg.eigvalsh(cov))
    if low <= 0:
        return np.array([np.inf if high > 0 else 1.0])
    return np.array([np.sqrt(high / low)])


FEATURES = {
    "radial_profile": radial_profile_features,
    "pixels": pixel_features,
    "aspect_ratio": aspect_ratio_features,
}


def image_features(images, features="radial_profile", **kwargs):
    """
    Feature matrix of a list of images.

    Parameters
    ----------
    images : sequence of numpy.ndarray
        Images (2-D, or stacks of frames that are summed).
    features : str or callable, optional
        "radial_profile" (default, rotation invariant), "aspect_ratio"
        (rotation and size invariant shape), "pixels", or a function
        image -> 1-D feature vector.
    **kwargs
        Passed to the feature function.

    Returns
    -------
    numpy.ndarray
        Array of shape (n_images, n_features).
    """
    function = FEATURES[features] if isinstance(features, str) else features
    return np.array([function(image, **kwargs) for image in images])


def _bootstrap_interval(labels, scores, metric, n_bootstrap, confidence, rng):
    values = []
    index_a = np.flatnonzero(labels == 0)
    index_b = np.flatnonzero(labels == 1)
    for _ in range(n_bootstrap):
        sample = np.concatenate(
            [rng.choice(index_a, index_a.size), rng.choice(index_b, index_b.size)]
        )
        values.append(metric(labels[sample], scores[sample]))
    alpha = (1 - confidence) / 2
    return tuple(np.quantile(values, [alpha, 1 - alpha]))


def _accuracy(labels, scores):
    return float(np.mean((scores >= 0.5) == labels))


def distinguishability(
    images_a,
    images_b,
    features="radial_profile",
    n_splits=5,
    n_bootstrap=1000,
    confidence=0.95,
    classifier=None,
    random_state=None,
    **feature_kwargs,
):
    """
    Accuracy and ROC AUC with which two sets of images are told apart.

    Each image is reduced to a feature vector (image_features). A
    classifier (default: standardised logistic regression) is evaluated by
    stratified k-fold cross-validation, so every image is scored by a model
    that did not see it. Accuracy (threshold 0.5) and ROC AUC are computed
    from these out-of-fold scores, and their intervals from bootstrap
    resamples of the realisations within each class.

    Parameters
    ----------
    images_a, images_b : sequence of numpy.ndarray
        Simulated images of structure A and of structure B, one per
        independent realisation (e.g. from run_replicates).
    features : str or callable, optional
        Feature extraction (see image_features). Default "radial_profile".
    n_splits : int, optional
        Number of cross-validation folds. Default is 5 (reduced if a class
        has fewer images).
    n_bootstrap : int, optional
        Number of bootstrap resamples for the intervals. Default is 1000.
    confidence : float, optional
        Confidence level of the intervals. Default is 0.95.
    classifier : sklearn estimator, optional
        Classifier with predict_proba. Default: StandardScaler +
        LogisticRegression.
    random_state : int, optional
        Seed for the folds and the bootstrap.
    **feature_kwargs
        Passed to the feature function.

    Returns
    -------
    dict
        accuracy, accuracy_interval, auc, auc_interval, scores (out-of-fold
        probability of class B per image), labels (0 for A, 1 for B),
        n_a, n_b, features (name or function), n_splits.
    """
    features_a = image_features(images_a, features, **feature_kwargs)
    features_b = image_features(images_b, features, **feature_kwargs)
    X = np.vstack([features_a, features_b])
    labels = np.concatenate([np.zeros(len(features_a)), np.ones(len(features_b))]).astype(int)
    n_splits = int(min(n_splits, len(features_a), len(features_b)))
    if n_splits < 2:
        raise ValueError("Each structure needs at least 2 images")
    if classifier is None:
        classifier = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    scores = cross_val_predict(classifier, X, labels, cv=folds, method="predict_proba")[:, 1]
    rng = np.random.default_rng(random_state)
    return dict(
        accuracy=_accuracy(labels, scores),
        accuracy_interval=_bootstrap_interval(labels, scores, _accuracy, n_bootstrap, confidence, rng),
        auc=float(roc_auc_score(labels, scores)),
        auc_interval=_bootstrap_interval(labels, scores, roc_auc_score, n_bootstrap, confidence, rng),
        scores=scores,
        labels=labels,
        n_a=len(features_a),
        n_b=len(features_b),
        features=features,
        n_splits=n_splits,
    )


def images_from_replicates(results, modality, channel="ch0", noiseless=False):
    """
    Images of one modality and channel from run_replicates results.

    Parameters
    ----------
    results : list of dict
        Output of run_replicates.
    modality : str
        Modality name.
    channel : str, optional
        Channel name. Default "ch0".
    noiseless : bool, optional
        Use the images without detector noise. Default False.

    Returns
    -------
    list of numpy.ndarray
    """
    key = "images_noiseless" if noiseless else "images"
    return [np.asarray(r[key][modality][channel]) for r in results]


def distinguishability_from_replicates(
    results_a, results_b, modality, channel="ch0", noiseless=False, **kwargs
):
    """
    Distinguishability of two structures from run_replicates results.

    Parameters
    ----------
    results_a, results_b : list of dict
        run_replicates output for structure A and structure B, with one
        particle per realisation.
    modality : str
        Modality to compare.
    channel : str, optional
        Channel name. Default "ch0".
    noiseless : bool, optional
        Use images without detector noise. Default False.
    **kwargs
        Passed to distinguishability.

    Returns
    -------
    dict
        See distinguishability.
    """
    return distinguishability(
        images_from_replicates(results_a, modality, channel, noiseless),
        images_from_replicates(results_b, modality, channel, noiseless),
        **kwargs,
    )

def distinguishability_from_dataframes(
        features_a, 
        features_b,
        n_splits=5,
        n_bootstrap=1000,
        confidence=0.95,
        classifier=None,
        random_state=None,
        ):
    """
    features_a is a list of features measured on model a
    features_b is a list of features measured on model b
    """
    X = np.vstack([features_a, features_b])
    labels = np.concatenate([np.zeros(len(features_a)), np.ones(len(features_b))]).astype(int)
    n_splits = int(min(n_splits, len(features_a), len(features_b)))
    if n_splits < 2:
        raise ValueError("Each structure needs at least 2 images")
    if classifier is None:
        classifier = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    scores = cross_val_predict(classifier, X, labels, cv=folds, method="predict_proba")[:, 1]
    rng = np.random.default_rng(random_state)
    return dict(
        accuracy=_accuracy(labels, scores),
        accuracy_interval=_bootstrap_interval(labels, scores, _accuracy, n_bootstrap, confidence, rng),
        auc=float(roc_auc_score(labels, scores)),
        auc_interval=_bootstrap_interval(labels, scores, roc_auc_score, n_bootstrap, confidence, rng),
        scores=scores,
        labels=labels,
        n_a=len(features_a),
        n_b=len(features_b),
        n_splits=n_splits,
    )


# columns of a distinguishability_sweep table that are not swept parameters
SWEEP_RESULT_COLUMNS = (
    "modality",
    "channel",
    "n_a",
    "n_b",
    "accuracy",
    "accuracy_low",
    "accuracy_high",
    "auc",
    "auc_low",
    "auc_high",
    "seed_a",
    "seed_b",
)


def _as_structure_kwargs(structure):
    if isinstance(structure, str):
        return dict(structure=structure)
    if isinstance(structure, dict):
        return dict(structure)
    raise TypeError("A structure is a structure ID / path or a dict of image_vsample arguments")


def _set_parameter(kwargs, name, value):
    # "Modality.parameter" sets a parameter of one modality, the way
    # image_vsample takes them: Modality={"parameter": value}
    if "." in name:
        modality, key = name.split(".", 1)
        kwargs[modality] = {**kwargs.get(modality, {}), key: value}
    else:
        kwargs[name] = value


def _point_seeds(random_seed, n_points):
    if random_seed is None:
        return [(None, None)] * n_points
    state = np.random.SeedSequence(random_seed).generate_state(2 * n_points)
    return [(int(state[2 * k]), int(state[2 * k + 1])) for k in range(n_points)]


def distinguishability_sweep(
    structure_a,
    structure_b,
    sweep,
    modalities,
    n_replicates=50,
    common=None,
    channel="ch0",
    noiseless=False,
    random_seed=None,
    plot=False,
    plot_parameter=None,
    metric="auc",
    progress=True,
    **kwargs,
):
    """
    Distinguishability of two structures at every point of a parameter sweep.

    At each combination of the swept parameter values, n_replicates
    independent realisations of each structure are simulated with
    run_replicates (a new experiment per structure and per point), and
    distinguishability is computed for every modality.

    Parameters
    ----------
    structure_a, structure_b : str or dict
        The two structures: a structure ID, or a dict of image_vsample
        arguments specific to that structure (e.g. structure,
        structure_is_path, probe_list). These override common.
    sweep : dict
        {parameter: list of values}. Every combination of values (Cartesian
        product) is a sweep point. Names are image_vsample arguments (e.g.
        "labelling_efficiency", "probe_DoL", "structural_integrity",
        "random_orientations") or "Modality.parameter" for a modality
        parameter (e.g. "SMLM.lateral_precision", "STED.pixelsize_nm").
        A swept value overrides the same argument in common and in the
        structure arguments.
    modalities : str or list of str
        Modalities to simulate and compare.
    n_replicates : int, optional
        Realisations per structure and per point. Default is 50.
    common : dict, optional
        image_vsample arguments shared by both structures (probes, virtual
        sample, modality parameters...).
    channel : str, optional
        Channel to compare. Default "ch0".
    noiseless : bool, optional
        Use images without detector noise. Default False.
    random_seed : int, optional
        Seed of the whole sweep. Each structure at each point gets its own
        seed derived from it (recorded in the table as seed_a, seed_b), so
        the sweep is reproducible and realisations are independent across
        points and structures.
    plot : bool, optional
        Also plot the score against a swept parameter (see
        plot_distinguishability_sweep). Default False.
    plot_parameter : str, optional
        Swept parameter on the x axis of the plot. Needed when more than one
        parameter is swept.
    metric : str, optional
        "auc" (default) or "accuracy", for the plot.
    progress : bool, optional
        Show a progress bar over the sweep points. Default True.
    **kwargs
        Passed to distinguishability (features, n_splits, n_bootstrap,
        confidence, classifier, random_state and feature arguments).

    Returns
    -------
    pandas.DataFrame
        One row per sweep point and modality: the swept parameter values,
        modality, channel, n_a, n_b, accuracy, accuracy_low, accuracy_high,
        auc, auc_low, auc_high, seed_a, seed_b.
    matplotlib.figure.Figure
        Only if plot is True, returned as (table, figure).
    """
    from vlab4mic import experiments

    structure_a = _as_structure_kwargs(structure_a)
    structure_b = _as_structure_kwargs(structure_b)
    modalities = [modalities] if isinstance(modalities, str) else list(modalities)
    if not sweep:
        raise ValueError("sweep needs at least one parameter with a list of values")
    names = list(sweep)
    for name in names:
        if name in SWEEP_RESULT_COLUMNS:
            raise ValueError(f"{name!r} cannot be a swept parameter name")
    points = list(itertools.product(*(list(sweep[name]) for name in names)))
    seeds = _point_seeds(random_seed, len(points))

    rows = []
    for point, point_seeds in tqdm(
        list(zip(points, seeds)), desc="sweep points", disable=not progress
    ):
        setting = dict(zip(names, point))
        replicates = []
        for structure, seed in zip((structure_a, structure_b), point_seeds):
            arguments = copy.deepcopy({**(common or {}), **structure})
            for name, value in setting.items():
                _set_parameter(arguments, name, value)
            arguments.update(
                multimodal=modalities,
                random_seed=seed,
                n_replicates=n_replicates,
                modality="All",
            )
            results, _ = experiments.run_replicates(**arguments)
            replicates.append(results)
        for modality in modalities:
            score = distinguishability_from_replicates(
                replicates[0], replicates[1], modality, channel, noiseless, **kwargs
            )
            rows.append(
                dict(
                    setting,
                    modality=modality,
                    channel=channel,
                    n_a=score["n_a"],
                    n_b=score["n_b"],
                    accuracy=score["accuracy"],
                    accuracy_low=score["accuracy_interval"][0],
                    accuracy_high=score["accuracy_interval"][1],
                    auc=score["auc"],
                    auc_low=score["auc_interval"][0],
                    auc_high=score["auc_interval"][1],
                    seed_a=point_seeds[0],
                    seed_b=point_seeds[1],
                )
            )
    table = pd.DataFrame(rows)
    if plot:
        figure = plot_distinguishability_sweep(table, plot_parameter, metric=metric)
        return table, figure
    return table


def swept_parameters(table):
    """
    Names of the swept parameters in a distinguishability_sweep table.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of distinguishability_sweep (also after a CSV round trip).

    Returns
    -------
    list of str
    """
    return [c for c in table.columns if c not in SWEEP_RESULT_COLUMNS]


def scores_vs_parameter(table, parameter=None, metric="auc"):
    """
    Score against one swept parameter, one series per curve.

    There is one series per modality and per combination of the values of
    the other swept parameters.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of distinguishability_sweep.
    parameter : str, optional
        Swept parameter for the x axis. Default: the only swept parameter.
    metric : str, optional
        "auc" (default) or "accuracy".

    Returns
    -------
    pandas.DataFrame
        Columns series, the parameter, the metric and its interval
        (metric_low, metric_high), sorted by series and parameter value.
    """
    if metric not in ("auc", "accuracy"):
        raise ValueError("metric is 'auc' or 'accuracy'")
    swept = swept_parameters(table)
    if parameter is None:
        if len(swept) != 1:
            raise ValueError(f"Several parameters were swept {swept}: choose one")
        parameter = swept[0]
    if parameter not in swept:
        raise ValueError(f"{parameter!r} is not a swept parameter {swept}")
    others = [p for p in swept if p != parameter]
    name_modality = table["modality"].nunique() > 1 or not others

    def label(row):
        parts = [str(row["modality"])] if name_modality else []
        parts += [f"{p}={row[p]}" for p in others]
        return ", ".join(parts)

    out = table[[parameter, metric, f"{metric}_low", f"{metric}_high"]].copy()
    out.insert(0, "series", table.apply(label, axis=1))
    return out.sort_values(["series", parameter], kind="stable").reset_index(drop=True)


def plot_distinguishability_sweep(table, parameter=None, metric="auc", ax=None):
    """
    Plot the distinguishability score against one swept parameter.

    One line per modality and per combination of the other swept
    parameters, with the bootstrap interval as a band and a dashed line at
    chance (0.5).

    Parameters
    ----------
    table : pandas.DataFrame
        Output of distinguishability_sweep.
    parameter : str, optional
        Swept parameter for the x axis. Default: the only swept parameter.
    metric : str, optional
        "auc" (default) or "accuracy".
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. Default: a new figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    data = scores_vs_parameter(table, parameter, metric)
    parameter = data.columns[1]
    if ax is None:
        figure, ax = plt.subplots(figsize=(5, 3.5))
    else:
        figure = ax.figure
    x_values = pd.to_numeric(data[parameter], errors="coerce")
    numeric = data[parameter].dtype != bool and not x_values.isna().any()
    if not numeric:
        # categorical values (booleans, strings, lists): evenly spaced
        labels = data[parameter].astype(str)
        categories = list(dict.fromkeys(labels))
        x_values = labels.map(categories.index)
    for k, (series, group) in enumerate(data.groupby("series", sort=False)):
        x = x_values[group.index]
        colour = f"C{k % 10}"
        ax.fill_between(x, group[f"{metric}_low"], group[f"{metric}_high"],
                        color=colour, alpha=0.2, linewidth=0)
        ax.plot(x, group[metric], "o-", color=colour, linewidth=2, markersize=5,
                label=series)
    if not numeric:
        ax.set_xticks(range(len(categories)), categories)
    ax.axhline(0.5, color="0.5", linestyle="--", linewidth=1)
    ax.set_ylim(min(0.4, data[f"{metric}_low"].min() - 0.02), 1.02)
    ax.set_xlabel(parameter)
    ax.set_ylabel("ROC AUC" if metric == "auc" else "accuracy")
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    if data["series"].nunique() > 1:
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.set_title(data["series"].iloc[0], fontsize=10)
    figure.tight_layout()
    return figure


def group_dataframes_1value(dataframe1, dataframe2,
                             main_category="modality_name",
                             param_name_1="labelling_efficiency",
                             metric_name = "image_aspect_ratio"):
    rows = []
    main_category_values = list(dataframe1[main_category].unique())
    for category_value in main_category_values:
        for param_1_value in list(dataframe1[param_name_1].unique()):
            features_tmp_1 = dataframe1.loc[
                (dataframe1[main_category] == category_value) & (dataframe1[param_name_1] == param_1_value),
                metric_name]
            features_tmp_2 = dataframe2.loc[
                (dataframe2[main_category] == category_value) & (dataframe2[param_name_1] == param_1_value),
                metric_name]
            score = distinguishability_from_dataframes(
                np.reshape(features_tmp_1, -1).reshape(-1,1), 
                np.reshape(features_tmp_2, -1).reshape(-1,1))
            tmp = dict(
                    modality=category_value,
                    metric_used = metric_name,
                    n_a=score["n_a"],
                    n_b=score["n_b"],
                    accuracy=score["accuracy"],
                    accuracy_low=score["accuracy_interval"][0],
                    accuracy_high=score["accuracy_interval"][1],
                    auc=score["auc"],
                    auc_low=score["auc_interval"][0],
                    auc_high=score["auc_interval"][1],
                )
            tmp[param_name_1] = param_1_value
            rows.append(tmp)
    return pd.DataFrame(rows)