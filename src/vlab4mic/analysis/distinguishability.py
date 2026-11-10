"""Structural distinguishability: how well two candidate structures can be told apart.

Structural distinguishability is the accuracy with which a classifier assigns
simulated images of two candidate structures to the right class, under
stated simulation parameters. Each image is reduced to a feature vector, a
classifier is evaluated by stratified k-fold cross-validation, and the
out-of-fold scores give the classification accuracy and the area under the
receiver operating characteristic curve (ROC AUC), with bootstrap intervals
over realisations. An AUC of 0.5 means the structures cannot be told apart;
an AUC of 1 means every image is assigned correctly.

The default features ("shape") describe the shape of the particle image
only: each image is background subtracted and scaled to unit total
intensity before the features are computed, so two structures that differ
only in brightness (e.g. in the number of labelled sites) are not told
apart. Brightness can be scored on its own with features="brightness", as a
control. The default classifier is standardisation, principal component
analysis (95% of the variance) and logistic regression.

Two controls are reported with every score:

- null_auc: the same analysis on two random halves of structure A, which
  should be close to 0.5. A value far from 0.5 means the score is biased.
- p_value: the permutation p-value of the AUC (n_permutations > 0), the
  fraction of label permutations that score at least as well.

A warning is raised when all realisations of a structure have the same
features, which happens when the realisations are not independent (e.g.
every replicate simulated with the same random seed); the score is then
meaningless.

Choose features that measure the difference in question and are blind to
the rest. "whitened" features are invariant to size, elongation and
orientation, and compare the internal distribution of labels (e.g. the
curvature of a lattice); "shape" features include size and outline. Simulate
the variability that real samples have (random orientation, and random
position relative to the pixels with position_jitter_nm), so the classifier
cannot learn features of the simulation that no experiment reproduces.

With a parameter sweep, features of every image are computed as sweep
metrics and scored afterwards from the feature table::

    metrics = feature_metrics("whitened", fov_nm=600)
    tables = []
    for name, path in {"dome": dome_path, "flat": flat_path}.items():
        sweep = run_parameter_sweep(structures=[path], structure_is_path=True,
                                    use_reference=False, custom_metrics=metrics,
                                    default_metrics=[], position_jitter_nm=100,
                                    labelling_efficiency=[0.1, 0.5, 1.0], ...,
                                    return_generator=True)
        tables.append(sweep.analysis["dataframes"].assign(model=name))
    features = pd.concat(tables)
    scores = distinguishability_from_features(
        features, "model", ("dome", "flat"), "whitened",
        by=["modality_name", "labelling_efficiency"], n_permutations=200)
    plot_distinguishability_heatmaps(scores, x="labelling_efficiency", ...)

Other typical use::

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
import inspect
import itertools
import warnings

import numpy as np
import pandas as pd
from scipy import ndimage
from sklearn.decomposition import PCA
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


def normalise_image(image, threshold_rel=0.1):
    """
    Background-subtracted image scaled to unit total intensity.

    The median is subtracted, pixels below threshold_rel of the maximum are
    set to 0 and the image is divided by its total. Multiplying the signal
    by a constant (a brighter particle) does not change the result. The
    threshold removes the background noise, which would otherwise make up a
    larger share of the total in a dimmer image and so let brightness leak
    into the shape features.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames (frames, x, y) that is summed.
    threshold_rel : float, optional
        Fraction of the maximum below which pixels are set to 0. Default
        0.1; 0 keeps every pixel above the median.

    Returns
    -------
    numpy.ndarray
        2-D image with non-negative values summing to 1 (or all zeros).
    """
    image = _as_2d(image)
    image = np.clip(image - np.median(image), 0, None)
    image = np.where(image >= threshold_rel * image.max(), image, 0)
    total = image.sum()
    return image / total if total > 0 else image


def same_field_of_view(image, pixelsize_nm, fov_nm):
    """
    Pad or crop an image, about its centre, to a field of view in nm.

    The image size of a simulation can differ between conditions or
    realisations of the same imaging method; features that depend on the
    image size (e.g. rings out to half the image) are only comparable on a
    common field of view. Padding uses the median of the image.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    pixelsize_nm : float
        Pixel size of the image in nm.
    fov_nm : float
        Side of the square output field of view in nm.

    Returns
    -------
    numpy.ndarray
        Square image of round(fov_nm / pixelsize_nm) pixels per side.
    """
    image = _as_2d(image)
    n = max(int(round(fov_nm / pixelsize_nm)), 1)
    out = np.full((n, n), np.median(image))
    src, dst = [], []
    for size in image.shape:
        if size >= n:
            start = (size - n) // 2
            src.append(slice(start, start + n))
            dst.append(slice(0, n))
        else:
            start = (n - size) // 2
            src.append(slice(0, size))
            dst.append(slice(start, start + size))
    out[tuple(dst)] = image[tuple(src)]
    return out


def prepare_image(image, pixelsize_nm=None, fov_nm=None, threshold_rel=0.1):
    """
    Common field of view (optional), background subtraction and unit total.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    pixelsize_nm, fov_nm : float, optional
        If both are given, the image is first brought to a field of view of
        fov_nm (see same_field_of_view).
    threshold_rel : float, optional
        Background threshold, see normalise_image. Default 0.1.

    Returns
    -------
    numpy.ndarray
        2-D image with non-negative values summing to 1 (or all zeros).
    """
    if pixelsize_nm is not None and fov_nm is not None:
        image = same_field_of_view(image, pixelsize_nm, fov_nm)
    return normalise_image(image, threshold_rel)


def shape_features(image, n_bins=8, threshold_rel=0.1):
    """
    Rotation-invariant, brightness-independent shape features.

    The image is normalised (normalise_image) and centred on its intensity
    centroid. Features are:

    - the radial intensity profile (n_bins rings out to half the smaller
      image dimension);
    - the log of the rotation-averaged power spectrum (n_bins rings out to
      the Nyquist frequency), which separates rings from filled discs;
    - the radius of gyration, the lengths of the major and minor principal
      axes (square roots of the eigenvalues of the second-moment matrix)
      and their ratio, minor / major (1 for a round particle).

    Lengths are in pixels, so compare images of the same modality.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    n_bins : int, optional
        Number of rings of the profile and of the spectrum. Default 8.
    threshold_rel : float, optional
        Background threshold, see normalise_image. Default 0.1.

    Returns
    -------
    numpy.ndarray
        Feature vector of length 2 * n_bins + 4.
    """
    image = normalise_image(image, threshold_rel)
    if image.sum() <= 0:
        return np.zeros(2 * n_bins + 4)
    x, y = np.indices(image.shape)
    dx = x - (image * x).sum()
    dy = y - (image * y).sum()
    r = np.hypot(dx, dy)
    profile = np.histogram(
        r, bins=np.linspace(0, min(image.shape) / 2, n_bins + 1), weights=image
    )[0]
    # rotation-averaged power spectrum; frequency in units of Nyquist
    power = np.abs(np.fft.fftshift(np.fft.fft2(image))) ** 2
    centre = np.array(image.shape) // 2
    k = np.hypot(x - centre[0], y - centre[1]) / (min(image.shape) / 2)
    edges = np.linspace(0, 1, n_bins + 1)
    sums = np.histogram(k, bins=edges, weights=power)[0]
    counts = np.histogram(k, bins=edges)[0]
    spectrum = np.log(sums / np.maximum(counts, 1) + 1e-12)
    gyration = np.sqrt((image * r**2).sum())
    cov = np.array(
        [
            [(image * dx * dx).sum(), (image * dx * dy).sum()],
            [(image * dx * dy).sum(), (image * dy * dy).sum()],
        ]
    )
    minor, major = np.sqrt(np.clip(np.sort(np.linalg.eigvalsh(cov)), 0, None))
    roundness = minor / major if major > 0 else 1.0
    return np.concatenate([profile, spectrum, [gyration, major, minor, roundness]])


def brightness_features(image):
    """
    Total intensity of an image above background (a control feature).

    Use it to measure how much of a difference between two structures is
    brightness alone (e.g. a different number of labelled sites), which the
    shape features do not see.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.

    Returns
    -------
    numpy.ndarray
        One-element feature vector [log(1 + total intensity above the median)].
    """
    image = _as_2d(image)
    return np.array([np.log1p(np.clip(image - np.median(image), 0, None).sum())])


def whitened_features(image, n_bins=8, rho_max=3.0, threshold_rel=0.1):
    """
    Size-, elongation- and orientation-invariant ("curvature") features.

    The image is normalised (normalise_image) and whitened: coordinates are
    centred on the intensity centroid, rotated to the principal axes and
    divided by the standard deviation along each axis, so every particle
    has unit spread along both axes. What remains is how intensity is
    distributed within that normalised outline, e.g. labels crowded towards
    the rim of a curved cap seen from above:

    - the radial density profile in whitened units (n_bins rings out to
      rho_max, intensity per unit area);
    - the moments <rho^4> and <rho^6> (<rho^2> = 2 by construction);
    - an indicator of images without a 2-D signal (all zeros, a single
      pixel or a line), for which the other features are 0.

    Use it to compare structures by their internal distribution of labels
    rather than by their size or outline.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    n_bins : int, optional
        Number of rings. Default 8.
    rho_max : float, optional
        Outer radius of the profile in whitened units. Default 3.
    threshold_rel : float, optional
        Background threshold, see normalise_image. Default 0.1.

    Returns
    -------
    numpy.ndarray
        Feature vector of length n_bins + 3.
    """
    image = normalise_image(image, threshold_rel)
    empty = np.zeros(n_bins + 3)
    empty[-1] = 1.0
    if image.sum() <= 0:
        return empty
    x, y = np.indices(image.shape)
    dx = x - (image * x).sum()
    dy = y - (image * y).sum()
    cov = np.array(
        [
            [(image * dx * dx).sum(), (image * dx * dy).sum()],
            [(image * dx * dy).sum(), (image * dy * dy).sum()],
        ]
    )
    values, vectors = np.linalg.eigh(cov)
    if values.min() <= 1e-9:
        return empty
    u = (dx * vectors[0, 1] + dy * vectors[1, 1]) / np.sqrt(values[1])
    v = (dx * vectors[0, 0] + dy * vectors[1, 0]) / np.sqrt(values[0])
    rho = np.hypot(u, v)
    edges = np.linspace(0, rho_max, n_bins + 1)
    mass = np.histogram(rho, bins=edges, weights=image)[0]
    density = mass / (np.pi * (edges[1:] ** 2 - edges[:-1] ** 2))
    moments = [(image * rho**4).sum(), (image * rho**6).sum()]
    return np.concatenate([density, moments, [0.0]])


def radial_profile_features(image, n_bins=16, max_radius=None):
    """
    Rotation-invariant features of a single-particle image.

    The image is centred on its intensity centroid. Features are the
    radial intensity profile (n_bins rings, normalised to total intensity),
    the radius of gyration, the two eigenvalues of the second-moment
    matrix (shape and elongation) and the total intensity. Because of the
    total intensity, and because the background is not subtracted, these
    features depend on brightness; use "shape" to compare shapes only.

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
    "shape": shape_features,
    "whitened": whitened_features,
    "brightness": brightness_features,
    "radial_profile": radial_profile_features,
    "pixels": pixel_features,
    "aspect_ratio": aspect_ratio_features,
}


def _call_feature(function, image, kwargs):
    # pass a feature function only the arguments it accepts, so one set of
    # arguments (e.g. n_bins) can serve several feature sets
    accepted = inspect.signature(function).parameters
    return function(image, **{k: v for k, v in kwargs.items() if k in accepted})


def feature_names(features, **kwargs):
    """
    Column names of a named feature set, as used in feature tables.

    Parameters
    ----------
    features : str
        "shape", "whitened", "brightness" or "aspect_ratio" (sets with a
        fixed number of features).
    **kwargs
        Feature arguments that change the number of features (n_bins).

    Returns
    -------
    list of str
    """
    if features == "shape":
        n = 2 * kwargs.get("n_bins", 8) + 4
        return [f"shape_{k:02d}" for k in range(n)]
    if features == "whitened":
        n_bins = kwargs.get("n_bins", 8)
        return ([f"whitened_{k:02d}" for k in range(n_bins)]
                + ["whitened_m4", "whitened_m6", "whitened_empty"])
    if features in ("brightness", "aspect_ratio"):
        return [features]
    raise ValueError(f"No fixed column names for features {features!r}")


def image_features(images, features="shape", pixelsize_nm=None, fov_nm=None, **kwargs):
    """
    Feature matrix of a list of images.

    Parameters
    ----------
    images : sequence of numpy.ndarray
        Images (2-D, or stacks of frames that are summed).
    features : str or callable, optional
        "shape" (default: rotation invariant, brightness independent),
        "whitened" (curvature: size, elongation and orientation invariant),
        "brightness" (total intensity, a control), "aspect_ratio"
        (rotation, size and brightness invariant), "radial_profile"
        (rotation invariant, includes brightness), "pixels", or a function
        image -> 1-D feature vector.
    pixelsize_nm, fov_nm : float, optional
        If both are given, every image is first brought to a field of view
        of fov_nm (see same_field_of_view), so that features of images of
        different sizes are comparable.
    **kwargs
        Passed to the feature function.

    Returns
    -------
    numpy.ndarray
        Array of shape (n_images, n_features).
    """
    function = FEATURES[features] if isinstance(features, str) else features
    if pixelsize_nm is not None and fov_nm is not None:
        images = [same_field_of_view(image, pixelsize_nm, fov_nm) for image in images]
    return np.array([function(image, **kwargs) for image in images])


def feature_metrics(features, fov_nm=None, **kwargs):
    """
    Features as custom metrics of a parameter sweep.

    Returns one metric function per feature column, to pass as
    custom_metrics to run_parameter_sweep (or add_custom_analysis_metrics).
    The sweep analysis dataframe then has one row per simulated image and
    one column per feature: a feature table for
    distinguishability_from_features. Each image is measured once, whatever
    the number of columns.

    Parameters
    ----------
    features : str or list of str
        Named feature sets with fixed columns (see feature_names), e.g.
        ["whitened"] or ["whitened", "brightness"].
    fov_nm : float, optional
        Common field of view in nm (see same_field_of_view), using the pixel
        size of each image. Recommended: image sizes can differ between
        conditions and realisations.
    **kwargs
        Passed to the feature functions.

    Returns
    -------
    list of callable
        Metric functions, named after their feature column.
    """
    features = [features] if isinstance(features, str) else list(features)
    columns = [(name, feature_names(name, **kwargs)) for name in features]
    last = {"image": None, "values": None}

    def measure(image, pixelsize_nm):
        # the sweep calls every metric in turn with the same image object;
        # keeping a reference to it means its id cannot be reused
        if last["image"] is not image:
            prepared = image
            if fov_nm is not None:
                prepared = same_field_of_view(image, pixelsize_nm, fov_nm)
            last["image"] = image
            last["values"] = np.concatenate(
                [np.atleast_1d(_call_feature(FEATURES[name], prepared, kwargs))
                 for name, _ in columns]
            )
        return last["values"]

    metrics = []
    index = 0
    for _, names in columns:
        for name in names:
            def metric(simulated_image=None, simulated_image_pixelsize_nm=None,
                       _index=index, **unused):
                return float(measure(simulated_image, simulated_image_pixelsize_nm)[_index])
            metric.__name__ = name
            metrics.append(metric)
            index += 1
    return metrics


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


def default_classifier(n_features):
    """
    Default classifier: standardisation, PCA and logistic regression.

    PCA keeps the components that explain 95% of the variance (only with
    more than two features), so correlated features (e.g. neighbouring
    rings of a profile) do not make the logistic regression unstable.

    Parameters
    ----------
    n_features : int
        Number of features.

    Returns
    -------
    sklearn.pipeline.Pipeline
    """
    steps = [StandardScaler()]
    if n_features > 2:
        steps.append(PCA(n_components=0.95, svd_solver="full"))
    steps.append(LogisticRegression(max_iter=1000))
    return make_pipeline(*steps)


def _out_of_fold_scores(X, labels, classifier, n_splits, random_state):
    n_splits = int(min(n_splits, np.bincount(labels).min()))
    if n_splits < 2:
        raise ValueError("Each structure needs at least 2 images")
    if np.all(np.ptp(X, axis=0) == 0):
        # every image has the same features: nothing to learn
        return np.full(len(labels), 0.5), n_splits
    if classifier is None:
        classifier = default_classifier(X.shape[1])
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    scores = cross_val_predict(classifier, X, labels, cv=folds, method="predict_proba")[:, 1]
    return scores, n_splits


def _as_feature_matrix(features):
    features = np.asarray(features, dtype=float)
    return features.reshape(len(features), -1)


def _score_features(
    features_a,
    features_b,
    n_splits=5,
    n_bootstrap=1000,
    confidence=0.95,
    classifier=None,
    random_state=None,
    n_permutations=0,
    null_control=True,
):
    # cross-validated scores of two feature matrices, with intervals and
    # controls (see distinguishability)
    features_a = _as_feature_matrix(features_a)
    features_b = _as_feature_matrix(features_b)
    n_unique = []
    for name, feats in (("A", features_a), ("B", features_b)):
        unique = len(np.unique(feats, axis=0))
        n_unique.append(unique)
        if len(feats) > 1 and unique < 2:
            warnings.warn(
                f"All {len(feats)} realisations of structure {name} have the "
                "same features: they are not independent (were they simulated "
                "with the same random seed?), so the score is not meaningful.",
                RuntimeWarning,
                stacklevel=3,
            )
    X = np.vstack([features_a, features_b])
    labels = np.concatenate([np.zeros(len(features_a)), np.ones(len(features_b))]).astype(int)
    scores, n_splits = _out_of_fold_scores(X, labels, classifier, n_splits, random_state)
    auc = float(roc_auc_score(labels, scores))
    rng = np.random.default_rng(random_state)
    p_value = np.nan
    if n_permutations:
        permuted = []
        for _ in range(int(n_permutations)):
            shuffled = rng.permutation(labels)
            perm_scores, _ = _out_of_fold_scores(X, shuffled, classifier, n_splits, random_state)
            permuted.append(roc_auc_score(shuffled, perm_scores))
        p_value = (1 + np.sum(np.array(permuted) >= auc)) / (1 + int(n_permutations))
    null_auc = np.nan
    if null_control and len(features_a) >= 4:
        # two random halves of structure A: should not be told apart
        order = rng.permutation(len(features_a))
        half = len(order) // 2
        null_labels = np.r_[np.zeros(half), np.ones(len(order) - half)].astype(int)
        null_scores, _ = _out_of_fold_scores(
            features_a[order], null_labels, classifier, n_splits, random_state
        )
        null_auc = float(roc_auc_score(null_labels, null_scores))
    return dict(
        accuracy=_accuracy(labels, scores),
        accuracy_interval=_bootstrap_interval(labels, scores, _accuracy, n_bootstrap, confidence, rng),
        auc=auc,
        auc_interval=_bootstrap_interval(labels, scores, roc_auc_score, n_bootstrap, confidence, rng),
        null_auc=null_auc,
        p_value=float(p_value),
        scores=scores,
        labels=labels,
        n_a=len(features_a),
        n_b=len(features_b),
        n_unique_a=n_unique[0],
        n_unique_b=n_unique[1],
        n_splits=n_splits,
    )


def distinguishability(
    images_a,
    images_b,
    features="shape",
    n_splits=5,
    n_bootstrap=1000,
    confidence=0.95,
    classifier=None,
    random_state=None,
    n_permutations=0,
    null_control=True,
    **feature_kwargs,
):
    """
    Accuracy and ROC AUC with which two sets of images are told apart.

    Each image is reduced to a feature vector (image_features; by default
    shape only, independent of brightness). A classifier (default:
    default_classifier) is evaluated by stratified k-fold cross-validation,
    so every image is scored by a model that did not see it. Accuracy
    (threshold 0.5) and ROC AUC are computed from these out-of-fold scores,
    and their intervals from bootstrap resamples of the realisations within
    each class.

    Parameters
    ----------
    images_a, images_b : sequence of numpy.ndarray
        Simulated images of structure A and of structure B, one per
        independent realisation (e.g. from run_replicates).
    features : str or callable, optional
        Feature extraction (see image_features). Default "shape".
    n_splits : int, optional
        Number of cross-validation folds. Default is 5 (reduced if a class
        has fewer images).
    n_bootstrap : int, optional
        Number of bootstrap resamples for the intervals. Default is 1000.
    confidence : float, optional
        Confidence level of the intervals. Default is 0.95.
    classifier : sklearn estimator, optional
        Classifier with predict_proba. Default: default_classifier.
    random_state : int, optional
        Seed for the folds, the bootstrap, the permutations and the null
        control.
    n_permutations : int, optional
        Number of label permutations for the p-value. Default 0 (no
        p-value).
    null_control : bool, optional
        Also score two random halves of structure A (null_auc, expected
        near 0.5). Default True; needs at least 4 images of A.
    **feature_kwargs
        Passed to the feature function.

    Returns
    -------
    dict
        accuracy, accuracy_interval, auc, auc_interval, null_auc (NaN if
        not computed), p_value (NaN if n_permutations is 0), scores
        (out-of-fold probability of class B per image), labels (0 for A,
        1 for B), n_a, n_b, n_unique_a, n_unique_b (distinct feature
        vectors per class), features (name or function), n_splits.
    """
    result = _score_features(
        image_features(images_a, features, **feature_kwargs),
        image_features(images_b, features, **feature_kwargs),
        n_splits=n_splits,
        n_bootstrap=n_bootstrap,
        confidence=confidence,
        classifier=classifier,
        random_state=random_state,
        n_permutations=n_permutations,
        null_control=null_control,
    )
    result["features"] = features
    return result


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
    n_permutations=0,
    null_control=True,
):
    """
    Distinguishability of two structures from precomputed features.

    Parameters
    ----------
    features_a, features_b : array-like
        Features measured on each realisation of structure A and of
        structure B: shape (n_realisations,) for one feature or
        (n_realisations, n_features).
    n_splits, n_bootstrap, confidence, classifier, random_state,
    n_permutations, null_control
        See distinguishability.

    Returns
    -------
    dict
        See distinguishability (without features).
    """
    return _score_features(
        features_a,
        features_b,
        n_splits=n_splits,
        n_bootstrap=n_bootstrap,
        confidence=confidence,
        classifier=classifier,
        random_state=random_state,
        n_permutations=n_permutations,
        null_control=null_control,
    )


def _feature_columns(table, features):
    # column names of a feature set in a feature table: a column name, a
    # named feature set (columns "<name>" or "<name>_..."), or a list
    if not isinstance(features, str):
        return list(features)
    if features in table.columns:
        return [features]
    columns = [c for c in table.columns if c.startswith(features + "_")]
    if not columns:
        raise ValueError(f"No columns for features {features!r} in the table")
    return columns


def _score_condition(condition, a, b, controls, options):
    # scores of one condition; a, b: feature matrices of the two classes
    score = _score_features(a, b, **options)
    row = dict(condition)
    row.update(
        n_a=score["n_a"], n_b=score["n_b"],
        n_unique_a=score["n_unique_a"], n_unique_b=score["n_unique_b"],
        accuracy=score["accuracy"],
        accuracy_low=score["accuracy_interval"][0],
        accuracy_high=score["accuracy_interval"][1],
        auc=score["auc"], auc_low=score["auc_interval"][0], auc_high=score["auc_interval"][1],
        null_auc=score["null_auc"], p_value=score["p_value"],
    )
    control_options = dict(options, n_permutations=0, null_control=False)
    for name, (control_a, control_b) in controls.items():
        control = _score_features(control_a, control_b, **control_options)
        row.update({
            f"{name}_auc": control["auc"],
            f"{name}_auc_low": control["auc_interval"][0],
            f"{name}_auc_high": control["auc_interval"][1],
        })
    return row, dict(scores=score["scores"], labels=score["labels"])


def distinguishability_from_features(
    table,
    label_column,
    classes,
    features,
    by=None,
    controls=None,
    replicate_column=None,
    n_splits=5,
    n_bootstrap=1000,
    confidence=0.95,
    classifier=None,
    random_state=None,
    n_permutations=0,
    null_control=True,
    n_workers=1,
    progress=True,
    return_scores=False,
):
    """
    Distinguishability of two classes in a feature table, per condition.

    The table has one row per image (e.g. the analysis dataframe of
    run_parameter_sweep with feature_metrics, for both structures, with a
    column naming the structure). For every condition (combination of the
    by columns), the images of the two classes are scored as in
    distinguishability: out-of-fold accuracy and ROC AUC with bootstrap
    intervals, the null control and, optionally, a permutation p-value.
    Nothing is simulated, so other features or classifiers can be tried on
    the same table in seconds.

    Parameters
    ----------
    table : pandas.DataFrame
        Feature table, one row per image.
    label_column : str
        Column naming the class of each image (e.g. "model").
    classes : tuple of 2
        The two class values to compare (A, B), e.g. ("dome", "flat").
    features : str or list of str
        Feature columns: a named feature set ("whitened" selects the
        columns whitened_...), a column name, or a list of columns.
    by : str or list of str, optional
        Columns defining a condition (e.g. ["modality_name",
        "labelling_efficiency"]). Default: one condition, the whole table.
    controls : dict, optional
        {name: features}: also score these features alone, as controls
        (e.g. {"brightness": "brightness"}), reported as <name>_auc with
        its interval.
    replicate_column : str, optional
        Column used to order the images of each class (e.g. "Replica"),
        for reproducible folds.
    n_splits, n_bootstrap, confidence, classifier, random_state,
    n_permutations, null_control
        See distinguishability.
    n_workers : int, optional
        Conditions scored in parallel processes. Default 1. In a script,
        call it under ``if __name__ == "__main__":`` when n_workers > 1.
    progress : bool, optional
        Show a progress bar. Default True.
    return_scores : bool, optional
        Also return the out-of-fold scores and labels of every condition.

    Returns
    -------
    pandas.DataFrame
        One row per condition: the by columns, features, n_a, n_b,
        n_unique_a, n_unique_b, accuracy, accuracy_low, accuracy_high, auc,
        auc_low, auc_high, null_auc, p_value and the control columns.
    dict
        Only if return_scores: {condition tuple: {"scores", "labels"}}.
    """
    columns = _feature_columns(table, features)
    control_columns = {name: _feature_columns(table, f) for name, f in (controls or {}).items()}
    by = [] if by is None else ([by] if isinstance(by, str) else list(by))
    options = dict(n_splits=n_splits, n_bootstrap=n_bootstrap, confidence=confidence,
                   classifier=classifier, random_state=random_state,
                   n_permutations=n_permutations, null_control=null_control)
    groups = table.groupby(by, sort=False) if by else [((), table)]
    jobs = []
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        a = group[group[label_column] == classes[0]]
        b = group[group[label_column] == classes[1]]
        if len(a) < 2 or len(b) < 2:
            warnings.warn(f"Skipping condition {dict(zip(by, key))}: fewer than 2 images "
                          f"of a class", RuntimeWarning, stacklevel=2)
            continue
        if replicate_column is not None:
            a, b = a.sort_values(replicate_column), b.sort_values(replicate_column)
        condition = dict(zip(by, key))
        condition["features"] = features if isinstance(features, str) else "+".join(columns)
        jobs.append((condition, a[columns].to_numpy(float), b[columns].to_numpy(float),
                     {name: (a[c].to_numpy(float), b[c].to_numpy(float))
                      for name, c in control_columns.items()}))
    outputs = [None] * len(jobs)
    bar = tqdm(total=len(jobs), desc="conditions", disable=not progress)
    if n_workers <= 1:
        for k, job in enumerate(jobs):
            outputs[k] = _score_condition(*job, options)
            bar.update()
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(min(n_workers, max(len(jobs), 1))) as pool:
            futures = {pool.submit(_score_condition, *job, options): k
                       for k, job in enumerate(jobs)}
            for future in as_completed(futures):
                outputs[futures[future]] = future.result()
                bar.update()
    bar.close()
    result = pd.DataFrame([row for row, _ in outputs])
    if return_scores:
        scores = {tuple(row[c] for c in by): s for row, s in outputs}
        return result, scores
    return result


def feature_ablation(table, label_column, classes, groups, by=None, **kwargs):
    """
    AUC of each group of features alone, per condition.

    Shows which features carry the distinguishability of two classes, e.g.
    whether a score comes from the size of the particles, their outline,
    their internal structure or their brightness.

    Parameters
    ----------
    table, label_column, classes, by
        See distinguishability_from_features.
    groups : dict
        {group name: features} (named feature set, column or list of
        columns), e.g. {"profile": [...], "moments": ["whitened_m4",
        "whitened_m6"], "brightness": "brightness"}.
    **kwargs
        Passed to distinguishability_from_features (by default without the
        null control and permutations).

    Returns
    -------
    pandas.DataFrame
        The rows of distinguishability_from_features for every group, with
        a feature_group column.
    """
    kwargs.setdefault("null_control", False)
    kwargs.setdefault("n_permutations", 0)
    tables = []
    for name, features in groups.items():
        result = distinguishability_from_features(table, label_column, classes, features,
                                                  by=by, **kwargs)
        result.insert(0, "feature_group", name)
        tables.append(result)
    return pd.concat(tables, ignore_index=True)


def plot_distinguishability_heatmaps(
    scores,
    x,
    y,
    facet=None,
    metric="auc",
    facet_order=None,
    highlight=None,
    vmin=0.5,
    vmax=1.0,
    cmap="viridis",
    annotate=True,
    axes=None,
):
    """
    Heatmaps of a score against two parameters, one panel per facet.

    Parameters
    ----------
    scores : pandas.DataFrame
        Output of distinguishability_from_features (or any table with the
        x, y, facet and metric columns), one row per condition.
    x, y : str
        Columns on the horizontal and vertical axes.
    facet : str, optional
        Column with one panel per value (e.g. "modality_name").
    metric : str, optional
        Column shown. Default "auc".
    facet_order : list, optional
        Order of the panels. Default: order of appearance.
    highlight : list of dict, optional
        Cells to outline, each {"x": value, "y": value, "facet": value
        (if facet), "colour": colour (default "red")}.
    vmin, vmax : float, optional
        Colour scale. Default 0.5 (chance) to 1.
    cmap : str, optional
        Colour map. Default "viridis".
    annotate : bool, optional
        Print the value in every cell. Default True.
    axes : sequence of matplotlib.axes.Axes, optional
        One axes per panel. Default: a new figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    facets = [None] if facet is None else (facet_order or list(dict.fromkeys(scores[facet])))
    if axes is None:
        figure, axes = plt.subplots(1, len(facets), figsize=(3.2 * len(facets), 3.0),
                                    sharey=True, constrained_layout=True, squeeze=False)
        axes = axes[0]
    figure = axes[0].figure
    xs = sorted(scores[x].unique())
    ys = sorted(scores[y].unique())
    image = None
    for ax, value in zip(axes, facets):
        data = scores if facet is None else scores[scores[facet] == value]
        grid = (data.pivot_table(index=y, columns=x, values=metric)
                .reindex(index=ys, columns=xs).to_numpy())
        image = ax.imshow(grid, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        if annotate:
            for i in range(len(ys)):
                for j in range(len(xs)):
                    if np.isfinite(grid[i, j]):
                        level = (grid[i, j] - vmin) / (vmax - vmin)
                        ax.text(j, i, f"{grid[i, j]:.2f}", ha="center", va="center",
                                fontsize=7, color="white" if level < 0.6 else "black")
        for cell in highlight or []:
            if facet is not None and cell.get("facet") != value:
                continue
            j, i = xs.index(cell["x"]), ys.index(cell["y"])
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, linewidth=2.5,
                                   edgecolor=cell.get("colour", "red")))
        ax.set_xticks(range(len(xs)), [f"{v:g}" for v in xs], rotation=45)
        ax.set_yticks(range(len(ys)), [f"{v:g}" for v in ys])
        ax.set_xlabel(x)
        if value is not None:
            ax.set_title(str(value))
    axes[0].set_ylabel(y)
    figure.colorbar(image, ax=list(axes), label=metric, shrink=0.9)
    return figure


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
    "null_auc",
    "p_value",
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
                    null_auc=score["null_auc"],
                    p_value=score["p_value"],
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