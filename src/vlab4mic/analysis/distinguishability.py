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
"""

import numpy as np
from scipy import ndimage
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


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


FEATURES = {"radial_profile": radial_profile_features, "pixels": pixel_features}


def image_features(images, features="radial_profile", **kwargs):
    """
    Feature matrix of a list of images.

    Parameters
    ----------
    images : sequence of numpy.ndarray
        Images (2-D, or stacks of frames that are summed).
    features : str or callable, optional
        "radial_profile" (default, rotation invariant), "pixels", or a
        function image -> 1-D feature vector.
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
