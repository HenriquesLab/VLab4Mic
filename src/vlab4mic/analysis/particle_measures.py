"""Measurements on single simulated particles.

Used for the comparisons with published measurements: ring radius and
width (nuclear pore), corner detection and apparent breaks (gaps), and the
number of resolved labelled sites (PCNA).
"""

import numpy as np
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation


def fit_circle(points):
    """
    Least-squares circle through 2-D points (algebraic Kasa fit).

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 (or Nx3, z ignored) coordinates, N >= 3.

    Returns
    -------
    centre : numpy.ndarray
        Circle centre (x, y).
    radius : float
        Circle radius.
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    if xy.shape[0] < 3:
        raise ValueError("A circle fit needs at least 3 points")
    A = np.c_[2 * xy, np.ones(len(xy))]
    b = (xy**2).sum(axis=1)
    cx, cy, c = np.linalg.lstsq(A, b, rcond=None)[0]
    radius = np.sqrt(c + cx**2 + cy**2)
    return np.array([cx, cy]), float(radius)


def ring_measures(points):
    """
    Radius and width of a ring of localisations or emitters.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 or Nx3 coordinates of one ring seen from above.

    Returns
    -------
    dict
        radius (circle fit), width (standard deviation of the distances to
        the centre), centre.
    """
    centre, radius = fit_circle(points)
    distances = np.linalg.norm(np.asarray(points, float)[:, :2] - centre, axis=1)
    return dict(radius=radius, width=float(distances.std()), centre=centre)


def angular_gaps(points, centre=None):
    """
    Gaps in azimuth between consecutive points around a ring.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 or Nx3 coordinates.
    centre : array-like, optional
        Ring centre. Default: circle fit.

    Returns
    -------
    numpy.ndarray
        Sorted angular gaps in degrees (they sum to 360).
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    if centre is None:
        centre, _ = fit_circle(xy)
    angles = np.sort(np.degrees(np.arctan2(*(xy - centre).T[::-1])) % 360)
    gaps = np.diff(np.r_[angles, angles[0] + 360])
    return np.sort(gaps)


def has_apparent_break(points, max_gap_deg, centre=None):
    """
    Whether a ring shows a gap wider than max_gap_deg.

    For a ring of n evenly spaced subunits (e.g. the eight corners of the
    nuclear pore, 45 degrees apart), a gap wider than the subunit spacing
    means at least one subunit carries no label and appears as a break.

    Parameters
    ----------
    points : numpy.ndarray
        Coordinates of the ring.
    max_gap_deg : float
        Largest gap of an intact ring, in degrees.
    centre : array-like, optional
        Ring centre. Default: circle fit.

    Returns
    -------
    bool
    """
    if len(points) < 3:
        return True
    return bool(angular_gaps(points, centre).max() > max_gap_deg)


def count_occupied_sectors(points, n_sectors, centre=None, offset_deg=None):
    """
    Number of angular sectors of a ring that contain at least one point.

    Used to count detected corners of an n-fold symmetric ring.

    Parameters
    ----------
    points : numpy.ndarray
        Coordinates of the ring.
    n_sectors : int
        Number of sectors (e.g. 8 for the nuclear pore).
    centre : array-like, optional
        Ring centre. Default: circle fit.
    offset_deg : float, optional
        Angle of the first sector boundary. Default: chosen so that the
        sectors are centred on the densest directions.

    Returns
    -------
    int
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    if centre is None:
        centre, _ = fit_circle(xy)
    angles = np.degrees(np.arctan2(*(xy - centre).T[::-1])) % 360
    width = 360 / n_sectors
    if offset_deg is None:
        # centre sectors on the circular mean of angles modulo the spacing
        phase = np.angle(np.exp(1j * np.radians(angles * n_sectors)).mean())
        offset_deg = np.degrees(phase) / n_sectors - width / 2
    sectors = np.floor(((angles - offset_deg) % 360) / width).astype(int)
    return int(np.unique(sectors).size)


def count_resolved_sites(image, pixelsize_nm, min_separation_nm, threshold_rel=0.2):
    """
    Number of separate peaks in a rendered single-particle image.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image (or stack of frames, which is summed).
    pixelsize_nm : float
        Pixel size in nm.
    min_separation_nm : float
        Minimum distance between peaks, in nm.
    threshold_rel : float, optional
        Minimum peak height relative to the maximum. Default 0.2.

    Returns
    -------
    int
    """
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        image = image.sum(axis=0)
    if image.max() <= 0:
        return 0
    min_distance = max(1, int(round(min_separation_nm / pixelsize_nm)))
    peaks = peak_local_max(
        image, min_distance=min_distance, threshold_rel=threshold_rel,
        exclude_border=False,
    )
    return int(len(peaks))


def sites_resolved(points, n_sites, random_state=0):
    """
    Whether n_sites labelled sites are resolved from their localisations.

    Gaussian mixtures with 1 to n_sites spherical components are fitted to
    the localisations (x, y), and the sites are resolved if the Bayesian
    information criterion (BIC) selects n_sites components, i.e. the data
    support n_sites separate spots rather than fewer, overlapping ones.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 or Nx3 localisation coordinates.
    n_sites : int
        Number of sites.
    random_state : int, optional
        Seed for the mixture fits.

    Returns
    -------
    resolved : bool
    n_components : int
        Number of components selected by BIC.
    """
    from sklearn.mixture import GaussianMixture

    xy = np.asarray(points, dtype=float)[:, :2]
    if xy.shape[0] < 2 * n_sites:
        return False, 0
    bic = [
        GaussianMixture(k, covariance_type="spherical", n_init=2, random_state=random_state)
        .fit(xy)
        .bic(xy)
        for k in range(1, n_sites + 1)
    ]
    n_components = int(np.argmin(bic)) + 1
    return n_components == n_sites, n_components


# Global methods return a single threshold value
_GLOBAL_METHODS = {
    "otsu": filters.threshold_otsu,
    "li": filters.threshold_li,
    "yen": filters.threshold_yen,
    "triangle": filters.threshold_triangle,
    "isodata": filters.threshold_isodata,
    "mean": filters.threshold_mean,
    "minimum": filters.threshold_minimum,
}

# Local (adaptive) methods return a threshold image with the same shape as the input
_LOCAL_METHODS = {
    "local": filters.threshold_local,       # e.g. block_size=35, offset=0
    "niblack": filters.threshold_niblack,   # e.g. window_size=15, k=0.2
    "sauvola": filters.threshold_sauvola,   # e.g. window_size=15, k=0.2
}

AVAILABLE_METHODS = sorted([*_GLOBAL_METHODS, *_LOCAL_METHODS, "multiotsu"])


def compute_threshold(image, method="otsu", method_kwargs=None, multiotsu_index=-1):
    """
    Compute a threshold for a 2D image.

    Parameters
    ----------
    image : 2D array
    method : str or callable
        Name of a skimage threshold method (see AVAILABLE_METHODS) or a callable
        ``f(image, **method_kwargs)`` returning a scalar or an array of image shape.
    method_kwargs : dict, optional
        Extra keyword arguments passed to the threshold function.
    multiotsu_index : int
        For ``method="multiotsu"``: which of the returned thresholds to use
        (-1 = highest, i.e. keep only the brightest class).

    Returns
    -------
    threshold : float or np.ndarray
    """
    kwargs = dict(method_kwargs or {})

    if callable(method):
        return method(image, **kwargs)

    name = method.lower()
    if name in _GLOBAL_METHODS:
        return _GLOBAL_METHODS[name](image, **kwargs)
    if name in _LOCAL_METHODS:
        if name == "local":
            kwargs.setdefault("block_size", 35)  # must be odd
        return _LOCAL_METHODS[name](image, **kwargs)
    if name == "multiotsu":
        thresholds = filters.threshold_multiotsu(image, **kwargs)
        return thresholds[multiotsu_index]

    raise ValueError(
        f"Unknown threshold method '{method}'. Choose from {AVAILABLE_METHODS} "
        f"or pass a callable."
    )


def threshold_mask_and_measure(
    image,
    *,
    method="otsu",           # threshold method name or callable
    method_kwargs=None,      # dict of kwargs for the threshold method
    multiotsu_index=-1,      # which multi-Otsu threshold to use
    threshold=None,          # explicit threshold; overrides `method`
    threshold_scale=1.0,     # multiply the computed threshold
    smooth_sigma=0.0,        # Gaussian pre-smoothing (0 = off)
    invert=False,            # True if objects are darker than background
    min_size=0,              # remove objects smaller than this (pixels)
    fill_holes=False,        # fill holes inside objects
    opening_radius=0,        # morphological opening radius (0 = off)
    clear_border=False,      # drop objects touching the image border
    connectivity=2,          # 1 = 4-connected, 2 = 8-connected
    pixel_size=1.0,          # physical size per pixel (e.g. nm/px)
    extra_properties=None,   # list of additional regionprops names
):
    """
    Create a binary mask from a 2D grey-value image using a selectable
    thresholding method, then label connected objects and measure each one.

    Returns
    -------
    mask : np.ndarray (bool)
    table : pd.DataFrame
        One row per object. Metadata in ``table.attrs``:
        ``method``, ``threshold`` (float, or None for local methods),
        ``n_objects``, ``pixel_size``.
    """
    img = np.asarray(image, dtype=float)
    if img.ndim != 2:
        raise ValueError(f"Expected a 2D array, got shape {img.shape}")

    # --- preprocessing ---
    work = filters.gaussian(img, sigma=smooth_sigma, preserve_range=True) \
        if smooth_sigma > 0 else img

    # --- threshold ---
    if threshold is None:
        threshold = compute_threshold(
            work, method=method, method_kwargs=method_kwargs,
            multiotsu_index=multiotsu_index,
        )
        method_used = method if isinstance(method, str) else getattr(method, "__name__", "custom")
    else:
        method_used = "manual"
    threshold = np.asarray(threshold, dtype=float) * threshold_scale

    mask = work < threshold if invert else work > threshold

    # --- mask cleanup ---
    if opening_radius > 0:
        mask = morphology.binary_opening(mask, morphology.disk(opening_radius))
    if fill_holes:
        mask = ndi.binary_fill_holes(mask)
    if min_size > 0:
        mask = morphology.remove_small_objects(
            mask, min_size=min_size, connectivity=connectivity
        )
    if clear_border:
        mask = segmentation.clear_border(mask)

    # --- labelling ---
    labels = measure.label(mask, connectivity=connectivity)

    # --- measurements ---
    properties = [
        "label",
        "area",
        "perimeter",
        "centroid",
        "bbox",
        "equivalent_diameter_area",
        "major_axis_length",
        "minor_axis_length",
        "orientation",
        "eccentricity",
        "solidity",
        "extent",
        "mean_intensity",
        "max_intensity",
        "min_intensity",
    ]
    if extra_properties:
        properties += [p for p in extra_properties if p not in properties]

    table = pd.DataFrame(
        measure.regionprops_table(labels, intensity_image=img, properties=properties)
    )

    if not table.empty:
        table["circularity"] = (
            4 * np.pi * table["area"] / np.square(table["perimeter"]).replace(0, np.nan)
        )
        table["integrated_intensity"] = table["mean_intensity"] * table["area"]

        table["area_physical"] = table["area"] * pixel_size**2
        table["perimeter_physical"] = table["perimeter"] * pixel_size
        table["equivalent_diameter_physical"] = (
            table["equivalent_diameter_area"] * pixel_size
        )
        table["major_axis_physical"] = table["major_axis_length"] * pixel_size
        table["minor_axis_physical"] = table["minor_axis_length"] * pixel_size

        table = table.rename(columns={"centroid-0": "centroid_y",
                                      "centroid-1": "centroid_x"})

    table.attrs["method"] = method_used
    table.attrs["threshold"] = float(threshold) if threshold.ndim == 0 else None
    table.attrs["n_objects"] = int(labels.max())
    table.attrs["pixel_size"] = pixel_size

    return mask, table, work