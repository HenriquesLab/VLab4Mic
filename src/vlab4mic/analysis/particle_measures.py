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


def measure_filament_width(
    image,
    pixel_size=1.0,
    search_halfwidth=None,
    row_binning=1,
    background_width=5,
    min_snr=2.0,
    mode="single",
    min_sum = 0.001,
):
    """
    Measure the apparent width (FWHM) of a vertical filament centred in a 2D image.

    Each row (or bin of rows) is treated as a transversal line profile. The FWHM
    is measured per profile with sub-pixel linear interpolation. The FWHM of the
    centre-aligned average profile along the whole filament is also reported.

    Parameters
    ----------
    image : 2D array
        Grayscale image with a bright, roughly vertical filament near the centre column.
    pixel_size : float
        Physical size of a pixel. Widths are returned in these units.
    search_halfwidth : int or None
        Half-width (in px) of the window around the centre column where the peak
        is searched. Defaults to a quarter of the image width.
    row_binning : int
        Number of consecutive rows averaged into one profile (improves SNR).
    background_width : int
        Number of pixels at each edge of the search window used to estimate
        the local background.
    min_snr : float
        Profiles whose (peak - background) / background_std is below this value
        are discarded (NaN).

    Returns
    -------
    dict with keys
        "widths"          : 1D array, FWHM per profile (NaN where invalid)
        "centers"         : 1D array, sub-pixel filament centre per profile (px)
        "rows"            : 1D array, central row index of each profile
        "mean", "median", "std" : statistics of valid per-profile widths
        "n_valid"         : number of valid profiles
        "fwhm_avg_profile": FWHM of the centre-aligned mean profile
        "avg_profile"     : (x, profile) tuple of the aligned mean profile
    """
    img = np.asarray(image, dtype=float)
    if img.ndim != 2:
        raise ValueError("image must be 2D")
    n_rows, n_cols = img.shape

    if image.sum() < min_sum:
        return  {
            "widths": None,
            "centers": None,
            "rows": None,
            "mean": None,
            "median": None,
            "std": None,
            "n_valid": None,
            "fwhm_avg_profile": None,
            "avg_profile": None
        }
    if search_halfwidth is None:
        search_halfwidth = max(n_cols // 4, 3)
    c0 = n_cols // 2
    x_lo = max(c0 - search_halfwidth, 0)
    x_hi = min(c0 + search_halfwidth + 1, n_cols)
    x = np.arange(x_lo, x_hi, dtype=float)

    # Bin rows
    n_bins = n_rows // row_binning
    binned = img[: n_bins * row_binning, x_lo:x_hi]
    binned = binned.reshape(n_bins, row_binning, -1).mean(axis=1)
    rows = np.arange(n_bins) * row_binning + (row_binning - 1) / 2.0

    widths = np.full(n_bins, np.nan)
    centers = np.full(n_bins, np.nan)
    aligned = []

    for i, prof in enumerate(binned):
        res = _fwhm_1d(prof, x, background_width, min_snr, mode)
        if res is None:
            continue
        widths[i], centers[i], prof_bg = res
        aligned.append((x - centers[i], prof_bg))

    valid = widths[~np.isnan(widths)] * pixel_size
    widths *= pixel_size

    # Centre-aligned mean profile over the whole filament
    fwhm_avg = np.nan
    avg_profile = (None, None)
    if aligned:
        step = 0.25  # sub-pixel sampling grid
        x_al = np.arange(-search_halfwidth, search_halfwidth + step, step)
        try:
            stack = np.array(
                [np.interp(x_al, xs, p, left=np.nan, right=np.nan) for xs, p in aligned]
            )
            mean_prof = np.nanmean(stack, axis=0)
            ok = ~np.isnan(mean_prof)
            x_al, mean_prof = x_al[ok], mean_prof[ok]
            # Profile is already background-subtracted -> baseline 0
            res = _fwhm_1d(mean_prof, x_al, background_width=0, min_snr=0, mode=mode)
            if res is not None:
                fwhm_avg = res[0] * pixel_size
            avg_profile = (x_al * pixel_size, mean_prof)
        except:
            avg_profile = (None, None)


    return {
        "widths": widths,
        "centers": centers,
        "rows": rows,
        "mean": np.mean(valid) if valid.size else np.nan,
        "median": np.median(valid) if valid.size else np.nan,
        "std": np.std(valid) if valid.size else np.nan,
        "n_valid": int(valid.size),
        "fwhm_avg_profile": fwhm_avg,
        "avg_profile": avg_profile,
    }


def _fwhm_1d(profile, x, background_width, min_snr, mode="single"):
    """
    FWHM of a 1D profile using linear interpolation at half maximum.

    mode="single": width of the brightest peak (walk outward from the maximum
                   until the profile first drops below half max).
    mode="total" : total width of the structure, from the outermost left
                   half-max crossing to the outermost right one. Several
                   peaks are treated as one structure.

    Returns (width, center, background_subtracted_profile) or None if invalid.
    """
    if mode not in ("single", "total"):
        raise ValueError("mode must be 'single' or 'total'")

    p = np.asarray(profile, dtype=float)
    n = p.size

    if background_width > 0 and 2 * background_width < n:
        edges = np.concatenate([p[:background_width], p[-background_width:]])
        bg, bg_std = np.median(edges), np.std(edges)
    else:
        bg, bg_std = 0.0, 0.0

    p = p - bg
    i_max = int(np.argmax(p))
    amp = p[i_max]
    if amp <= 0:
        return None, None, None
    if bg_std > 0 and amp / bg_std < min_snr:
        return None, None, None

    half = amp / 2.0

    if mode == "single":
        # Walk outward from the peak until the profile drops below half max
        i = i_max
        while i > 0 and p[i] >= half:
            i -= 1
        if p[i] >= half:
            return None
        j = i_max
        while j < n - 1 and p[j] >= half:
            j += 1
        if p[j] >= half:
            return None
        i_in, i_out = i + 1, i      # left crossing between i and i+1
        j_in, j_out = j - 1, j      # right crossing between j-1 and j
    else:
        # Outermost samples at or above half max
        above = np.flatnonzero(p >= half)
        i_first, i_last = above[0], above[-1]
        if i_first == 0 or i_last == n - 1:
            return None  # structure touches the window edge
        i_in, i_out = i_first, i_first - 1
        j_in, j_out = i_last, i_last + 1

    x_left = x[i_out] + (half - p[i_out]) * (x[i_in] - x[i_out]) / (p[i_in] - p[i_out])
    x_right = x[j_in] + (half - p[j_in]) * (x[j_out] - x[j_in]) / (p[j_out] - p[j_in])

    width = x_right - x_left
    center = 0.5 * (x_left + x_right)
    return width, center, p