"""Measurements on single simulated particles.

Used for the comparisons with published measurements: ring radius and
width (nuclear pore), corner detection and apparent breaks (gaps), and the
number of resolved labelled sites (PCNA).
"""

import numpy as np
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from scipy import optimize
from skimage import filters, morphology, measure, segmentation
import pandas as pd


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


def fit_circle_geometric(points):
    """
    Least-squares circle through 2-D points (geometric fit).

    Minimises the spread of the distances to the centre, starting from the
    mean of the points. Less biased than fit_circle for noisy or partial rings.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 (or Nx3, z ignored) coordinates, N >= 3.

    Returns
    -------
    centre : numpy.ndarray
        Circle centre (x, y).
    radius : float
        Mean distance of the points to the centre.
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    if xy.shape[0] < 3:
        raise ValueError("A circle fit needs at least 3 points")

    def residuals(centre):
        distances = np.linalg.norm(xy - centre, axis=1)
        return distances - distances.mean()

    centre, _ = optimize.leastsq(residuals, x0=xy.mean(axis=0))
    radius = np.linalg.norm(xy - centre, axis=1).mean()
    return centre, float(radius)


def fit_circle_fixed_radius(points, radius):
    """
    Least-squares centre of a circle of known radius through 2-D points.

    The starting point is the centre of fit_circle, so the fit does not
    drift away from the ring when its radius differs from `radius`.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 (or Nx3, z ignored) coordinates.
    radius : float
        Circle radius.

    Returns
    -------
    numpy.ndarray
        Circle centre (x, y).
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    start = fit_circle(xy)[0] if len(xy) >= 3 else xy.mean(axis=0)
    fit = optimize.least_squares(lambda c: np.linalg.norm(xy - c, axis=1) - radius, start)
    return fit.x


def split_by_midpoint_z(points, split=True):
    """
    Points above the mean z, e.g. one of the two rings of a nuclear pore.

    Parameters
    ----------
    points : numpy.ndarray or None
        Nx3 coordinates.
    split : bool, optional
        If False, points are returned unchanged. Default True.

    Returns
    -------
    numpy.ndarray or None
    """
    if points is None or not split:
        return points
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Splitting by z needs Nx3 coordinates")
    return points[points[:, 2] > points[:, 2].mean()]


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


def fit_spots(points, max_spots, random_state=0):
    """
    Separate spots supported by a set of localisations.

    Gaussian mixtures with 1 to max_spots spherical components are fitted to
    the localisations (x, y), and the Bayesian information criterion (BIC)
    selects the number of components. Sites closer than the localisation
    precision allows to separate are fitted as a single spot.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 or Nx3 localisation coordinates, N >= 1.
    max_spots : int
        Largest number of spots tried (limited to N // 2).
    random_state : int, optional
        Seed for the mixture fits.

    Returns
    -------
    means : numpy.ndarray
        kx2 spot centres (x, y).
    counts : numpy.ndarray
        Localisations assigned to each spot.
    """
    from sklearn.mixture import GaussianMixture

    xy = np.asarray(points, dtype=float)[:, :2]
    k_max = max(1, min(max_spots, len(xy) // 2))
    fits = [
        GaussianMixture(k, covariance_type="spherical", n_init=2, random_state=random_state).fit(xy)
        for k in range(1, k_max + 1)
    ]
    best = min(fits, key=lambda fit: fit.bic(xy))
    counts = np.bincount(best.predict(xy), minlength=best.n_components)
    return best.means_, counts


def count_ring_corners(
    points, radius, r_min, r_max, n_corners, n_steps=1000, min_per_section=1,
    resolve_spots=False, random_state=0,
):
    """
    Number of detected corners of a ring-shaped particle.

    1. A circle of fixed radius is fitted and coordinates converted to (phi_i, r_i).
    2. Localisations with r_i < r_min or r_i > r_max are discarded as background.
    3. phi_rot minimises the summed squared angular distance between phi_rot and
       phi_i modulo 2*pi/n_corners (wrapped, so both sides of an edge count as close).
    4. Corners = segments of width 2*pi/n_corners, starting at phi_rot - pi/n_corners,
       that contain at least min_per_section localisations.
    5. If resolve_spots, the kept localisations must also form separate spots
       (fit_spots, at most n_corners): a corner is an occupied segment that
       contains the centre of a spot with at least min_per_section
       localisations and inside [r_min, r_max]. A single cloud spread over
       several segments then counts as at most one corner.

    Unlike count_occupied_sectors, background is removed with a radial window
    and a minimum number of localisations per segment can be required.

    Parameters
    ----------
    points : numpy.ndarray or None
        Nx2 or Nx3 localisations of one particle, in nm.
    radius : float
        Fixed circle radius in nm (mean radius of the data set).
    r_min, r_max : float
        Radial window kept around the centre, in nm.
    n_corners : int
        Rotational symmetry (8 for the nuclear pore, 3 for PCNA).
    n_steps : int, optional
        Candidate rotations tried within one segment. Default 1000.
    min_per_section : int, optional
        Localisations a segment needs to count as a corner. Default 1.
    resolve_spots : bool, optional
        Also require each corner to be a separate spot (step 5). Default False.
    random_state : int, optional
        Seed for the spot fits.

    Returns
    -------
    dict with
        n_corners      : int, number of detected corners
        n_segments_occupied : int, occupied segments (step 4 only)
        spot_means     : kx2 spot centres, original frame (empty if not resolve_spots)
        spot_counts    : localisations per spot
        spot_is_corner : bool mask, spots that count as a corner
        points         : input coordinates
        centre         : fitted centre (x, y), None if no points
        radius, r_min, r_max, symmetry : parameters used
        r, phi         : polar coordinates of all localisations
        kept           : bool mask, localisations inside [r_min, r_max]
        phi_rot        : ring rotation in radians (None if nothing kept)
        aligned        : Nx2 coords of all localisations, centred and rotated by
                         -phi_rot (segment k is centred on angle k*2*pi/n_corners)
        segment        : segment index of each kept localisation
        segment_counts : localisations per segment
        segment_edges  : n_corners+1 segment edges, radians, original frame
    """
    if points is None:
        points = np.empty((0, 2))
    points = np.asarray(points, dtype=float).reshape(-1, np.shape(points)[-1] if np.ndim(points) == 2 else 3)
    xy = points[:, :2]
    period = 2 * np.pi / n_corners
    result = dict(
        n_corners=0, points=points, centre=None, radius=radius, r_min=r_min,
        r_max=r_max, symmetry=n_corners, r=np.empty(0), phi=np.empty(0),
        kept=np.zeros(len(xy), bool), phi_rot=None, aligned=np.empty((0, 2)),
        segment=np.empty(0, int), segment_counts=np.zeros(n_corners, int),
        segment_edges=None, n_segments_occupied=0, spot_means=np.empty((0, 2)),
        spot_counts=np.empty(0, int), spot_is_corner=np.empty(0, bool),
    )
    if len(xy) == 0:
        return result

    # 1. fixed-radius circle fit, polar coordinates
    centre = fit_circle_fixed_radius(xy, radius)
    dx, dy = (xy - centre).T
    r = np.hypot(dx, dy)
    phi = np.arctan2(dy, dx)
    # 2. radial filter
    kept = (r >= r_min) & (r <= r_max)
    result.update(centre=centre, r=r, phi=phi, kept=kept)
    if not kept.any():
        return result

    # 3. rotation of the ring
    phi_k = phi[kept]
    candidates = np.linspace(0, period, n_steps, endpoint=False)
    diff = (phi_k[None, :] - candidates[:, None] + period / 2) % period - period / 2
    phi_rot = candidates[np.argmin((diff**2).sum(axis=1))]

    # 4. occupied segments
    start = phi_rot - period / 2
    segment = np.floor(((phi_k - start) % (2 * np.pi)) / period).astype(int) % n_corners
    counts = np.bincount(segment, minlength=n_corners)

    c, s = np.cos(-phi_rot), np.sin(-phi_rot)
    aligned = np.c_[dx, dy] @ np.array([[c, s], [-s, c]])
    n_occupied = int((counts >= min_per_section).sum())
    result.update(
        n_corners=n_occupied, n_segments_occupied=n_occupied, phi_rot=float(phi_rot),
        aligned=aligned, segment=segment, segment_counts=counts,
        segment_edges=start + period * np.arange(n_corners + 1),
    )

    # 5. corners must be separate spots
    if resolve_spots:
        spot_means, spot_counts = fit_spots(xy[kept], n_corners, random_state=random_state)
        sx, sy = (spot_means - centre).T
        spot_r = np.hypot(sx, sy)
        spot_segment = np.floor(((np.arctan2(sy, sx) - start) % (2 * np.pi)) / period).astype(int) % n_corners
        is_corner = (
            (spot_r >= r_min) & (spot_r <= r_max)
            & (spot_counts >= min_per_section)
            & (counts[spot_segment] >= min_per_section)
        )
        result.update(
            n_corners=int(np.unique(spot_segment[is_corner]).size),
            spot_means=spot_means, spot_counts=spot_counts, spot_is_corner=is_corner,
        )
    return result


def ring_corners_from_replicates(
    replicates,
    radius,
    r_min,
    r_max,
    n_corners,
    modality="SMLM",
    channel="ch0",
    split_z=False,
    min_per_section=1,
    min_locs_to_circlefit=10,
    max_radius=None,
    resolve_spots=False,
):
    """
    Corners (count_ring_corners) and radius (fit_circle_geometric) of the
    single particle simulated in each replicate.

    Parameters
    ----------
    replicates : list of dict
        Output of run_replicates.
    radius, r_min, r_max, n_corners, min_per_section, resolve_spots
        Passed to count_ring_corners.
    modality : str, optional
        Modality name. Default "SMLM".
    channel : str, optional
        Channel name. Default "ch0".
    split_z : bool, optional
        Keep only localisations above the mean z (split_by_midpoint_z).
    min_locs_to_circlefit : int, optional
        Localisations needed to fit the radius (at least 3). Default 10.
    max_radius : float, optional
        Fitted radii at or above this value are set to None.

    Returns
    -------
    table : pandas.DataFrame
        One row per replicate: replicate, n_localisations, radius_fit and
        the count_ring_corners output.
    images : list of numpy.ndarray
        Image of each replicate.
    rows : list of dict
        Rows of the table.
    """
    rows = []
    images = []
    for result in replicates:
        record = list(result["positions"][modality][channel].values())[0]
        locs = split_by_midpoint_z(record["localisations"], split=split_z)
        n_locs = 0 if locs is None else len(locs)
        corners = count_ring_corners(
            locs, radius, r_min, r_max, n_corners, min_per_section=min_per_section,
            resolve_spots=resolve_spots,
        )
        radius_fit = None
        if n_locs >= max(min_locs_to_circlefit, 3):
            radius_fit = fit_circle_geometric(locs)[1]
            if max_radius is not None and radius_fit >= max_radius:
                radius_fit = None
        rows.append(dict(
            replicate=result["replicate"],
            n_localisations=n_locs,
            radius_fit=radius_fit,
            **corners,
        ))
        images.append(result["images"][modality][channel])
    return pd.DataFrame(rows), images, rows


def plot_ring_corners(result, axes=None):
    """
    Steps of count_ring_corners for one particle.

    (a) input, fixed-radius circle and radial window (excluded in grey);
    (b) kept localisations centred and rotated by -phi_rot;
    (c) kept localisations coloured by segment, segment edges, detected corners
        (filled circle = occupied, light circle = empty; numbers = counts).
    With resolve_spots, spot centres are drawn in (b) and (c) (black cross =
    counted as a corner, grey cross = not counted).

    Parameters
    ----------
    result : dict
        Output of count_ring_corners (or a row of ring_corners_from_replicates).
    axes : sequence of 3 matplotlib.axes.Axes, optional
        Axes to draw on. Default: a new figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if axes is None:
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    else:
        fig = axes[0].figure
    xy = result["points"][:, :2]
    kept = result["kept"]
    n = result["symmetry"]
    period = 2 * np.pi / n
    t = np.linspace(0, 2 * np.pi, 200)
    circle = np.c_[np.cos(t), np.sin(t)]

    ax = axes[0]
    ax.scatter(*xy[~kept].T, s=6, c="0.7", label="excluded")
    ax.scatter(*xy[kept].T, s=6, c="C0", label="kept")
    if result["centre"] is not None:
        c = result["centre"]
        ax.plot(*(c + result["radius"] * circle).T, "k-", lw=1, label="fixed-radius fit")
        for rr in (result["r_min"], result["r_max"]):
            ax.plot(*(c + rr * circle).T, "k:", lw=0.8)
        ax.plot(*c, "k+", ms=10)
    ax.set_title("(a) input and circle fit")
    ax.legend(fontsize=7, loc="upper right")

    al = result["aligned"]
    r_edge = result["r_max"]
    ax = axes[1]
    if len(al):
        ax.scatter(*al[kept].T, s=6, c="C0")
        ax.plot(*(result["radius"] * circle).T, "k-", lw=1)
        ax.set_title(f"(b) aligned, phi_rot = {np.degrees(result['phi_rot']):.1f} deg")
    else:
        ax.set_title("(b) aligned (no localisations kept)")

    ax = axes[2]
    if len(al):
        ax.scatter(*al[kept].T, s=8, c=result["segment"], cmap="tab10", vmin=-0.5, vmax=9.5)
        for k in range(n):
            edge = (k + 0.5) * period
            ax.plot([0, r_edge * np.cos(edge)], [0, r_edge * np.sin(edge)], "k-", lw=0.5)
            mid = k * period
            occupied = result["segment_counts"][k] > 0
            ax.plot(result["radius"] * np.cos(mid), result["radius"] * np.sin(mid),
                    marker="o", ms=12, ls="none",
                    mfc=(0, 0, 0, 0.15) if occupied else "none",
                    mec="k" if occupied else "0.6")
            ax.text(1.15 * r_edge * np.cos(mid), 1.15 * r_edge * np.sin(mid),
                    str(result["segment_counts"][k]), ha="center", va="center", fontsize=8)
    ax.set_title(f"(c) segments: {result['n_corners']} corners")

    spots = result.get("spot_means", np.empty((0, 2)))
    if len(spots) and len(al):
        c, s = np.cos(-result["phi_rot"]), np.sin(-result["phi_rot"])
        spots_aligned = (spots - result["centre"]) @ np.array([[c, s], [-s, c]])
        is_corner = result["spot_is_corner"]
        for ax in axes[1:]:
            ax.scatter(*spots_aligned[is_corner].T, marker="x", s=80, c="k", lw=2)
            ax.scatter(*spots_aligned[~is_corner].T, marker="x", s=80, c="0.6", lw=2)

    lim = 1.3 * r_edge
    for ax in axes[1:]:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
    fig.tight_layout()
    return fig


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
    xy = np.asarray(points, dtype=float)[:, :2]
    if xy.shape[0] < 2 * n_sites:
        return False, 0
    n_components = len(fit_spots(xy, n_sites, random_state=random_state)[0])
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

def simple_width(image):
    profile_ = image.sum(axis=0)
    x = np.arange(0, 500, dtype=float)
    width, center, p = _fwhm_1d(profile_, x, 5, 0, mode="total")
    return width

def filament_strip_profiles(
    image,
    pixel_size=1.0,
    strip_length=10.0,
    half_width=50.0,
    background_width=10.0,
    threshold_rel=0.05,
    flat_fraction=0.1,
    reference_percentile=90,
):
    """
    Line profiles across a straight filament, in strips along its central axis.

    The central axis is the principal axis of the pixels above threshold_rel
    of the maximum (intensity weighted), through their intensity centroid,
    so the filament can have any orientation. From the first to the last of
    these pixels along the axis, the filament is divided into strips of
    strip_length (as many as fit, centred on that extent). The pixels of each strip within half_width of the axis are
    projected on the direction across the axis: one line profile per strip
    (mean intensity per bin of one pixel). The amplitude of a profile is its
    maximum above the background (median of the outer background_width at
    each side). A strip is flat, no filament crosses it, when its amplitude
    is below flat_fraction of the reference amplitude (reference_percentile
    of the amplitudes of all strips). A run of flat strips between two strips
    that are not flat is an apparent break.

    Parameters
    ----------
    image : 2D array
        Image of a single, roughly straight filament.
    pixel_size : float
        Physical size of a pixel. strip_length, half_width, background_width
        and the lengths returned are in the same units.
    strip_length : float
        Length of a strip along the axis.
    half_width : float
        Half the length of the profiles across the axis.
    background_width : float
        Length at each end of a profile used for its background.
    threshold_rel : float
        Fraction of the maximum above which pixels define the axis and its
        extent.
    flat_fraction : float
        Fraction of the reference amplitude below which a strip is flat.
    reference_percentile : float
        Percentile of the strip amplitudes used as reference amplitude.

    Returns
    -------
    dict with keys
        "centroid", "direction" : point (px) and unit vector (row, col) of the axis
        "positions"   : centre of each strip along the axis, from the centroid
        "offsets"     : centre of each profile bin across the axis
        "profiles"    : 2D array (n_strips, n_bins), one line profile per strip
        "amplitudes"  : amplitude of each profile above its background
        "flat"        : bool per strip
        "n_strips"    : number of strips (0 for an empty image)
        "fraction_flat"  : fraction of flat strips
        "n_breaks"       : number of apparent breaks
        "longest_break"  : length of the longest apparent break (0 if none)
        "min_amplitude"  : smallest amplitude of the strips, without the first
                           and last, relative to the reference amplitude
    """
    img = np.clip(np.asarray(image, dtype=float), 0, None)
    if img.ndim != 2:
        raise ValueError("image must be 2D")
    result = dict(centroid=None, direction=None, positions=None, offsets=None,
                  profiles=None, amplitudes=None, flat=None, n_strips=0,
                  fraction_flat=np.nan, n_breaks=0, longest_break=0.0,
                  min_amplitude=np.nan)
    if img.max() <= 0:
        return result

    # central axis: principal axis of the bright pixels
    rows, cols = np.indices(img.shape)
    support = img >= threshold_rel * img.max()
    if support.sum() < 2:
        return result
    weights = img[support]
    points = np.column_stack([rows[support], cols[support]]).astype(float)
    centroid = np.average(points, axis=0, weights=weights)
    centred = points - centroid
    cov = (centred * weights[:, None]).T @ centred / weights.sum()
    direction = np.linalg.eigh(cov)[1][:, -1]
    normal = np.array([-direction[1], direction[0]])
    support_along = centred @ direction * pixel_size
    start, end = support_along.min(), support_along.max()

    # strips along the axis, centred on the extent, profiles across it
    n_strips = max(int(round((end - start) / strip_length)), 1)
    along_edges = 0.5 * (start + end - n_strips * strip_length) + strip_length * np.arange(n_strips + 1)
    n_bins = max(int(round(2 * half_width / pixel_size)), 1)
    across_edges = np.linspace(-half_width, half_width, n_bins + 1)
    offsets_px = np.column_stack([rows.ravel(), cols.ravel()]) - centroid
    along = offsets_px @ direction * pixel_size
    across = offsets_px @ normal * pixel_size
    inside = ((along >= along_edges[0]) & (along <= along_edges[-1])
              & (np.abs(across) <= half_width))
    bins = [along_edges, across_edges]
    sums = np.histogram2d(along[inside], across[inside], bins=bins,
                          weights=img.ravel()[inside])[0]
    counts = np.histogram2d(along[inside], across[inside], bins=bins)[0]
    profiles = sums / np.maximum(counts, 1)

    # amplitude above the background of each profile
    n_background = min(max(int(round(background_width / pixel_size)), 1), n_bins // 2)
    edges = np.concatenate([profiles[:, :n_background], profiles[:, -n_background:]], axis=1)
    amplitudes = np.clip(profiles.max(axis=1) - np.median(edges, axis=1), 0, None)
    reference = np.percentile(amplitudes, reference_percentile)
    flat = amplitudes < flat_fraction * reference if reference > 0 else np.ones(n_strips, bool)

    # apparent breaks: runs of flat strips with a strip that is not flat at each side
    padded = np.concatenate([[False], flat, [False]]).astype(int)
    changes = np.flatnonzero(np.diff(padded))
    run_starts, run_stops = changes[::2], changes[1::2]
    interior = (run_starts > 0) & (run_stops < n_strips)
    break_lengths = (run_stops - run_starts)[interior] * strip_length
    inner = amplitudes[1:-1] if n_strips > 2 else amplitudes

    result.update(
        centroid=centroid,
        direction=direction,
        positions=0.5 * (along_edges[:-1] + along_edges[1:]),
        offsets=0.5 * (across_edges[:-1] + across_edges[1:]),
        profiles=profiles,
        amplitudes=amplitudes,
        flat=flat,
        n_strips=n_strips,
        fraction_flat=float(flat.mean()),
        n_breaks=int(interior.sum()),
        longest_break=float(break_lengths.max()) if break_lengths.size else 0.0,
        min_amplitude=float(inner.min() / reference) if reference > 0 else np.nan,
    )
    return result
