"""Measurements on single simulated particles.

Used for the comparisons with published measurements: ring radius and
width (nuclear pore), corner detection and apparent breaks (gaps), and the
number of resolved labelled sites (PCNA).
"""

import numpy as np
from skimage.feature import peak_local_max


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
