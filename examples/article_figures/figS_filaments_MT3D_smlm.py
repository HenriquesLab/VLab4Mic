"""Fig. S: can a break be detected in an experimental microtubule, and at which labelling efficiency?

The filament is a 500 nm segment of an isolated microtubule cropped from
experimental 3D SMLM localisations (ShareLoc.xyz, Zenodo 5541063) and aligned
to the x axis. Every localisation of the segment is a target, labelled
directly (Linker probe, no normals) and imaged by SMLM. The localisations of
this segment and of the intact experimental segments (below) are fitted with
a hollow cylinder each (fit_anisotropic_hollow_cylinder): the axial spread of
the localisations is about 2.5 times the lateral one, so the shell is blurred
differently across the filament (sigma_u) and out of the plane (sigma_v).
The isotropic fit of fig2_filaments.py (fit_hollow_cylinder), which reads the
axial blur as a larger shell, is kept for comparison.

The analysis is that of fig2_filament_break_distinguishability.py with this
filament: a second model has an artificial break in the middle (every
localisation within gap/2 of the centre of the segment is removed). Each
model is simulated with run_replicates at every labelling efficiency and
every image is measured with measure_filaments (widths, mask continuity,
number of segments, eccentricities, strip profiles), giving a feature table
with one row per image. A break is detected in an image when its mask splits
into MIN_SEGMENTS_BROKEN or more segments; the share of broken filaments
where it is detected is given with a 95% Wilson interval, and the same share
of intact filaments is its false alarm rate. The table is scored per break
size and labelling efficiency with distinguishability_from_features
(out-of-fold ROC AUC with a 95% bootstrap interval, null control and
permutation p-value, for each feature set). A gap of 0 nm is a second,
independently simulated intact filament: a negative control.

Experimental segments, without simulation: segments of the top right of the
dataset, where the microtubules are isolated and continuous, are taken as
intact (EXPERIMENTAL_INTACT_CENTERS), and segments of the bottom left, where
they look broken or incomplete, as broken (EXPERIMENTAL_BROKEN_CENTERS). The
localisations of each segment are rendered as an image as the simulated ones
(2 nm pixels, Gaussian of the 5 nm lateral precision) and measured with
measure_filaments. The share of segments where a break is detected (mask
split, or a break in the strip profiles) is given per class, and the classes
are scored with distinguishability_from_features.

Outputs, in the folder <date>figS_filaments_MT3D of the VLab4Mic output
directory (created if it does not exist): segment_overview.png (the dataset
and the reference segment), filament.png (the targets of the reference
segment), hollow_cylinder_fit.png (the hollow cylinder fit of its
localisations), hollow_cylinder_fits.png (the fit of every segment, one row
each), hollow_cylinder_fit.csv (one row per segment), filament_figure.png,
the main figure:
(a) one filament per break of FIGURE_EXAMPLE_GAPS_NM imaged at every
labelling efficiency, (b) the percentage of broken filaments whose break is
detected and the false alarm rate of the intact filaments, (c) the
distinguishability (AUC) of intact and broken filaments from the number of
segments alone, and (d) the hollow cylinder fits of the reference and intact
experimental segments (pooled cross-section and lateral profile, and the
diameter of every segment with the mean +- s.d. of the fits that converged);
features.csv, distinguishability.csv, break_detection.csv,
filament_figure_values.csv, settings.json, strip_profiles.png,
distinguishability_heatmap.png, feature_distributions.png,
auc_vs_labelling.png, examples_gap<break>.png, models.png, models_grid.png
and examples.npz, as in fig2_filament_break_distinguishability.py; and for
the experimental segments experimental_segments_figure.png ((a) the segments
on the dataset, (b) their images, (c) the percentage of segments with a
detected break per class and (d) the distinguishability of the classes),
experimental_features.csv, experimental_break_detection.csv,
experimental_distinguishability.csv and experimental_segments.npz.

Run:   python figS_filaments_MT3D_smlm.py
Score an existing feature table again, without simulating:
       python figS_filaments_MT3D_smlm.py --features <folder>/features.csv
(the outputs are written next to the feature table)

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and
condition; e.g. 5 for a quick run), VLAB4MIC_N_PERMUTATIONS (200),
VLAB4MIC_N_WORKERS (1; simulations and scoring in parallel processes),
VLAB4MIC_SAVE_IMAGES (0; 1 saves the images of every model and condition as
a TIFF stack in <folder>/images), VLAB4MIC_SMLM_PATH (the .smlm file of the
localisations, local or URL; default: Zenodo).
"""

import copy
import io
import json
import os
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import tifffile
from matplotlib.patches import Polygon
from scipy import ndimage, optimize, stats
from tqdm.auto import tqdm

from vlab4mic import experiments
from vlab4mic.analysis import particle_measures
from vlab4mic.analysis.distinguishability import (
    distinguishability_from_features,
    plot_distinguishability_heatmaps,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24
OUTPUT_NAME = "figS_filaments_MT3D"  # base name of the output folder (after the date)

# Experimental 3D SMLM localisations of microtubules (ShareLoc.xyz, Zenodo 5541063)
SMLM_EXPERIMENTAL_PATH = os.environ.get(
    "VLAB4MIC_SMLM_PATH",
    "https://zenodo.org/records/5541063/files/MT3D_160528_C1C2_1424K/data.smlm")
# isolated, straight microtubule; 500 nm long as the simulated filament in fig2_filaments.py
SEGMENT_CENTER = (26500, 15200)  # nm
SEGMENT_LENGTH_NM = 500
SEGMENT_WIDTH_NM = 150
REFERENCE_DIAMETER = 25  # nm, microtubule
# localisations fitted with a hollow cylinder: within these half-widths (nm) of its axis,
# across the filament in the plane (the crop) and out of it (axial, worse precision)
CYLINDER_WINDOW_NM = (SEGMENT_WIDTH_NM / 2, 200)

STRUCTURE = "5SYF"
PROBE_MODEL = "1IGT"
SAMPLE_DIMENSIONS = [1000, 1000, 10]  # nm

LABELLING_EFFICIENCIES = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1]
GAPS_NM = [0, 5, 10, 25, 50, 100, 200]  # length of the break; 0: intact (negative control)
# breaks of the main figure, one line each in b and c; 0: the intact filament
FIGURE_EXAMPLE_GAPS_NM = [0, 50, 100]

# a break is detected when the mask splits into this many segments or more
MIN_SEGMENTS_BROKEN = 2

# example images: broken filaments whose break is detected first, then others
N_EXAMPLES = 8  # per break size and labelling efficiency
# crop of the example images around the image centre, shown with the filament horizontal
EXAMPLE_LENGTH_NM = 620  # along the filament
EXAMPLE_WIDTH_NM = 160  # across it

# experimental segments (centres in nm, as SEGMENT_CENTER), without simulation:
# isolated, continuous microtubules of the top right of the dataset ...
EXPERIMENTAL_INTACT_CENTERS = [
    (39619, 17210), (33227, 26713), (38253, 25983), (34292, 18301), (39570, 24841),
    (24532, 25702), (40022, 15971), (40319, 20646), (34518, 16741), (26789, 27774),
    (37304, 16491), (37040, 18012),
]
# ... and broken or incomplete ones of the bottom left
EXPERIMENTAL_BROKEN_CENTERS = [
    (7441, 37973), (5432, 38599), (5252, 33711), (11380, 34701), (7970, 37190),
    (4054, 28329), (10134, 35524), (2609, 30983), (11188, 33752), (2953, 37721),
    (3287, 39554), (4676, 37023),
]
RENDERING_SIGMA_NM = 5  # as the simulated images: the lateral precision

# measure_filaments columns used as features (the others are settings)
MASK_FEATURES = ["mask_continuity", "mask_major_lenght", "mask_number_of_segments",
                 "min_eccentricity", "max_eccentricity", "mean_eccentricity"]
STRIP_FEATURES = ["strip_fraction_flat", "strip_n_breaks", "strip_longest_break",
                  "strip_min_amplitude"]
WIDTHS = ["mean_filament_width", "std_filament_width", "median_filament_width", "s_width"]
FEATURES = WIDTHS + MASK_FEATURES + STRIP_FEATURES
FEATURE_SETS = {"segments": ["mask_number_of_segments"], "mask": MASK_FEATURES,
                "strip profiles": STRIP_FEATURES, "all": FEATURES}
CONTROLS = {feature: feature for feature in FEATURES}  # each alone
# images without a measurable filament (e.g. too few labels)
FILL_VALUES = dict({width: 0 for width in WIDTHS}, strip_fraction_flat=1, strip_min_amplitude=0)
CONDITIONS = ["gap_nm", "labelling_efficiency"]

COMMON = dict(
    structure=STRUCTURE,
    sample_dimensions=SAMPLE_DIMENSIONS,
    number_of_particles=1,
    modality="SMLM",
    SMLM=dict(
        lateral_precision=5,  # nm
        axial_precision=5,  # nm
        n_localisations=30,
        precision_model="fixed",  # not photon-limited
    ),
    clear_experiment=True,
    random_seed=RANDOM_SEED,
)
FILAMENTS_PARAMS = dict(
    pixel_size=2,
    max_length_nm=SEGMENT_LENGTH_NM,
    modality="SMLM",
    channel="ch0",
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=1,  # px (2 nm), Gaussian before the triangle threshold
    threshold_method="triangle",
    strip_kwargs=dict(
        strip_length=10,  # nm, along the filament axis
        half_width=50,  # nm, profiles across the axis
        background_width=10,  # nm at each end of a profile
        flat_fraction=0.1,  # flat: amplitude below this fraction of the reference
    ),
)


def load_smlm_table(path):
    """Read a .smlm file (zip with manifest.json and a binary table) into a DataFrame."""
    if path.startswith("http"):
        response = requests.get(path)
        response.raise_for_status()
        archive = zipfile.ZipFile(io.BytesIO(response.content))
    else:
        archive = zipfile.ZipFile(path)
    manifest = json.loads(archive.read("manifest.json"))
    table_file = manifest["files"][0]
    table_format = manifest["formats"][table_file["format"]]
    dtype = np.dtype([(h, "<" + np.dtype(t).str[1:]) for h, t in zip(table_format["headers"], table_format["dtype"])])
    table = np.frombuffer(archive.read(table_file["name"]), dtype=dtype, count=table_file["rows"])
    return pd.DataFrame(table)


def crop_filament_segment(localisations, center, length=500, width=150, search_radius=400):
    """
    Crop a straight filament segment and align it to the x axis.

    The filament direction is the principal axis of the localisations within
    search_radius of center (xy, nm). Localisations within length/2 along and
    width/2 across this axis are kept, rotated so the filament is horizontal
    and centred at the origin (z centred on its median).
    """
    xyz = localisations[["x", "y", "z"]].to_numpy(dtype=float)
    near = np.hypot(xyz[:, 0] - center[0], xyz[:, 1] - center[1]) < search_radius
    segment_center = xyz[near, :2].mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov((xyz[near, :2] - segment_center).T))
    direction = eigenvectors[:, -1]
    rotation = np.array([[direction[0], direction[1]], [-direction[1], direction[0]]])
    along_across = (xyz[:, :2] - segment_center) @ rotation.T
    keep = (np.abs(along_across[:, 0]) <= length / 2) & (np.abs(along_across[:, 1]) <= width / 2)
    segment = np.zeros((keep.sum(), 3))
    segment[:, :2] = along_across[keep]
    segment[:, 2] = xyz[keep, 2] - np.median(xyz[keep, 2])
    return segment, keep, segment_center, direction


def experimental_template(fluorophore_positions):
    """Targets of the experimental filament, as particle source data: without normals the
    targets are labelled directly, each localisation becomes a fluorophore at its own
    position."""
    return dict(
        targets={"ATubulin_epitopes": dict(coordinates=fluorophore_positions.copy(), normals=None)},
        scale=1e-9,
        reference_point=np.array([0.0, 0.0, 0.0]),
        axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0])),
    )


### Hollow cylinder fit of the 3D localisations

def orthonormal_basis(direction):
    """Two unit vectors orthogonal to direction and to each other."""
    direction = direction / np.linalg.norm(direction)
    helper = np.array([0.0, 0.0, 1.0]) if abs(direction[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(direction, helper)
    u /= np.linalg.norm(u)
    v = np.cross(direction, u)
    return u, v


def cylinder_coordinates(points, pivot, direction):
    """Position along the axis and (u, v) position across it, for each point."""
    u, v = orthonormal_basis(direction)
    centred = points - pivot
    return centred @ direction, centred @ u, centred @ v


def fit_hollow_cylinder(points, reference_diameter=25, background_radius_factor=3, max_points=20000, random_seed=0):
    """
    Fit a hollow cylinder to 3D coordinates (nm).

    Model: points lie on a cylindrical shell of radius R, blurred with an
    isotropic Gaussian of width sigma (labelling offset, localisation
    precision). Their distance to the axis then follows a Rice distribution
    Rice(R, sigma). Points further than background_radius_factor *
    reference_diameter from the axis are ignored, and a uniform background
    in the cross-section accounts for localisations off the filament. The
    axis is first fitted geometrically (robust least squares of distance to
    axis minus R, starting from the principal axis of the points and the
    reference radius), then axis, R, sigma and the background fraction are
    refined by maximum likelihood.

    Returns a dict with the axis (pivot, direction), radius, sigma,
    background fraction, diameter and the coordinates of every point in the
    cylinder frame (along, u, v, radial).
    """
    points = np.asarray(points, dtype=float)
    rng = np.random.default_rng(random_seed)
    if points.shape[0] > max_points:
        subset = points[rng.choice(points.shape[0], max_points, replace=False)]
    else:
        subset = points
    centroid = subset.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov((subset - centroid).T))
    direction0 = eigenvectors[:, -1]
    u0, v0 = orthonormal_basis(direction0)
    shift_limit = 2 * reference_diameter

    def axis_from_params(params):
        tilt_u, tilt_v, shift_u, shift_v = params[:4]
        direction = direction0 + tilt_u * u0 + tilt_v * v0
        direction /= np.linalg.norm(direction)
        pivot = centroid + shift_u * u0 + shift_v * v0
        return pivot, direction

    def radial_distances(params, coordinates):
        pivot, direction = axis_from_params(params)
        __, across_u, across_v = cylinder_coordinates(coordinates, pivot, direction)
        return np.hypot(across_u, across_v)

    # geometric fit: distance to the axis equal to the radius
    geometric = optimize.least_squares(
        lambda params: radial_distances(params, subset) - params[4],
        x0=[0, 0, 0, 0, reference_diameter / 2],
        bounds=([-0.5, -0.5, -shift_limit, -shift_limit, 0], [0.5, 0.5, shift_limit, shift_limit, np.inf]),
        loss="soft_l1", f_scale=5,
    )

    # maximum likelihood on the points within max_radial of the geometric
    # axis: distances to the axis follow Rice(R, sigma) plus a uniform
    # background in the cross-section
    max_radial = background_radius_factor * reference_diameter
    near_axis = subset[radial_distances(geometric.x, subset) < max_radial]
    residual_sd = np.std(radial_distances(geometric.x, near_axis) - geometric.x[4])

    def negative_log_likelihood(axis_params, radius, sigma, background):
        radial = radial_distances(axis_params, near_axis)
        shell = stats.rice.pdf(radial, radius / sigma, scale=sigma) / stats.rice.cdf(max_radial, radius / sigma, scale=sigma)
        uniform = 2 * radial / max_radial ** 2
        return -np.sum(np.log((1 - background) * shell + background * uniform + 1e-300))

    shell_bounds = [(0, max_radial), (0.5, max_radial), (0, 0.99)]
    # shell and background with the geometric axis, then all together
    shell_fit = optimize.minimize(
        lambda params: negative_log_likelihood(geometric.x[:4], *params),
        x0=[geometric.x[4], max(residual_sd, 1), 0.05],
        method="L-BFGS-B",
        bounds=shell_bounds,
    )
    likelihood = optimize.minimize(
        lambda params: negative_log_likelihood(params[:4], *params[4:]),
        x0=np.append(geometric.x[:4], shell_fit.x),
        method="L-BFGS-B",
        bounds=[(-0.5, 0.5), (-0.5, 0.5), (-shift_limit, shift_limit), (-shift_limit, shift_limit)] + shell_bounds,
    )
    pivot, direction = axis_from_params(likelihood.x)
    along, across_u, across_v = cylinder_coordinates(points, pivot, direction)
    radius, sigma, background = likelihood.x[4:]
    return dict(
        pivot=pivot,
        direction=direction,
        radius=radius,
        sigma=sigma,
        background=background,
        max_radial=max_radial,
        diameter=2 * radius,
        reference_diameter=reference_diameter,
        geometric_diameter=2 * geometric.x[4],
        along=along,
        across_u=across_u,
        across_v=across_v,
        radial=np.hypot(across_u, across_v),
        n_points=points.shape[0],
        success=bool(likelihood.success),
    )


def fit_anisotropic_hollow_cylinder(points, reference_diameter=25, window=(75, 200), n_angles=48):
    """
    Fit a hollow cylinder to 3D coordinates (nm), with a different blur across the filament
    in the plane (u) and out of it (v, close to z).

    Model: points lie on a cylindrical shell of radius R, blurred with a Gaussian of width
    sigma_u along u and sigma_v along v (labelling offset, localisation precision; the axial
    precision of 3D SMLM is worse than the lateral one, and an isotropic blur would inflate
    R). A uniform background in the cross-section accounts for localisations off the
    filament. Only the points within window = (u, v) half-widths of the initial axis (the
    principal axis of the points, through their centroid) are fitted, and the same points
    for every axis: axis (tilt and shift), R, sigma_u, sigma_v and the background fraction
    are fitted by maximum likelihood, from a few starting radii.

    Returns a dict with the axis (pivot, direction), radius, sigma_u, sigma_v, background
    fraction, diameter, the coordinates of every point in the cylinder frame (along, u, v)
    and whether the fit converged (not at a bound of radius, widths or background).
    """
    points = np.asarray(points, dtype=float)
    half_u, half_v = window
    centroid = points.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov((points - centroid).T))
    direction0 = eigenvectors[:, -1]
    u0, v0 = orthonormal_basis(direction0)
    __, across_u0, across_v0 = cylinder_coordinates(points, centroid, direction0)
    inside = (np.abs(across_u0) <= half_u) & (np.abs(across_v0) <= half_v)
    fitted = points[inside]
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)

    def axis_from_params(params):
        tilt_u, tilt_v, shift_u, shift_v = params[:4]
        direction = direction0 + tilt_u * u0 + tilt_v * v0
        direction /= np.linalg.norm(direction)
        pivot = centroid + shift_u * u0 + shift_v * v0
        return pivot, direction

    def negative_log_likelihood(params):
        radius, sigma_u, sigma_v, background = params[4:]
        __, across_u, across_v = cylinder_coordinates(fitted, *axis_from_params(params))
        shell = shell_density(across_u, across_v, radius, sigma_u, sigma_v, window, angles)
        return -np.sum(np.log((1 - background) * shell + background / (4 * half_u * half_v) + 1e-300))

    bounds = [(-0.3, 0.3), (-0.3, 0.3), (-reference_diameter, reference_diameter),
              (-4 * reference_diameter, 4 * reference_diameter),
              (0, half_u), (0.5, half_u), (0.5, half_v), (0, 0.99)]
    best = None
    for radius0 in [reference_diameter / 2, reference_diameter]:
        x0 = [0, 0, 0, 0, radius0, max(across_u0[inside].std() / 2, 1),
              max(across_v0[inside].std(), 1), 0.05]
        result = optimize.minimize(negative_log_likelihood, x0, method="L-BFGS-B", bounds=bounds)
        if best is None or result.fun < best.fun:
            best = result
    pivot, direction = axis_from_params(best.x)
    along, across_u, across_v = cylinder_coordinates(points, pivot, direction)
    radius, sigma_u, sigma_v, background = best.x[4:]
    converged = (bool(best.success) and background < 0.9 and radius < 0.95 * half_u
                 and sigma_u < 0.95 * half_u and sigma_v < 0.95 * half_v)
    return dict(
        pivot=pivot,
        direction=direction,
        radius=radius,
        sigma_u=sigma_u,
        sigma_v=sigma_v,
        background=background,
        window=window,
        diameter=2 * radius,
        reference_diameter=reference_diameter,
        along=along,
        across_u=across_u,
        across_v=across_v,
        n_points=int(inside.sum()),
        success=bool(best.success),
        converged=converged,
    )


def shell_density(across_u, across_v, radius, sigma_u, sigma_v, window, angles):
    """Density of a shell of radius radius blurred by sigma_u and sigma_v at (across_u,
    across_v), normalised within the window (u, v half-widths)."""
    half_u, half_v = window
    cos, sin = radius * np.cos(angles), radius * np.sin(angles)
    density = (stats.norm.pdf(np.asarray(across_u)[:, None], cos, sigma_u)
               * stats.norm.pdf(np.asarray(across_v)[:, None], sin, sigma_v)).mean(axis=1)
    norm = ((stats.norm.cdf((half_u - cos) / sigma_u) - stats.norm.cdf((-half_u - cos) / sigma_u))
            * (stats.norm.cdf((half_v - sin) / sigma_v) - stats.norm.cdf((-half_v - sin) / sigma_v))).mean()
    return density / norm


def cylinder_profile(fit, positions, axis="u", n_angles=48):
    """Fitted density of the points along u (across the filament, in the plane) or v (out of
    it), at positions (nm): the shell and the background, within the window."""
    half_u, half_v = fit["window"]
    radius, sigma_u, sigma_v, background = (fit[key] for key in
                                            ["radius", "sigma_u", "sigma_v", "background"])
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    cos, sin = radius * np.cos(angles), radius * np.sin(angles)
    in_u = stats.norm.cdf((half_u - cos) / sigma_u) - stats.norm.cdf((-half_u - cos) / sigma_u)
    in_v = stats.norm.cdf((half_v - sin) / sigma_v) - stats.norm.cdf((-half_v - sin) / sigma_v)
    if axis == "u":
        shell = (stats.norm.pdf(positions[:, None], cos, sigma_u) * in_v).mean(axis=1)
        half = half_u
    else:
        shell = (stats.norm.pdf(positions[:, None], sin, sigma_v) * in_u).mean(axis=1)
        half = half_v
    return (1 - background) * shell / (in_u * in_v).mean() + background / (2 * half)


def in_window(fit):
    """Points within the window of the fit."""
    half_u, half_v = fit["window"]
    return (np.abs(fit["across_u"]) <= half_u) & (np.abs(fit["across_v"]) <= half_v)


def cross_section_limits(fit):
    """Half-widths (u, v) of the cross-section shown: the shell and its blur, within the
    window."""
    half_u, half_v = fit["window"]
    return (min(fit["radius"] + 4 * fit["sigma_u"], half_u),
            min(fit["radius"] + 2 * fit["sigma_v"], half_v))


def plot_cross_section(ax, across_u, across_v, limits, radius, color, label):
    """2D histogram of the cross-section with the fitted shell and the reference diameter."""
    angles = np.linspace(0, 2 * np.pi, 100)
    limit_u, limit_v = limits
    bins = [np.arange(-limit_u, limit_u + 2.5, 2.5), np.arange(-limit_v, limit_v + 2.5, 2.5)]
    cross_section, __, __ = np.histogram2d(across_u, across_v, bins=bins)
    ax.imshow(cross_section.T, cmap="grey", origin="lower", extent=[-limit_u, limit_u, -limit_v, limit_v])
    ax.plot(radius * np.cos(angles), radius * np.sin(angles), c=color, label=label)
    reference_radius = REFERENCE_DIAMETER / 2
    ax.plot(reference_radius * np.cos(angles), reference_radius * np.sin(angles), "--", c="white",
            label=f"reference: {REFERENCE_DIAMETER:.0f} nm")
    ax.set_xlabel("u, across (nm)")
    ax.set_ylabel("v, axial (nm)")
    ax.legend(loc="upper right", fontsize="x-small")


def plot_hollow_cylinder_fit(fit, axs, color="#F58B19", title=""):
    """
    Draw an anisotropic hollow cylinder fit on three axes: the cross-section of its points
    with the fitted and reference circles, and the profiles of the points across the
    filament in the plane (u) and out of it (v, axial), with the fitted model.
    """
    inside = in_window(fit)
    plot_cross_section(axs[0], fit["across_u"][inside], fit["across_v"][inside],
                       cross_section_limits(fit), fit["radius"], color,
                       f"fit: {fit['diameter']:.1f} nm")
    axs[0].set_title(f"Cross-section{' (' + title + ')' if title else ''}")
    for ax, axis, sigma in [(axs[1], "u", "sigma_u"), (axs[2], "v", "sigma_v")]:
        half = fit["window"][0 if axis == "u" else 1]
        positions = np.linspace(-half, half, 300)
        ax.hist(fit["across_" + axis][inside], bins=np.linspace(-half, half, 60), density=True,
                color=color, alpha=0.6)
        ax.plot(positions, cylinder_profile(fit, positions, axis), c="black",
                label=f"fit: R = {fit['radius']:.1f} nm, $\\sigma_{axis}$ = {fit[sigma]:.1f} nm, "
                      f"bg = {100 * fit['background']:.0f}%")
        for sign in (-1, 1):
            ax.axvline(sign * fit["radius"], c=color)
            ax.axvline(sign * REFERENCE_DIAMETER / 2, ls="--", c="grey",
                       label=f"reference R = {REFERENCE_DIAMETER / 2:.1f} nm" if sign > 0 else None)
        ax.set_xlabel(f"{axis}, {'across the filament' if axis == 'u' else 'axial'} (nm)")
        ax.set_ylabel("density")
        ax.legend(fontsize="x-small")
    axs[1].set_title("Profile across the filament (lateral)")
    axs[2].set_title("Profile out of the plane (axial)")


### Break distinguishability (as fig2_filament_break_distinguishability.py)

def position_along_filament(coordinates):
    """Position of each point along the principal axis, from the median (nm)."""
    centred = coordinates - coordinates.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov(centred.T))
    along = centred @ eigenvectors[:, -1]
    return along - np.median(along)


def break_filament(template, gap_nm):
    """Copy of the template without the targets within gap_nm / 2 of its centre."""
    broken = copy.deepcopy(template)
    if gap_nm <= 0:
        return broken
    target = broken["targets"]["ATubulin_epitopes"]
    along = position_along_filament(target["coordinates"])
    span = np.ptp(along)
    assert 0.9 * SEGMENT_LENGTH_NM < span < 1.2 * SEGMENT_LENGTH_NM, \
        f"filament spans {span:.0f} along its axis: coordinates not in nm?"
    keep = np.abs(along) >= gap_nm / 2
    target["coordinates"] = target["coordinates"][keep]
    if target["normals"] is not None:
        target["normals"] = target["normals"][keep]
    return broken


def measure_filaments(replicates,
                      max_length_nm=500,
                      pixel_size=1,
                      modality="SMLM",
                      channel="ch0",
                      fill_holes=True,
                      smooth_sigma=1,
                      row_binning=1,
                      min_snr=2,
                      mode="single",
                      threshold_method="otsu",
                      threshold_method_kwargs=None,
                      strip_kwargs=None,):
    """Features of the first frame of every realisation, one row each: filament widths, mask
    continuity, number of mask segments and their eccentricities, and the strip profile
    measures. Also returns the width metrics, masks, smoothed images, segment tables and strip
    profiles of every image."""
    tables = []
    masks = []
    eccentricities = []
    major_axis_lengths = []
    filament_metrics_list = []
    rows = []
    blurred_list = []
    strips_list = []
    for result in replicates:
        mask, table, blurred = particle_measures.threshold_mask_and_measure(result["images"][modality][channel][0], method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = particle_measures.measure_filament_width(result["images"][modality][channel][0], pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = particle_measures.simple_width(result["images"][modality][channel][0])
        strips = particle_measures.filament_strip_profiles(result["images"][modality][channel][0], pixel_size=pixel_size, **(strip_kwargs or {}))
        strips_list.append(strips)
        filament_metrics_list.append(filament_metrics)
        blurred_list.append(blurred)
        tables.append(table)
        masks.append(mask)
        n_segments = table.shape[0]
        if n_segments == 0:
            min_eccentricity = 0
            max_eccentricity = 0
            mean_eccentricity = 0
            major_lenght = 0
            major_axis_lengths.append(major_lenght)
            apparent_continuity = 0
        else:
            min_eccentricity = np.min(table["eccentricity"])
            max_eccentricity = np.max(table["eccentricity"])
            if len(list(table["eccentricity"])) > 1:
                mean_eccentricity = np.mean(table["eccentricity"])
            else:
                mean_eccentricity = table["eccentricity"][0]
            major_lenght = table.iloc[table['major_axis_length'].idxmax()]["major_axis_length"]
            major_axis_lengths.append(major_lenght)
            apparent_continuity = (major_lenght / max_length_nm) / n_segments

        rows.append(dict(
            replicate=result["replicate"],
            mean_filament_width=filament_metrics["mean"],
            std_filament_width=filament_metrics["std"],
            median_filament_width=filament_metrics["median"],
            mask_continuity=apparent_continuity,
            mask_major_lenght=major_lenght,
            mask_number_of_segments=n_segments,
            mask_prefilter_smooth_sigma=smooth_sigma,
            pixelsize=pixel_size,
            row_binning=row_binning,
            s_width=s_width,
            min_eccentricity=min_eccentricity,
            max_eccentricity=max_eccentricity,
            mean_eccentricity=mean_eccentricity,
            strip_fraction_flat=strips["fraction_flat"],
            strip_n_breaks=strips["n_breaks"],
            strip_longest_break=strips["longest_break"],
            strip_min_amplitude=strips["min_amplitude"],
        ))
    df = pd.DataFrame(rows)
    return df, filament_metrics_list, masks, blurred_list, tables, strips_list


def simulate_replicates(labelling_efficiency, seed, particle_data, n_replicates=N_REALISATIONS):
    """Independent realisations of a filament (run_replicates output): its targets labelled
    directly (Linker) with labelling_efficiency."""
    primary = dict(
        probe_template="Linker",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope=0,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        probe_steric_hindrance=0,  # localisations are not excluded by probe size
        labelling_efficiency=labelling_efficiency,
    )
    __, __, experiment = experiments.image_vsample(
        probe_list=[primary],
        **dict(COMMON, random_seed=seed),
    )
    experiment.particle.source = dict()
    experiment.particle.load_source(**copy.deepcopy(particle_data))
    experiment.set_virtualsample_params(
        sample_dimensions=SAMPLE_DIMENSIONS,
        update_mode=False,
    )
    experiment.particle.generate_instance()
    experiment.build(modules=["coordinate_field", "imager"])
    return experiment.run_replicates(n_replicates=n_replicates, modality="SMLM")


def break_detected(table):
    """Rows where the mask splits into MIN_SEGMENTS_BROKEN or more segments."""
    return table.mask_number_of_segments >= MIN_SEGMENTS_BROKEN


def crop_filament(image):
    """The image centre, EXAMPLE_LENGTH_NM along the filament by EXAMPLE_WIDTH_NM across,
    with the filament horizontal: it runs along the rows of the image, so the crop is
    transposed."""
    image = np.asarray(image, dtype=float)
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    rows, cols = np.array(image.shape) // 2
    half_length = int(round(EXAMPLE_LENGTH_NM / pixel_size / 2))
    half_width = int(round(EXAMPLE_WIDTH_NM / pixel_size / 2))
    crop = image[rows - half_length:rows + half_length, cols - half_width:cols + half_width]
    return crop.T.astype(np.float32)


EXTENT = [-EXAMPLE_LENGTH_NM / 2, EXAMPLE_LENGTH_NM / 2, -EXAMPLE_WIDTH_NM / 2,
          EXAMPLE_WIDTH_NM / 2]


def show(ax, image):
    """A cropped filament image, without ticks."""
    ax.imshow(image, cmap="magma", extent=EXTENT, origin="lower")
    ax.set_xticks([])
    ax.set_yticks([])


def scale_bar(ax, colour="w"):
    """100 nm, at the lower left of a cropped filament image."""
    x0, y0 = -EXAMPLE_LENGTH_NM / 2 + 15, -EXAMPLE_WIDTH_NM / 2 + 15
    ax.plot([x0, x0 + 100], [y0, y0], "-", color=colour, lw=2)


def model_key(model, gap, le):
    return f"model_{model}_gap_{gap:g}_labelling_efficiency_{le:g}"


def epitope_key(model, gap):
    return f"model_{model}_gap_{gap:g}_epitopes"


def example_key(gap, le):
    return f"gap_{gap:g}_labelling_efficiency_{le:g}"


def model_label(model, gap):
    return "intact" if model == "intact" else f"break {gap:g} nm"


def wilson_interval(k, n, z=1.96):
    """95% interval of a proportion k / n."""
    p = k / n
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return max(centre - half, 0.0), min(centre + half, 1.0)


def error_bars(value, low, high):
    """Lengths of the error bars below and above value, never negative (rounding can put
    value a hair outside [low, high], which matplotlib refuses)."""
    value, low, high = (np.asarray(x, dtype=float) for x in (value, low, high))
    return [np.clip(value - low, 0, None), np.clip(high - value, 0, None)]


def labelling_axis(ax):
    """Labelling efficiency on a log scale, ticks at the simulated values."""
    les = sorted(LABELLING_EFFICIENCIES)
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks(les, [f"{le:g}" for le in les], fontsize=8)
    ax.set_xlabel("Labelling efficiency")
    ax.spines[["top", "right"]].set_visible(False)


def simulate(job):
    """One model at one labelling efficiency: its feature table (one row per image), the
    first image (cropped) and, for the broken filament, N_EXAMPLES images (cropped), those
    whose break is detected first."""
    model, gap_nm, labelling_efficiency, seed, particle_data, images_dir = job
    replicates = simulate_replicates(labelling_efficiency, seed, particle_data)
    table = measure_filaments(replicates=replicates, **FILAMENTS_PARAMS)[0]
    table = table.assign(model=model, gap_nm=gap_nm, labelling_efficiency=labelling_efficiency,
                         seed=seed)
    # the frame measured by measure_filaments, one per realisation
    images = [np.asarray(result["images"]["SMLM"]["ch0"][0]) for result in replicates]
    if images_dir is not None:
        os.makedirs(images_dir, exist_ok=True)
        tifffile.imwrite(
            os.path.join(images_dir, f"{model}_gap_{gap_nm}_labelling_efficiency_"
                                     f"{labelling_efficiency}.tif"),
            np.array(images, dtype=np.float32))
    examples = {model_key(model, gap_nm, labelling_efficiency): crop_filament(images[0])}
    if model == "broken":
        replicate = table["replicate"].astype(int)
        detected = break_detected(table)
        replicas = pd.concat([replicate[detected], replicate[~detected]]).to_numpy()[:N_EXAMPLES]
        key = example_key(gap_nm, labelling_efficiency)
        examples[key] = np.array([crop_filament(images[r]) for r in replicas])
        examples[key + "_replicas"] = replicas
    return table, examples


def model_epitopes(models):
    """Targets of every model in the plane, from the centre of the intact filament (nm)."""
    centre = models[0][2]["targets"]["ATubulin_epitopes"]["coordinates"].mean(axis=0)
    return {epitope_key(model, gap):
            data["targets"]["ATubulin_epitopes"]["coordinates"][:, :2] - centre[:2]
            for model, gap, data in models}


def models_shown():
    """(model, break) of the intact filament and of every broken one."""
    return [("intact", 0)] + [("broken", gap) for gap in GAPS_NM]


def plot_models(examples, prefix):
    """The models, one per column: their targets from above (removed ones open) and an
    image at every labelling efficiency (rows)."""
    models = models_shown()
    les = sorted(LABELLING_EFFICIENCIES)
    all_epitopes = examples[epitope_key("intact", 0)]
    fig, axes = plt.subplots(1 + len(les), len(models), squeeze=False,
                             figsize=[2.4 * len(models), 0.75 * (1 + len(les)) + 0.6])
    for column, (model, gap) in zip(axes.T, models):
        epitopes = examples[epitope_key(model, gap)]
        column[0].scatter(*all_epitopes.T, s=0.5, facecolors="none", edgecolors="0.7",
                          linewidths=0.3)
        column[0].scatter(*epitopes.T, s=0.5, color="k")
        column[0].set(xlim=EXTENT[:2], ylim=EXTENT[2:], aspect="equal", xticks=[], yticks=[])
        column[0].set_title(f"{model_label(model, gap)}\n{len(epitopes)} localisations", fontsize=8)
        for ax, le in zip(column[1:], les):
            show(ax, examples[model_key(model, gap, le)])
    scale_bar(axes[0, 0], colour="k")
    axes[0, 0].set_ylabel("Targets", fontsize=8)
    for ax, le in zip(axes[1:, 0], les):
        scale_bar(ax)
        ax.set_ylabel(f"LE {le:g}", fontsize=8)
    fig.suptitle("Experimental microtubule models: targets (experimental localisations; black, "
                 "kept; open, removed) and one image at every labelling efficiency; scale bar "
                 "100 nm", x=0.02, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(prefix + "models.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strip_examples(template, prefix, labelling_efficiency=1, gaps=(0, 50, 100)):
    """Mask and strip profiles of one image of the intact and broken filaments."""
    flat_fraction = FILAMENTS_PARAMS["strip_kwargs"]["flat_fraction"]
    fig, axes = plt.subplots(len(gaps), 4, figsize=(16, 3.2 * len(gaps)),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.4, 2]))
    for row, gap in zip(axes, gaps):
        replicates = simulate_replicates(labelling_efficiency, RANDOM_SEED + gap,
                                         break_filament(template, gap), n_replicates=1)
        table, __, masks, __, __, strips = measure_filaments(replicates, **FILAMENTS_PARAMS)
        image, mask, strip = replicates[0]["images"]["SMLM"]["ch0"][0], masks[0], strips[0]
        plot_strip_row(row, image, mask, strip, table, f"gap {gap} nm, LE {labelling_efficiency}",
                       flat_fraction)
    fig.tight_layout()
    fig.savefig(prefix + "strip_profiles.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strip_row(row, image, mask, strip, table, title, flat_fraction):
    """Image and mask around the filament, with the axis and the flat strips, the line
    profile of every strip and their amplitude along the filament, on four axes."""
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    rows, cols = np.nonzero(image > 0.05 * image.max())
    crop = (slice(max(rows.min() - 30, 0), rows.max() + 30),
            slice(max(cols.min() - 30, 0), cols.max() + 30))
    origin = np.array([crop[0].start, crop[1].start])
    row[0].imshow(image[crop], cmap="grey", vmax=np.percentile(image[image > 0], 99.5))
    axis_points = (strip["centroid"] - origin
                   + np.outer(strip["positions"] / pixel_size, strip["direction"]))
    row[0].plot(axis_points[:, 1], axis_points[:, 0], c="C1", lw=0.8)
    row[0].scatter(axis_points[strip["flat"], 1], axis_points[strip["flat"], 0], c="red", s=6)
    row[0].set_title(title, loc="left", fontsize=9)
    row[1].imshow(mask[crop], cmap="grey")
    row[1].set_title(f"mask: {int(table.mask_number_of_segments[0])} segment(s)",
                     loc="left", fontsize=9)
    for ax in row[:2]:
        ax.axis("off")
    # one line profile per strip
    extent = [strip["offsets"][0], strip["offsets"][-1],
              strip["positions"][-1], strip["positions"][0]]
    row[2].imshow(strip["profiles"], cmap="grey", aspect="auto", extent=extent)
    row[2].set_xlabel("across the axis (nm)")
    row[2].set_ylabel("along the axis (nm)")
    row[2].set_title("line profile of each strip", loc="left", fontsize=9)
    # amplitude of the profiles along the filament
    reference = strip["amplitudes"].max() if strip["amplitudes"].max() > 0 else 1
    row[3].plot(strip["positions"], strip["amplitudes"] / reference, "o-", c="k", ms=2, lw=0.8)
    row[3].scatter(strip["positions"][strip["flat"]],
                   strip["amplitudes"][strip["flat"]] / reference, c="red", s=12, zorder=3,
                   label="flat")
    row[3].axhline(flat_fraction * np.percentile(strip["amplitudes"], 90) / reference,
                   ls="--", c="grey", lw=0.8, label="flat threshold")
    row[3].set_xlabel("along the axis (nm)")
    row[3].set_ylabel("profile amplitude (norm.)")
    row[3].set_title(f"{strip['n_breaks']} break(s), longest {strip['longest_break']:.0f} nm",
                     loc="left", fontsize=9)
    row[3].legend(frameon=False, fontsize=8)


def feature_table(tables):
    """Intact and broken rows for every break size."""
    features = pd.concat(tables, ignore_index=True)
    intact = features[features.model == "intact"]
    broken = features[features.model == "broken"]
    # the intact images are shared by every break size
    return pd.concat([broken] + [intact.assign(gap_nm=gap) for gap in GAPS_NM],
                     ignore_index=True)


def break_detection(features):
    """Per model, break size and labelling efficiency: the share of filaments whose mask splits
    into MIN_SEGMENTS_BROKEN or more segments (broken: detected breaks; intact: false alarms),
    and the median, mean and s.d. of their number of segments."""
    intact = features[features.model == "intact"]
    # the intact filaments appear once per break size: keep one copy
    intact = intact[intact.gap_nm == intact.gap_nm.min()].assign(gap_nm=0)
    rows = []
    for (model, gap, le), data in pd.concat([intact, features[features.model == "broken"]]) \
            .groupby(["model", "gap_nm", "labelling_efficiency"]):
        split = int(break_detected(data).sum())
        low, high = wilson_interval(split, len(data))
        rows.append(dict(model=model, gap_nm=gap, labelling_efficiency=le,
                         n_filaments=len(data), n_split=split, fraction_split=split / len(data),
                         fraction_low=low, fraction_high=high,
                         median_segments=data.mask_number_of_segments.median(),
                         mean_segments=data.mask_number_of_segments.mean(),
                         std_segments=data.mask_number_of_segments.std()))
    return pd.DataFrame(rows).sort_values(["model", "gap_nm", "labelling_efficiency"],
                                          ascending=[False, True, True], ignore_index=True)


def figure_values(detection, scores):
    """The values of b and c of the main figure, one row per model (the intact filament and every
    break of FIGURE_EXAMPLE_GAPS_NM) and labelling efficiency: b, the percentage of filaments whose
    mask splits into MIN_SEGMENTS_BROKEN or more segments (95% Wilson interval), with the mean and
    s.d. of their number of segments; c, the AUC of intact against broken filaments from the
    number of segments alone (95% bootstrap interval), its null control (intact against intact)
    and permutation p-value."""
    breaks = [gap for gap in FIGURE_EXAMPLE_GAPS_NM if gap > 0]
    shown = detection[(detection.model == "intact")
                      | ((detection.model == "broken") & detection.gap_nm.isin(breaks))]
    values = pd.DataFrame(dict(
        model=shown.model, gap_nm=shown.gap_nm, labelling_efficiency=shown.labelling_efficiency,
        n_filaments=shown.n_filaments, b_percent_detected=100 * shown.fraction_split,
        b_percent_low=100 * shown.fraction_low, b_percent_high=100 * shown.fraction_high,
        b_mean_segments=shown.mean_segments, b_std_segments=shown.std_segments))
    segments = scores.loc[scores.feature_set == "segments",
                          CONDITIONS + ["auc", "auc_low", "auc_high", "null_auc", "p_value"]]
    segments = segments.rename(columns=lambda column: column if column in CONDITIONS
                               else "c_" + column)
    # c compares the intact filament with each break: no AUC on the intact rows
    return values.merge(segments.assign(model="broken"), how="left", on=["model"] + CONDITIONS)


def plot_feature_distributions(features, prefix):
    """Mean +- s.d. of the continuity features against labelling efficiency, one line per break
    size (a break of 0 nm: the second intact filament)."""
    broken = features[features.model == "broken"]
    gaps = sorted(broken.gap_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(gaps)))
    columns = ["mask_number_of_segments", "mask_continuity", "strip_n_breaks",
               "strip_longest_break", "strip_fraction_flat", "mean_filament_width"]
    fig, axes = plt.subplots(1, len(columns), figsize=[3.5 * len(columns), 3])
    for ax, feature in zip(axes, columns):
        for gap, colour in zip(gaps, colours):
            summary = broken[np.isclose(broken.gap_nm, gap)] \
                .groupby("labelling_efficiency")[feature].agg(["mean", "std"])
            ax.errorbar(summary.index, summary["mean"], yerr=summary["std"], fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{gap:g}")
        labelling_axis(ax)
        ax.set_ylabel(feature.replace("_", " "))
    axes[-1].legend(title="Break (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "feature_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_auc_vs_labelling(scores, prefix):
    """AUC (95% interval) against labelling efficiency, one line per break size."""
    gaps = sorted(scores.gap_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(gaps)))
    fig, axes = plt.subplots(1, len(FEATURE_SETS), figsize=[3.5 * len(FEATURE_SETS), 3],
                             sharey=True)
    for ax, name in zip(axes, FEATURE_SETS):
        data = scores[scores.feature_set == name]
        for gap, colour in zip(gaps, colours):
            sub = data[np.isclose(data.gap_nm, gap)].sort_values("labelling_efficiency")
            ax.errorbar(sub.labelling_efficiency, sub.auc,
                        yerr=error_bars(sub.auc, sub.auc_low, sub.auc_high), fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{gap:g}")
        ax.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
        ax.set_ylim(0.3, 1.05)
        labelling_axis(ax)
        ax.set_title(name, loc="left")
    axes[0].set_ylabel("ROC AUC (intact vs broken)")
    axes[-1].legend(title="Break (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "auc_vs_labelling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_examples(features, examples, prefix):
    """Per break size: example broken filaments, one column per labelling efficiency, those
    whose break is detected first (green frame), each with its segments and longest break."""
    broken = features[features.model == "broken"]
    les = sorted(LABELLING_EFFICIENCIES)
    for gap in GAPS_NM:
        # fewer rows when there are fewer realisations than N_EXAMPLES
        n_rows = max([len(examples.get(example_key(gap, le), [])) for le in les] + [1])
        fig, axes = plt.subplots(n_rows, len(les), squeeze=False,
                                 figsize=[2.4 * len(les), 0.75 * n_rows + 0.8])
        for column, le in zip(axes.T, les):
            key = example_key(gap, le)
            images = examples.get(key, [])
            replicas = examples.get(key + "_replicas", [])
            data = broken[np.isclose(broken.gap_nm, gap)
                          & np.isclose(broken.labelling_efficiency, le)]
            measured = data.set_index(data.replicate.astype(int))
            for k, ax in enumerate(column):
                if k >= len(images):
                    ax.axis("off")
                    continue
                show(ax, images[k])
                if k < len(replicas) and replicas[k] in measured.index:
                    annotate_break(ax, measured.loc[replicas[k]])
            if len(images):
                scale_bar(column[0])
            column[0].set_title(f"LE {le:g}\n{int(break_detected(data).sum())} / {len(data)} "
                                f"split", fontsize=9)
        fig.suptitle(f"Experimental microtubule, break of {gap:g} nm: labelling efficiency per "
                     f"column; green, break detected ({MIN_SEGMENTS_BROKEN} or more mask "
                     f"segments); scale bar 100 nm", x=0.02, ha="left", fontsize=10)
        fig.tight_layout()
        fig.savefig(prefix + f"examples_gap{gap:g}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def annotate_break(ax, found):
    """Segments and longest break of a filament image, framed green when its break is
    detected."""
    detected = found.mask_number_of_segments >= MIN_SEGMENTS_BROKEN
    ax.text(0.02, 0.95, f"{found.mask_number_of_segments:.0f} segment(s), "
                        f"longest break {found.strip_longest_break:.0f} nm",
            color="C2" if detected else "w", fontsize=5, ha="left", va="top",
            transform=ax.transAxes)
    if detected:
        ax.spines[:].set_color("C2")
        ax.spines[:].set_linewidth(2.5)


def plot_models_grid(examples, scores, prefix):
    """The pair of images compared in every cell of the distinguishability map (intact above
    broken), laid out as the map, with the AUC of all features."""
    gaps = sorted(GAPS_NM, reverse=True)  # as the map: smallest break at the bottom
    les = sorted(LABELLING_EFFICIENCIES)
    ratios = ([1, 1, 0.35] * len(gaps))[:-1]  # a pair of panels per break, then a gap
    fig, axes = plt.subplots(len(ratios), len(les), squeeze=False,
                             figsize=[2.0 * len(les) + 0.8, 0.52 * sum(ratios) + 0.8],
                             gridspec_kw=dict(height_ratios=ratios, wspace=0.05, hspace=0.08))
    auc = scores[scores.feature_set == "all"]
    for i, gap in enumerate(gaps):
        intact_row, broken_row = axes[3 * i], axes[3 * i + 1]
        for k, le in enumerate(les):
            show(intact_row[k], examples[model_key("intact", 0, le)])
            show(broken_row[k], examples[model_key("broken", gap, le)])
            cell = auc[np.isclose(auc.gap_nm, gap) & np.isclose(auc.labelling_efficiency, le)]
            if len(cell):
                broken_row[k].text(0.98, 0.9, f"AUC {cell.auc.iloc[0]:.2f}", color="w",
                                   fontsize=6, ha="right", va="top",
                                   transform=broken_row[k].transAxes)
            if i == 0:
                intact_row[k].set_title(f"LE {le:g}", fontsize=8)
        intact_row[0].set_ylabel("intact", fontsize=6)
        broken_row[0].set_ylabel(f"{gap:g} nm", fontsize=6)
    for row in axes[2::3]:
        for ax in row:
            ax.axis("off")
    scale_bar(axes[-1, 0])
    fig.supylabel("Break length (nm)", fontsize=9)
    fig.supxlabel("Labelling efficiency", fontsize=9)
    fig.suptitle("Experimental microtubule: filaments compared in each cell of the "
                 "distinguishability map (intact above broken; AUC, all features); scale bar "
                 "100 nm", x=0.02, ha="left", fontsize=10)
    fig.savefig(prefix + "models_grid.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_figure(examples, detection, scores, fits, prefix):
    """The main figure: (a) one filament per break of FIGURE_EXAMPLE_GAPS_NM (rows; 0 nm: the intact
    filament) at every labelling efficiency, (b) the share of broken filaments whose break is
    detected, and the false alarm rate of intact filaments, (c) the AUC of intact against
    broken filaments from the number of mask segments alone, both against labelling
    efficiency, and (d) the hollow cylinder fits of the reference and intact experimental
    segments."""
    les = sorted(LABELLING_EFFICIENCIES)
    breaks = [gap for gap in FIGURE_EXAMPLE_GAPS_NM if gap > 0]  # one line each in b and c
    colours = {gap: f"C{k}" for k, gap in enumerate(breaks)}
    colours[0] = "0.5"  # the intact filament, as in b
    fig = plt.figure(figsize=[11, 6.8 + 0.55 * len(FIGURE_EXAMPLE_GAPS_NM)])
    grid = fig.add_gridspec(3, 1, height_ratios=[0.25 * len(FIGURE_EXAMPLE_GAPS_NM) + 0.15, 1, 1],
                            hspace=0.55)
    top = grid[0].subgridspec(len(FIGURE_EXAMPLE_GAPS_NM), len(les), wspace=0.05, hspace=0.08)
    for i, gap in enumerate(FIGURE_EXAMPLE_GAPS_NM):
        model = "broken" if gap > 0 else "intact"
        for k, le in enumerate(les):
            ax = fig.add_subplot(top[i, k])
            show(ax, examples[model_key(model, gap, le)])
            if i == 0:
                ax.set_title(f"LE {le:g}", fontsize=9)
            if k == 0:
                ax.set_ylabel(f"{gap:g} nm" if gap > 0 else "no break", fontsize=8,
                              color=colours[gap])
            if i == 0 and k == 0:
                ax.set_title("a", loc="left", fontweight="bold", x=-0.12)
            if i == len(FIGURE_EXAMPLE_GAPS_NM) - 1 and k == 0:
                scale_bar(ax)
    bottom = grid[1].subgridspec(1, 2, wspace=0.3)

    ax = fig.add_subplot(bottom[0])
    for gap in breaks:
        data = detection[(detection.model == "broken") & np.isclose(detection.gap_nm, gap)] \
            .sort_values("labelling_efficiency")
        percent = 100 * data.fraction_split
        ax.errorbar(data.labelling_efficiency, percent,
                    yerr=error_bars(percent, 100 * data.fraction_low, 100 * data.fraction_high),
                    fmt="o-", color=colours[gap], capsize=3, label=f"{gap:g} nm break")
    intact = detection[detection.model == "intact"].sort_values("labelling_efficiency")
    percent = 100 * intact.fraction_split
    ax.errorbar(intact.labelling_efficiency, percent,
                yerr=error_bars(percent, 100 * intact.fraction_low, 100 * intact.fraction_high),
                fmt="s:", color="0.5", capsize=3, markersize=4, label="intact (false positive)")
    ax.set_ylim(-3, 103)
    ax.set_ylabel("Break detected (% of filaments)")
    ax.set_title("b", loc="left", fontweight="bold", x=-0.22)
    ax.set_title(f"Mask split into {MIN_SEGMENTS_BROKEN}+ segments "
                 f"(n = {detection.n_filaments.max()} per point)", fontsize=9)
    ax.legend(frameon=False, fontsize=7)

    ax_c = fig.add_subplot(bottom[1])
    data = scores[scores.feature_set == "segments"]
    for k, gap in enumerate(breaks):
        sub = data[np.isclose(data.gap_nm, gap)].sort_values("labelling_efficiency")
        ax_c.errorbar(sub.labelling_efficiency, sub.auc,
                      yerr=error_bars(sub.auc, sub.auc_low, sub.auc_high), fmt="o-",
                      color=colours[gap], capsize=3, label=f"{gap:g} nm break")
        ax_c.plot(sub.labelling_efficiency, sub.null_auc, "x", color=colours[gap],
                  label="null (intact vs intact)" if k == 0 else None)
    ax_c.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_c.set_ylim(min(0.3, np.nanmin(data.null_auc) - 0.05), 1.05)
    ax_c.set_ylabel("ROC AUC")
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.22)
    ax_c.set_title("Intact vs broken, number of mask segments", fontsize=9)
    ax_c.legend(frameon=False, fontsize=7, loc="lower right")
    for axis in (ax, ax_c):
        labelling_axis(axis)

    cylinder = grid[2].subgridspec(1, 3, wspace=0.35)
    axs = [fig.add_subplot(cylinder[k]) for k in range(3)]
    plot_cylinder_summary(fits, axs, color="#F58B19")
    for axis in axs:
        axis.title.set_fontsize(9)
    axs[0].set_title("d", loc="left", fontweight="bold", x=-0.3)
    fig.savefig(prefix + "filament_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


### Experimental segments: intact (top right) vs broken (bottom left), without simulation

def render_segment(segment):
    """Image of the localisations of an aligned segment, as the simulated images: a field of
    SAMPLE_DIMENSIONS with FILAMENTS_PARAMS["pixel_size"] pixels, the filament centred along
    the rows, each localisation a Gaussian of RENDERING_SIGMA_NM."""
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    edges = [np.arange(-size / 2, size / 2 + pixel_size / 2, pixel_size)
             for size in SAMPLE_DIMENSIONS[:2]]
    image, __, __ = np.histogram2d(segment[:, 0], segment[:, 1], bins=edges)
    return ndimage.gaussian_filter(image, RENDERING_SIGMA_NM / pixel_size)


def segment_outline(center, direction, length=SEGMENT_LENGTH_NM, width=SEGMENT_WIDTH_NM):
    """Corners (xy, nm) of a cropped segment in the dataset."""
    across = np.array([-direction[1], direction[0]])
    return np.array([center + a * length / 2 * direction + b * width / 2 * across
                     for a, b in [(-1, -1), (1, -1), (1, 1), (-1, 1)]])


def measure_experimental_segments(localisations_all):
    """The intact and broken experimental segments: their features (measure_filaments, one row
    per segment), and their images (cropped), centres, directions, number of localisations
    and aligned localisations (points, with the index of their segment in point_segment)."""
    rows, images, centres, directions, n_localisations, models = [], [], [], [], [], []
    points, point_segment = [], []
    for model, segment_centres in [("intact", EXPERIMENTAL_INTACT_CENTERS),
                                   ("broken", EXPERIMENTAL_BROKEN_CENTERS)]:
        replicates = []
        for k, centre in enumerate(segment_centres):
            segment, __, segment_centre, direction = crop_filament_segment(
                localisations_all, center=centre, length=SEGMENT_LENGTH_NM, width=SEGMENT_WIDTH_NM)
            image = render_segment(segment)
            replicates.append(dict(replicate=k, images={"SMLM": {"ch0": [image]}}))
            images.append(crop_filament(image))
            point_segment.append(np.full(len(segment), len(centres)))
            points.append(segment)
            centres.append(segment_centre)
            directions.append(direction)
            n_localisations.append(len(segment))
            models.append(model)
        table = measure_filaments(replicates=replicates, **FILAMENTS_PARAMS)[0]
        rows.append(table.assign(model=model, segment_x=[c[0] for c in segment_centres],
                                 segment_y=[c[1] for c in segment_centres],
                                 n_localisations=n_localisations[-len(segment_centres):]))
    features = pd.concat(rows, ignore_index=True)
    segments = dict(images=np.array(images), centres=np.array(centres),
                    directions=np.array(directions), n_localisations=np.array(n_localisations),
                    models=np.array(models), points=np.vstack(points),
                    point_segment=np.concatenate(point_segment))
    return features, segments


def experimental_break_detection(features):
    """Per class of experimental segments: the share whose mask splits into MIN_SEGMENTS_BROKEN
    or more segments, and the share with a break in the strip profiles (95% Wilson
    intervals), with the mean and s.d. of their segments, strip breaks and longest break."""
    rows = []
    for model in ("intact", "broken"):
        data = features[features.model == model]
        row = dict(model=model, n_segments=len(data))
        for name, found in [("mask_split", break_detected(data)),
                            ("strip_break", data.strip_n_breaks >= 1)]:
            k = int(found.sum())
            low, high = wilson_interval(k, len(data))
            row.update({f"n_{name}": k, f"fraction_{name}": k / len(data),
                        f"fraction_{name}_low": low, f"fraction_{name}_high": high})
        for column in ["mask_number_of_segments", "strip_n_breaks", "strip_longest_break",
                       "n_localisations"]:
            row[f"mean_{column}"] = data[column].mean()
            row[f"std_{column}"] = data[column].std()
        rows.append(row)
    return pd.DataFrame(rows)


def experimental_distinguishability(features):
    """Distinguishability of the intact and broken experimental segments, per feature set."""
    scores = []
    for name, columns in FEATURE_SETS.items():
        table = distinguishability_from_features(
            features, "model", ("intact", "broken"), columns,
            controls=CONTROLS if name == "all" else None,
            replicate_column="replicate", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, progress=False,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    return pd.concat(scores, ignore_index=True)


def plot_experimental_figure(features, segments, detection, scores, prefix):
    """(a) The intact and broken experimental segments (and the reference segment) on the
    dataset, (b) their images, with their segments and longest break, framed green when the
    break is detected, (c) the share of segments with a detected break per class and (d) the
    distinguishability of the classes per feature set."""
    colours = {"intact": "#2b8cbe", "broken": "#d7301f", "reference": "#F58B19"}
    models = segments["models"]
    fig = plt.figure(figsize=[15, 11])
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.5], height_ratios=[1.6, 1],
                            hspace=0.25, wspace=0.15)

    ax = fig.add_subplot(grid[0, 0])
    overview, extent = segments["overview"], segments["overview_extent"]
    ax.imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
              extent=extent)
    outlines = [("reference", 0, segments["reference_centre"], segments["reference_direction"])]
    counters = {"intact": 0, "broken": 0}
    for model, centre, direction in zip(models, segments["centres"], segments["directions"]):
        outlines.append((model, counters[model], centre, direction))
        counters[model] += 1
    for model, k, centre, direction in outlines:
        ax.add_patch(Polygon(segment_outline(centre, direction), closed=True, fill=False,
                             edgecolor=colours[model], linewidth=1.2))
        ax.text(centre[0] + 400, centre[1] - 400, "ref" if model == "reference" else f"{k}",
                color=colours[model], fontsize=6)
    for model in colours:
        ax.plot([], [], "s", mfc="none", color=colours[model], label=model)
    ax.legend(frameon=False, fontsize=7, labelcolor="w", loc="lower right")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.tick_params(labelsize=7)
    ax.set_title("a", loc="left", fontweight="bold")
    ax.set_title("Experimental segments (500 x 150 nm)", fontsize=9)

    # intact segments above, broken below
    n_columns = 3
    n_rows = int(np.ceil(max(counters.values()) / n_columns))
    gallery = grid[0, 1].subgridspec(2 * n_rows + 1, n_columns, wspace=0.05, hspace=0.25,
                                     height_ratios=[1] * n_rows + [0.4] + [1] * n_rows)
    for model, offset in [("intact", 0), ("broken", n_rows + 1)]:
        index = np.flatnonzero(models == model)
        measured = features[features.model == model].set_index("replicate")
        for k, i in enumerate(index):
            ax_k = fig.add_subplot(gallery[offset + k // n_columns, k % n_columns])
            show(ax_k, segments["images"][i])
            annotate_break(ax_k, measured.loc[k])
            ax_k.text(0.98, 0.05, f"{k}", color=colours[model], fontsize=6, ha="right",
                      transform=ax_k.transAxes)
            if k == 0:
                ax_k.set_title(f"{model} ({len(index)})", fontsize=9, loc="left",
                               color=colours[model])
                scale_bar(ax_k)
                if model == "intact":
                    ax_k.text(-0.08, 1.15, "b", fontweight="bold", fontsize=12, ha="right",
                              transform=ax_k.transAxes)

    bottom = grid[1, :].subgridspec(1, 2, wspace=0.3)
    ax_c = fig.add_subplot(bottom[0])
    width = 0.38
    for k, (name, label) in enumerate([("mask_split", f"mask split ({MIN_SEGMENTS_BROKEN}+ segments)"),
                                       ("strip_break", "strip profiles (1+ breaks)")]):
        positions = np.arange(len(detection)) + (k - 0.5) * width
        percent = 100 * detection[f"fraction_{name}"]
        ax_c.bar(positions, percent, width, color=[colours[m] for m in detection.model],
                 alpha=1 if k == 0 else 0.5, hatch=None if k == 0 else "//",
                 edgecolor="k", linewidth=0.5, label=label,
                 yerr=error_bars(percent, 100 * detection[f"fraction_{name}_low"],
                                 100 * detection[f"fraction_{name}_high"]), capsize=3)
    ax_c.set_xticks(np.arange(len(detection)),
                    [f"{m}\n(n = {n})" for m, n in zip(detection.model, detection.n_segments)])
    ax_c.set_ylim(0, 105)
    ax_c.set_ylabel("Break detected (% of segments)")
    ax_c.legend(frameon=False, fontsize=7, loc="upper left")
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.15)

    ax_d = fig.add_subplot(bottom[1])
    positions = np.arange(len(scores))
    ax_d.errorbar(positions, scores.auc, yerr=error_bars(scores.auc, scores.auc_low, scores.auc_high),
                  fmt="o", color="k", capsize=3, label="intact vs broken")
    ax_d.plot(positions, scores.null_auc, "x", color="0.5", label="null (intact vs intact)")
    for x, (auc, p) in enumerate(zip(scores.auc, scores.p_value)):
        ax_d.text(x + 0.08, auc, f"p = {p:.3g}", fontsize=6, va="center")
    ax_d.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_d.set_xticks(positions, scores.feature_set)
    ax_d.set_xlim(-0.5, len(scores) - 0.3)
    ax_d.set_ylim(min(0.3, np.nanmin(scores.null_auc) - 0.05), 1.05)
    ax_d.set_ylabel("ROC AUC")
    ax_d.legend(frameon=False, fontsize=7, loc="lower right")
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.set_title("d", loc="left", fontweight="bold", x=-0.15)
    ax_d.set_title("Intact vs broken experimental segments", fontsize=9)
    fig.savefig(prefix + "experimental_segments_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def experimental_analysis(features, segments, prefix):
    """Break detection and distinguishability of the experimental segments, with their
    figure."""
    features = features.fillna(FILL_VALUES)
    detection = experimental_break_detection(features)
    detection.to_csv(prefix + "experimental_break_detection.csv", index=False)
    scores = experimental_distinguishability(features)
    scores.to_csv(prefix + "experimental_distinguishability.csv", index=False)
    plot_experimental_figure(features, segments, detection, scores, prefix)
    print(detection.round(3).to_string(index=False))
    print(scores[["feature_set", "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc",
                  "p_value"]].round(3).to_string(index=False))


### Reference segment: overview, targets and hollow cylinder fit

def plot_overview(localisations_all, in_segment, fluorophore_positions, prefix):
    """The dataset with the reference segment, and the cropped segment aligned to x."""
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    overview, xedges, yedges = np.histogram2d(
        localisations_all.x, localisations_all.y, bins=int(localisations_all.x.max() / 40)
    )
    axs[0].imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
                  extent=[xedges[0], xedges[-1], yedges[-1], yedges[0]])
    axs[0].scatter(localisations_all.x[in_segment], localisations_all.y[in_segment], s=0.5, c="#F58B19")
    axs[0].set_title("Experimental localisations (segment in orange)")
    zoom, xedges, yedges = np.histogram2d(
        fluorophore_positions[:, 0], fluorophore_positions[:, 1],
        bins=[np.arange(-SEGMENT_LENGTH_NM / 2, SEGMENT_LENGTH_NM / 2 + 5, 5),
              np.arange(-SEGMENT_WIDTH_NM / 2, SEGMENT_WIDTH_NM / 2 + 5, 5)]
    )
    axs[1].imshow(zoom.T, cmap="grey", origin="lower", extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    axs[1].set_title("Cropped segment, aligned to x (nm)")
    fig.savefig(prefix + "segment_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_targets(particle_data, prefix):
    """The targets of the reference segment, from above and from the side."""
    primary = dict(
        probe_template="Antibody",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope=0,
        probe_model=PROBE_MODEL,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        probe_steric_hindrance=0,  # localisations are not excluded by probe size
        labelling_efficiency=1,
    )
    __, __, experiment = experiments.image_vsample(probe_list=[primary], **COMMON)
    experiment.particle.source = dict()
    experiment.particle.load_source(**copy.deepcopy(particle_data))
    experiment.set_virtualsample_params(sample_dimensions=SAMPLE_DIMENSIONS, update_mode=False)
    experiment.particle.generate_instance()
    experiment.build(modules=["coordinate_field", "imager"])
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    experiment.coordinate_field.show_field(view_init=[90, 0, 0], emitters_plotsize=1, axis_object=ax1)
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    experiment.coordinate_field.show_field(view_init=[0, 0, 0], emitters_plotsize=1, axis_object=ax2)
    fig.savefig(prefix + "filament.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def segment_cylinder_fit(points):
    """Anisotropic hollow cylinder fit of a segment, with the diameter of the isotropic fit
    (fit_hollow_cylinder, as in fig2_filaments.py) for comparison."""
    fit = fit_anisotropic_hollow_cylinder(points, reference_diameter=REFERENCE_DIAMETER,
                                          window=CYLINDER_WINDOW_NM)
    fit["isotropic_diameter"] = fit_hollow_cylinder(
        points, reference_diameter=REFERENCE_DIAMETER,
        background_radius_factor=SEGMENT_WIDTH_NM / 2 / REFERENCE_DIAMETER)["diameter"]
    return fit


def cylinder_fits(reference_points, segments=None):
    """Hollow cylinder fits of the reference segment ("ref") and of every intact experimental
    segment ("intact k", numbered as in the experimental figure), as (name, fit) pairs."""
    fits = [("ref", segment_cylinder_fit(reference_points))]
    if segments is not None and "points" in segments:
        intact = np.flatnonzero(segments["models"] == "intact")
        for k, i in enumerate(intact):
            fits.append((f"intact {k}",
                         segment_cylinder_fit(segments["points"][segments["point_segment"] == i])))
    return fits


def save_cylinder_fits(fits, prefix):
    """The hollow cylinder fits: their values (one row per segment), the fit of the reference
    segment and the fit of every segment (one row each)."""
    table = pd.DataFrame([dict(segment=name, **{key: fit[key] for key in
                               ["diameter", "sigma_u", "sigma_v", "background",
                                "isotropic_diameter", "reference_diameter", "n_points",
                                "success", "converged"]})
                          for name, fit in fits])
    table.to_csv(prefix + "hollow_cylinder_fit.csv", index=False)
    print(table.round(2).to_string(index=False))
    converged = table[table.converged]
    print(f"Hollow cylinder fits (n = {len(converged)} of {len(table)} converged): diameter "
          f"{converged.diameter.mean():.1f} ± {converged.diameter.std():.1f} nm, sigma_u "
          f"{converged.sigma_u.mean():.1f} ± {converged.sigma_u.std():.1f} nm, sigma_v "
          f"{converged.sigma_v.mean():.1f} ± {converged.sigma_v.std():.1f} nm (mean ± s.d.; "
          f"isotropic fit {converged.isotropic_diameter.mean():.1f} ± "
          f"{converged.isotropic_diameter.std():.1f} nm; reference {REFERENCE_DIAMETER} nm)")
    for name, shown in [("hollow_cylinder_fit.png", fits[:1]), ("hollow_cylinder_fits.png", fits)]:
        fig, axes = plt.subplots(len(shown), 3, squeeze=False, figsize=(18, 4.5 * len(shown)))
        for row, (segment, fit) in zip(axes, shown):
            plot_hollow_cylinder_fit(fit, row, color="#F58B19", title=segment)
        fig.tight_layout()
        fig.savefig(prefix + name, dpi=300 if len(shown) == 1 else 100, bbox_inches="tight")
        plt.close(fig)


def plot_cylinder_summary(fits, axs, color="#F58B19"):
    """The hollow cylinder fits of several segments on three axes: the cross-section and the
    profile across the filament of all their localisations, each segment in the frame of its
    own fitted axis, with the fit of each, and the fitted diameter of every segment with the
    mean +- s.d. of those that converged."""
    converged = np.array([fit["converged"] for __, fit in fits])
    diameters = np.array([fit["diameter"] for __, fit in fits])
    mean_diameter, std_diameter = diameters[converged].mean(), diameters[converged].std(ddof=1)
    shown = [fit for __, fit in fits if fit["converged"]]
    inside = [in_window(fit) for fit in shown]
    # cross-section, pooled
    limits = np.max([cross_section_limits(fit) for fit in shown], axis=0)
    plot_cross_section(axs[0], np.concatenate([fit["across_u"][k] for fit, k in zip(shown, inside)]),
                       np.concatenate([fit["across_v"][k] for fit, k in zip(shown, inside)]),
                       limits, mean_diameter / 2, color, f"mean fit: {mean_diameter:.1f} nm")
    axs[0].set_title(f"Cross-section, {len(shown)} segments")
    # profile across the filament, pooled, with the fit of every segment
    half_u = CYLINDER_WINDOW_NM[0]
    positions = np.linspace(-half_u, half_u, 300)
    axs[1].hist(np.concatenate([fit["across_u"][k] for fit, k in zip(shown, inside)]),
                bins=np.linspace(-half_u, half_u, 60), density=True, color=color, alpha=0.6)
    for k, fit in enumerate(shown):
        axs[1].plot(positions, cylinder_profile(fit, positions, "u"), c="0.3", lw=0.5, alpha=0.7,
                    label="fit of each segment" if k == 0 else None)
    for sign in (-1, 1):
        axs[1].axvline(sign * mean_diameter / 2, c=color)
        axs[1].axvline(sign * REFERENCE_DIAMETER / 2, ls="--", c="grey",
                       label=f"reference R = {REFERENCE_DIAMETER / 2:.1f} nm" if sign > 0 else None)
    axs[1].set_xlabel("u, across the filament (nm)")
    axs[1].set_ylabel("density")
    axs[1].legend(fontsize="x-small")
    axs[1].set_title("Profile across the filament (lateral)")
    # diameter of every segment
    jitter = np.random.default_rng(0).uniform(-0.15, 0.15, len(diameters))
    is_reference = np.array([name == "ref" for name, __ in fits])
    for selection, marker, size, colour, label in [
            (~is_reference, "o", 15, "0.4", "intact segments"),
            (is_reference, "D", 30, color, "reference segment")]:
        axs[2].scatter(jitter[selection & converged], diameters[selection & converged], s=size,
                       c=colour, marker=marker, label=label)
        if np.any(selection & ~converged):
            axs[2].scatter(jitter[selection & ~converged], diameters[selection & ~converged],
                           s=size, facecolors="none", edgecolors=colour, marker=marker,
                           label="not converged (excluded)" if not is_reference[selection][0] else None)
    axs[2].errorbar(0.4, mean_diameter, yerr=std_diameter, fmt="o", c="k", capsize=4,
                    label=f"{mean_diameter:.1f} ± {std_diameter:.1f} nm (mean ± s.d., "
                          f"n = {converged.sum()})")
    axs[2].axhline(REFERENCE_DIAMETER, ls="--", c="grey", label=f"reference: {REFERENCE_DIAMETER:.0f} nm")
    axs[2].set_xlim(-0.6, 0.8)
    axs[2].set_xticks([])
    axs[2].set_ylim(0, max(diameters.max(), REFERENCE_DIAMETER) * 1.6)
    axs[2].set_ylabel("Hollow cylinder diameter (nm)")
    axs[2].legend(fontsize="x-small", loc="upper left", frameon=False)
    axs[2].spines[["top", "right"]].set_visible(False)
    axs[2].set_title("Diameter per segment")


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
        examples = {}
        if os.path.exists(prefix + "examples.npz"):
            examples = dict(np.load(prefix + "examples.npz"))
        segments = None
        if os.path.exists(prefix + "experimental_features.csv") \
                and os.path.exists(prefix + "experimental_segments.npz"):
            segments = dict(np.load(prefix + "experimental_segments.npz"))
            experimental_analysis(pd.read_csv(prefix + "experimental_features.csv"),
                                  segments, prefix)
        fits = None
        if "fluorophore_positions" in examples:
            fits = cylinder_fits(examples["fluorophore_positions"], segments)
            save_cylinder_fits(fits, prefix)
    else:
        output_dir = os.path.join(str(output_path),
                                  ExperimentParametrisation().date_as_string + OUTPUT_NAME)
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, "")  # the outputs are named <output_dir>/<name>
        images_dir = prefix + "images" if SAVE_IMAGES else None

        localisations_all = load_smlm_table(SMLM_EXPERIMENTAL_PATH)
        fluorophore_positions, in_segment, reference_centre, reference_direction = crop_filament_segment(
            localisations_all, center=SEGMENT_CENTER, length=SEGMENT_LENGTH_NM, width=SEGMENT_WIDTH_NM
        )
        print(f"{fluorophore_positions.shape[0]} localisations used as fluorophore positions")
        plot_overview(localisations_all, in_segment, fluorophore_positions, prefix)
        template = experimental_template(fluorophore_positions)
        plot_targets(template, prefix)

        # experimental segments, without simulation
        experimental_features, segments = measure_experimental_segments(localisations_all)
        overview, xedges, yedges = np.histogram2d(
            localisations_all.x, localisations_all.y, bins=int(localisations_all.x.max() / 40))
        segments.update(overview=overview.astype(np.float32),
                        overview_extent=np.array([xedges[0], xedges[-1], yedges[-1], yedges[0]]),
                        reference_centre=reference_centre, reference_direction=reference_direction)
        experimental_features.to_csv(prefix + "experimental_features.csv", index=False)
        np.savez_compressed(prefix + "experimental_segments.npz", **segments)
        experimental_analysis(experimental_features, segments, prefix)
        fits = cylinder_fits(fluorophore_positions, segments)
        save_cylinder_fits(fits, prefix)

        plot_strip_examples(template, prefix)
        models = [("intact", 0, template)]
        n_targets = len(fluorophore_positions)
        for gap in GAPS_NM:
            broken = break_filament(template, gap)
            n = len(broken["targets"]["ATubulin_epitopes"]["coordinates"])
            print(f"gap {gap} nm: {n} of {n_targets} localisations")
            models.append(("broken", gap, broken))
        blocks = [(model, gap, le) for model, gap, __ in models for le in LABELLING_EFFICIENCIES]
        data = {(model, gap): particle_data for model, gap, particle_data in models}
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [(model, gap, le, int(seed), data[(model, gap)], images_dir)
                for (model, gap, le), seed in zip(blocks, seeds)]
        print(f"{len(jobs)} simulations, {N_REALISATIONS} realisations each, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                outputs = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="simulations"))
        else:
            outputs = [simulate(job) for job in tqdm(jobs, desc="simulations")]
        features = feature_table([table for table, __ in outputs])
        features.to_csv(prefix + "features.csv", index=False)
        examples = {key: images for __, found in outputs for key, images in found.items()}
        examples.update(model_epitopes(models))
        examples["fluorophore_positions"] = fluorophore_positions
        np.savez_compressed(prefix + "examples.npz", **examples)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict, tuple))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    plot_feature_distributions(features, prefix)
    features = features.fillna(FILL_VALUES)
    if examples:
        plot_examples(features, examples, prefix)
    scores = []
    for name, columns in FEATURE_SETS.items():
        print(f"features: {name}")
        table = distinguishability_from_features(
            features, "model", ("intact", "broken"), columns, by=CONDITIONS,
            controls=CONTROLS if name == "all" else None,
            replicate_column="replicate", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    scores = pd.concat(scores, ignore_index=True)
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    fig = plot_distinguishability_heatmaps(scores, x="labelling_efficiency", y="gap_nm",
                                           facet="feature_set", facet_order=list(FEATURE_SETS))
    n_le = scores.labelling_efficiency.nunique()
    gaps = sorted(scores.gap_nm.unique())
    panels = fig.axes[:len(FEATURE_SETS)]
    for ax in panels:
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Break length (nm)")
    fig.suptitle("Experimental microtubule: intact vs broken filament (ROC AUC)", x=0.02, ha="left")
    fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * len(gaps))
    fig.savefig(prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    plot_auc_vs_labelling(scores, prefix)
    detection = break_detection(features)
    detection.to_csv(prefix + "break_detection.csv", index=False)
    figure_values(detection, scores).to_csv(prefix + "filament_figure_values.csv", index=False)
    print(detection.round(3).to_string(index=False))
    if epitope_key("intact", 0) in examples:  # the targets of the models, from examples.npz
        plot_models(examples, prefix)
        plot_models_grid(examples, scores, prefix)
        if fits is not None:
            plot_figure(examples, detection, scores, fits, prefix)
    print(scores.pivot_table(index=["feature_set", "gap_nm"], columns="labelling_efficiency",
                             values="auc", sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
