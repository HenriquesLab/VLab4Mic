"""Fig. 2: can an under-labelled PCNA trimer be told apart from a fully labelled one, and at which resolution?

PCNA (4ZTD) is labelled at SER-186 of each chain (site-specific, no
linkage), i.e. three sites about 6 nm apart on a ring of ~3.3 nm radius
(as in figS4_pcna_sites.py), and imaged by SMLM, seen from above, always
in the same orientation (no random rotation), so that images can be
compared pixel by pixel. The fully labelled trimer (labelling efficiency
1, three labelled sites) is compared with trimers labelled at a lower
labelling efficiency, where each site carries a label with that
probability. The trimers are simulated with run_parameter_sweep over
labelling efficiency x localisation precision: every repetition is one
trimer with its own labelling and localisations. The localisations are
rendered with a Gaussian of RENDERING_KERNEL_NM at every precision, so
the image changes only through the spread of the localisations.

The reference image is one fully labelled trimer (REFERENCE_LE) imaged at
a localisation precision of REFERENCE_PRECISION_NM, simulated separately
(its own seed, so it is not one of the trimers it is compared with).

Every image is measured from its pixels only (the sweep gives images only):
- generic descriptors (vlab4mic.analysis.distinguishability): "shape"
  (radial profile, power spectrum, radius of gyration and principal axes,
  independent of brightness) and "brightness" (total intensity above
  background, which follows the number of labelled sites);
- PCNA descriptors, from the pixels within REGION_NM of the image centre
  above an Otsu threshold, taken as points of the trimer:
  - n_corners: occupied segments of the 3-fold ring (count_ring_corners,
    fixed radius MEAN_RADIUS_PCNA_NM; a segment needs at least
    MIN_SEGMENT_AREA_NM2 of bright pixels);
  - n_peaks: separate intensity peaks (count_resolved_sites), at most 3;
  - max_gap_deg: the largest gap in azimuth between ring points
    (angular_gaps); a missing site leaves a gap of ~240 degrees;
  - radius_nm: radius of a geometric circle fit to these points;
- triangle_score: how strongly the image repeats three times per turn
  around its centre of intensity (the magnitude of the 3-fold angular
  Fourier component within TRIANGLE_R_NM, over the total intensity there):
  0 for a round blob, 1 for three sharp points; independent of the
  orientation of the trimer and of brightness. It is measured on the image
  as SMLM images are usually rendered, each localisation a Gaussian as wide
  as the localisation precision: the image is smoothed with
  sqrt(precision^2 - RENDERING_KERNEL_NM^2) on top of its rendering kernel.
  This suppresses the speckle of a finite number of localisations without
  a free parameter. The triangle is "seen" when
  the score of a fully labelled trimer beats 1 - TRIANGLE_FALSE_ALARM of the
  scores of round blobs of the same size (N_NULL rings of the radius of the
  labelled sites, with as many localisations as the trimer at random angles,
  the same precision and rendering kernel, drawn with numpy);
- similarity to the reference image, within 2 * REGION_NM of the image
  centre (vlab4mic.analysis.metrics): "pearson_reference" (Pearson
  correlation of the pixels; 0 for an image without labels).

The feature table (one row per image) is scored per localisation precision
and labelling efficiency with distinguishability_from_features: out-of-fold
ROC AUC of a classifier (fully labelled vs reduced labelling) with a 95%
bootstrap interval, a null control (two random halves of the fully labelled
trimers, expected 0.5) and a permutation p-value, for each feature set. The
fully labelled trimers are simulated separately and shared by every
labelling efficiency: a labelling efficiency of 1 compares two independent
sets of fully labelled trimers, a negative control (expected 0.5).
auc_ideal = 1 - LE^3 / 2 is the AUC of a perfect count of labelled sites: a
trimer labelled at LE still has all three sites labelled with probability
LE^3, so no image feature can score higher.

Outputs, in the VLab4Mic output directory (<date>fig2_pcna_labelling_...):
pcna_pcna_figure.png, the main figure: (a) one fully labelled trimer imaged at
every localisation precision, (b) the percentage of fully labelled trimers whose
triangle is seen, with a 95% Wilson interval, at every precision, and (c)
the distinguishability (AUC) of labelling efficiencies 1 and FIGURE_LE at
every precision;
features.csv (one row per image), distinguishability.csv (one row per
localisation precision, labelling efficiency and feature set),
triangle_recovery.csv (one row per precision), settings.json,
reference.png (the reference image),
feature_distributions.png, distinguishability_heatmap.png,
auc_vs_labelling.png, and examples_LE<labelling efficiency>.png:
N_EXAMPLES trimers per localisation precision (DNA-PAINT and dSTORM named,
METHODS), those whose complete structure is detected (n_corners and
n_peaks both 3) first and framed in green, with the number of such
trimers, each with its Pearson correlation; their images and the
reference are kept in examples.npz, so --features plots them again (a
feature table from before the similarity descriptors cannot be scored).

The competing models are shown with one example trimer each (fully
labelled, and reduced at every labelling efficiency; a reduced example
misses at least one site when one of N_MODEL_CANDIDATES realisations
does), simulated apart from the scored trimers and imaged at every
precision: models.png, its virtual sample (labelled sites filled, the
unlabelled sites of the fully labelled trimer open) and its images, and
models_grid.png, the pair of images (fully labelled | reduced) compared in
each cell of the distinguishability map, laid out as the map and with the
AUC of all features. These are kept in examples.npz too.

Run:   python fig2_pcna_labelling_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_pcna_labelling_distinguishability.py --features <...features.csv>

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and
condition; e.g. 5 for a quick run), VLAB4MIC_N_PERMUTATIONS (200),
VLAB4MIC_N_WORKERS (1; sweeps and scoring in parallel processes),
VLAB4MIC_SAVE_IMAGES (0; 1 saves the simulated images as TIFF, one folder
per sweep, in <...>images: the scored trimers, the reference image and the
example of each competing model at every precision).
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage
from tqdm.auto import tqdm

from vlab4mic.analysis import particle_measures
from vlab4mic.analysis.metrics import pearson_correlation
from vlab4mic.analysis.distinguishability import (
    FEATURES as GENERIC_FUNCTIONS,
    distinguishability_from_features,
    feature_names,
    plot_distinguishability_heatmaps,
    same_field_of_view,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path
from vlab4mic.sweep_generator import run_parameter_sweep

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in (
    "0",
    "",
    "false",
    "False",
)
RANDOM_SEED = 24

STRUCTURE = "4ZTD"
MODALITY = "SMLM"
PIXEL_NM = 0.5
FOV_NM = 30  # common field of view of the generic descriptors

LABELLING_EFFICIENCIES = [0.3, 0.5, 0.8, 1.0]  # 1: negative control
# localisation precision (nm); 3 ~ DNA-PAINT, 6.5 ~ dSTORM (figS4_pcna_sites.py)
PRECISIONS_NM = [1, 2, 3, 4.5, 6.5]
N_LOCALISATIONS = 50  # per labelled site
RENDERING_KERNEL_NM = (
    1.0  # same for every precision: the localisations set the resolution
)

# reference image of the similarity descriptors: one fully labelled trimer
REFERENCE_LE = 1.0
REFERENCE_PRECISION_NM = 1

# measurement of the trimer
REGION_NM = 10  # the trimer is at the image centre: ignore pixels further out
N_SEGMENTS = 3  # PCNA trimer: three labelled sites
MEAN_RADIUS_PCNA_NM = 3.3  # SER-186 sites, ~6 nm apart
R_MIN_NM = 1.5  # ring points: between R_MIN_NM and R_MAX_NM of the centre
R_MAX_NM = 5
MIN_SEGMENT_AREA_NM2 = 2  # bright pixels a segment needs to count as detected
MIN_SEPARATION_NM = 3  # between separate peaks
# triangle score: annulus around the centre of intensity, and its test against round blobs;
# scored on the image rendered with a Gaussian as wide as the localisation precision
TRIANGLE_R_NM = (0.5, 8)
N_NULL = 1000  # round blobs per precision
TRIANGLE_FALSE_ALARM = 0.05

FIGURE_LE = (
    0.5  # labelling efficiency compared with full labelling in the main figure
)

# example images at every precision: trimers whose complete structure is
# detected (n_corners and n_peaks both 3) first, then others
METHODS = {3: "DNA-PAINT", 6.5: "dSTORM"}  # precision of each SMLM method (nm)
N_EXAMPLES = 8  # per precision and labelling efficiency
# example of each competing model: realisations to pick a reduced trimer from
N_MODEL_CANDIDATES = 20

GENERIC = ["shape", "brightness"]
SPECIFIC = ["n_corners", "n_peaks", "max_gap_deg", "radius_nm"]
TRIANGLE = ["triangle_score"]
# against the reference image (not "pearson": the sweep keeps that name for
# its own metric, against a reference image of its own)
SIMILARITY = ["pearson_reference"]
FEATURE_SETS = {
    "shape": "shape",
    "brightness": "brightness",
    "corners": ["n_corners"],
    "peaks": ["n_peaks"],
    "largest gap": ["max_gap_deg"],
    "specific": SPECIFIC,
    "triangle": TRIANGLE,
    "similarity": SIMILARITY,
    "all": None,
}  # all: filled in from the table
# each alone, with "all"
CONTROLS = {
    feature: feature for feature in GENERIC + SPECIFIC + TRIANGLE + SIMILARITY
}
# images without a measurable trimer (e.g. no labels): nothing detected, uncorrelated
FILL_VALUES = dict(
    n_corners=0,
    n_peaks=0,
    max_gap_deg=360,
    radius_nm=0,
    triangle_score=0,
    pearson_reference=0,
)
CONDITIONS = ["localisation_precision_nm", "labelling_efficiency"]

COMMON = dict(  # arguments of every sweep
    structures=[STRUCTURE],
    modalities=[MODALITY],
    probe_templates=["Linker"],
    probe_target_type=["Atom_residue"],
    probe_target_value=[
        {"atoms": ["CA"], "residues": ["SER"], "position": 186}
    ],
    probe_distance_to_epitope=[0],
    sample_dimensions=[[100, 100, 50]],
    # seen from above, all in the same orientation: compared pixel by pixel with the reference
    random_orientations=False,
    random_rotations=False,
    pixelsize_nm=[PIXEL_NM],
    psf_voxel_nm=[PIXEL_NM],
    sweep_repetitions=N_REALISATIONS,
    use_reference=False,
    default_metrics=[],
    save_analysis_results=False,
    analysis_plots=False,
    save_sweep_images=False,
    run_analysis=False,
    return_generator=True,
)
MODELS = {"full": [1.0], "reduced": LABELLING_EFFICIENCIES}


def measure_pcna(image, pixelsize_nm):
    """Corners, peaks, largest angular gap and radius of the trimer at the image centre."""
    missing = dict(
        n_corners=np.nan, n_peaks=np.nan, max_gap_deg=np.nan, radius_nm=np.nan
    )
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        image = image.sum(axis=0)  # frames, x, y
    image = same_field_of_view(image, pixelsize_nm, 2 * REGION_NM)
    rows_px, cols_px = np.indices(image.shape)
    xy = (
        np.c_[cols_px.ravel() + 0.5, rows_px.ravel() + 0.5] * pixelsize_nm
    )  # pixel centres, nm
    image_centre = np.array(image.shape[::-1]) * pixelsize_nm / 2
    region = np.linalg.norm(xy - image_centre, axis=1) <= REGION_NM
    image = np.where(region.reshape(image.shape), image, image.min())
    xy, values = xy[region], image.ravel()[region]
    if np.ptp(values) == 0:
        return missing
    points = xy[
        values > particle_measures.compute_threshold(values, method="otsu")
    ]
    if len(points) < 3:
        return missing
    _, radius = particle_measures.fit_circle_geometric(points)
    min_pixels = max(1, int(round(MIN_SEGMENT_AREA_NM2 / pixelsize_nm**2)))
    corners = particle_measures.count_ring_corners(
        points,
        MEAN_RADIUS_PCNA_NM,
        R_MIN_NM,
        R_MAX_NM,
        N_SEGMENTS,
        n_steps=360,
        min_per_section=min_pixels,
    )
    ring = points[corners["kept"]]
    max_gap = 360.0
    if len(ring) > 1:
        max_gap = particle_measures.angular_gaps(ring, corners["centre"]).max()
    n_peaks = particle_measures.count_resolved_sites(
        image, pixelsize_nm, MIN_SEPARATION_NM
    )
    return dict(
        n_corners=corners["n_corners"],
        n_peaks=min(n_peaks, N_SEGMENTS),
        max_gap_deg=max_gap,
        radius_nm=radius,
    )


def as_rendered_at_precision(image, pixelsize_nm, precision_nm):
    """The image centre (2 * REGION_NM) as if each localisation were rendered with a
    Gaussian as wide as the localisation precision: smoothed with what the rendering
    kernel lacks."""
    image = np.clip(
        same_field_of_view(image, pixelsize_nm, 2 * REGION_NM), 0, None
    )
    extra_nm = np.sqrt(max(precision_nm**2 - RENDERING_KERNEL_NM**2, 0))
    if extra_nm > 0:
        image = ndimage.gaussian_filter(
            image, extra_nm / pixelsize_nm, mode="constant"
        )
    return image


def triangle_score(image, pixelsize_nm, precision_nm):
    """Strength of the 3-fold pattern of the image, rendered at the localisation precision,
    around its centre of intensity: 0 for a round blob, 1 for three sharp points."""
    image = as_rendered_at_precision(image, pixelsize_nm, precision_nm)
    if image.sum() <= 0:
        return np.nan
    centre = (np.array(image.shape) - 1) / 2
    rows, cols = (
        np.indices(image.shape) - centre[:, None, None]
    ) * pixelsize_nm
    weights = image / image.sum()
    z = (cols - (weights * cols).sum()) + 1j * (rows - (weights * rows).sum())
    r = np.abs(z)
    ring = (r > TRIANGLE_R_NM[0]) & (r < TRIANGLE_R_NM[1])
    if image[ring].sum() <= 0:
        return np.nan
    return float(
        np.abs((image[ring] * (z[ring] / r[ring]) ** 3).sum())
        / image[ring].sum()
    )


def round_blob_scores(precision, radius_nm, rng, n=N_NULL):
    """Triangle scores of round blobs as large as the trimer: its localisations (a Poisson
    number per site, as the sweep draws them) at random angles on the ring of its labelled
    sites, with the same precision, rendering kernel and pixel size."""
    n_px = int(round(2 * REGION_NM / PIXEL_NM))
    edges = np.linspace(-REGION_NM, REGION_NM, n_px + 1)
    scores = np.empty(n)
    for k in range(n):
        n_loc = rng.poisson(N_LOCALISATIONS, N_SEGMENTS).sum()
        angles = rng.uniform(0, 2 * np.pi, n_loc)
        points = radius_nm * np.c_[np.cos(angles), np.sin(angles)]
        points += rng.normal(0, precision, (n_loc, 2))
        image = np.histogram2d(
            points[:, 1], points[:, 0], bins=[edges, edges]
        )[0]
        image = ndimage.gaussian_filter(image, RENDERING_KERNEL_NM / PIXEL_NM)
        scores[k] = triangle_score(image, PIXEL_NM, precision)
    return scores


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


def triangle_recovery(features, radius_nm):
    """Per precision: the share of fully labelled trimers whose triangle is seen, i.e. whose
    triangle score beats 1 - TRIANGLE_FALSE_ALARM of the round blobs."""
    rng = np.random.default_rng(RANDOM_SEED)
    # the fully labelled trimers appear once per labelling efficiency: keep one copy
    full = features[
        (features.model == "full")
        & np.isclose(features.labelling_efficiency, 1.0)
    ]
    rows = []
    for precision, data in full.groupby("localisation_precision_nm"):
        threshold = np.quantile(
            round_blob_scores(precision, radius_nm, rng),
            1 - TRIANGLE_FALSE_ALARM,
        )
        seen = int((data.triangle_score > threshold).sum())
        low, high = wilson_interval(seen, len(data))
        rows.append(
            dict(
                localisation_precision_nm=precision,
                n_trimers=len(data),
                n_seen=seen,
                fraction_seen=seen / len(data),
                fraction_low=low,
                fraction_high=high,
                threshold=threshold,
                median_score=data.triangle_score.median(),
            )
        )
    return pd.DataFrame(rows)


def site_radius(examples):
    """Radius of the ring of labelled sites, from the fully labelled example trimer."""
    sites = examples.get(model_key("full", 1.0) + "_sites")
    if sites is None or len(sites) < N_SEGMENTS:
        return MEAN_RADIUS_PCNA_NM
    return float(np.linalg.norm(sites - sites.mean(axis=0), axis=1).mean())


def generic_columns():
    return [(name, feature_names(name)) for name in GENERIC]


def measure_similarity(image, pixelsize_nm, reference):
    """Pearson correlation of the image centre against the reference image."""
    image = same_field_of_view(image, pixelsize_nm, 2 * REGION_NM)
    # a uniform image (no labels) has no correlation
    pearson = np.nan
    if np.ptp(image) > 0:
        pearson = pearson_correlation(
            reference_image=reference,
            reference_image_pixelsize_nm=PIXEL_NM,
            simulated_image=image,
            simulated_image_pixelsize_nm=pixelsize_nm,
        )
    return dict(pearson_reference=float(pearson))


def measure_image(image, pixelsize_nm, reference, precision_nm):
    """Generic, PCNA and similarity descriptors of one image, as {column: value}."""
    values = {}
    prepared = same_field_of_view(image, pixelsize_nm, FOV_NM)
    for name, columns in generic_columns():
        values.update(
            zip(columns, np.atleast_1d(GENERIC_FUNCTIONS[name](prepared)))
        )
    values.update(measure_pcna(image, pixelsize_nm))
    values["triangle_score"] = triangle_score(
        image, pixelsize_nm, precision_nm
    )
    values.update(measure_similarity(image, pixelsize_nm, reference))
    return values


def pcna_metrics(reference, precision_nm):
    """measure_image as custom metrics of a sweep, one per column; each image is measured once."""
    last = {"image": None, "values": None}

    def measure(image, pixelsize_nm):
        # the sweep calls every metric in turn with the same image object
        if last["image"] is not image:
            last["image"] = image
            last["values"] = measure_image(
                image, float(pixelsize_nm), reference, precision_nm
            )
        return last["values"]

    columns = (
        [c for _, names in generic_columns() for c in names]
        + SPECIFIC
        + TRIANGLE
        + SIMILARITY
    )
    metrics = []
    for name in columns:

        def metric(
            simulated_image=None,
            simulated_image_pixelsize_nm=None,
            _name=name,
            **unused,
        ):
            return float(
                measure(simulated_image, simulated_image_pixelsize_nm)[_name]
            )

        metric.__name__ = name
        metrics.append(metric)
    return metrics


def complete_detection(table):
    """Rows where corners and peaks both find all three sites of the trimer."""
    return (table.n_corners == N_SEGMENTS) & (table.n_peaks == N_SEGMENTS)


def example_key(precision, le):
    return f"precision_{precision:g}_labelling_efficiency_{le:g}"


def precision_label(precision):
    method = f" ({METHODS[precision]})" if precision in METHODS else ""
    return f"{precision:g} nm{method}"


def scale_bar(ax, colour="w"):
    """5 nm, at the lower left of an image of 2 * REGION_NM."""
    ax.plot(
        [-REGION_NM + 2, -REGION_NM + 7],
        [-REGION_NM + 2] * 2,
        "-",
        color=colour,
        lw=2,
    )


def set_emitters(sweep, precision):
    # the sweep adds the modality from these parameters before imaging; the
    # SMLM image renders the localisations, so its resolution is the
    # localisation precision (not the PSF width of lateral_resolution_nm)
    emitters = sweep.experiment.local_modalities_parameters[MODALITY][
        "emitters"
    ]
    emitters.update(
        lateral_precision=precision,
        axial_precision=precision,
        nlocalisations=N_LOCALISATIONS,
        rendering_kernel_nm=RENDERING_KERNEL_NM,
    )


def save_images(sweep, images_dir, name):
    """The images of the sweep as TIFF, in images_dir/name (nothing if images_dir is None)."""
    if images_dir is None:
        return
    output_directory = os.path.join(images_dir, name)
    os.makedirs(output_directory, exist_ok=True)
    sweep.save_images(output_directory=output_directory)


def simulate_reference(seed, images_dir=None):
    """One fully labelled trimer at REFERENCE_PRECISION_NM (field of view 2 * REGION_NM)."""
    sweep = run_parameter_sweep(
        labelling_efficiency=[REFERENCE_LE],
        random_seed=seed,
        **dict(COMMON, sweep_repetitions=1),
    )
    set_emitters(sweep, REFERENCE_PRECISION_NM)
    sweep.generate_acquisitions()
    save_images(
        sweep,
        images_dir,
        f"reference_precision_{REFERENCE_PRECISION_NM}"
        f"_labelling_efficiency_{REFERENCE_LE}",
    )
    image = next(iter(sweep.acquisition_outputs.values()))[0]
    return same_field_of_view(image, PIXEL_NM, 2 * REGION_NM)


def competing_models():
    """(model, labelling efficiency) of every model compared in the distinguishability map."""
    return [("full", 1.0)] + [("reduced", le) for le in LABELLING_EFFICIENCIES]


def model_key(model, le):
    return f"model_{model}_labelling_efficiency_{le:g}"


def labelled_sites(virtual_sample):
    """Labelled sites of a trimer (nm, from its centre), as (horizontal, vertical) in the
    images: the rows of an image run along x and its columns along y."""
    emitters = [
        np.reshape(e, (-1, 3))
        for e in virtual_sample["field_emitters"].values()
    ]
    emitters = np.vstack(emitters) if emitters else np.empty((0, 3))
    centre = np.asarray(virtual_sample["particle_positions"]).reshape(-1, 3)[0]
    return (emitters[:, :2] - centre[:2])[:, ::-1]


def simulate_model_examples(seed, images_dir=None):
    """One example trimer per competing model: its labelled sites, and its images at every
    precision (field of view 2 * REGION_NM). A reduced example misses some but not all sites
    when one of N_MODEL_CANDIDATES realisations does."""
    models = competing_models()
    seeds = np.random.SeedSequence(seed).generate_state(len(models))
    examples = {}
    for (model, le), model_seed in zip(models, seeds):
        sweep = run_parameter_sweep(
            labelling_efficiency=[le],
            random_seed=int(model_seed),
            **dict(COMMON, sweep_repetitions=N_MODEL_CANDIDATES),
        )
        sweep.generate_virtual_samples()
        combination = next(iter(sweep.virtual_samples))
        candidates = sweep.virtual_samples[combination]
        counts = np.array([len(labelled_sites(v)) for v in candidates])
        incomplete = np.flatnonzero((counts > 0) & (counts < N_SEGMENTS))
        pick = (
            int(incomplete[0]) if model == "reduced" and len(incomplete) else 0
        )
        # image only the example, once per precision
        sweep.virtual_samples[combination] = [candidates[pick]]
        images = []
        for precision in PRECISIONS_NM:
            set_emitters(sweep, precision)
            sweep.generate_acquisitions()
            save_images(
                sweep,
                images_dir,
                f"example_{model}_precision_{precision}_labelling_efficiency_{le}",
            )
            image = next(iter(sweep.acquisition_outputs.values()))[0]
            images.append(same_field_of_view(image, PIXEL_NM, 2 * REGION_NM))
        key = model_key(model, le)
        examples[key + "_sites"] = labelled_sites(candidates[pick])
        examples[key + "_images"] = np.array(images)
    return examples


def simulate(job):
    """One model, localisation precision and labelling efficiency: its feature table, and
    N_EXAMPLES images (field of view 2 * REGION_NM), complete detections first.
    """
    model, precision, le, seed, images_dir, reference = job
    sweep = run_parameter_sweep(
        labelling_efficiency=[le],
        custom_metrics=pcna_metrics(reference, precision),
        random_seed=seed,
        **COMMON,
    )
    set_emitters(sweep, precision)
    sweep.run_analysis(save=False, plots=False)
    save_images(
        sweep,
        images_dir,
        f"{model}_precision_{precision}_labelling_efficiency_{le}",
    )
    table = sweep.analysis["dataframes"].copy()
    table["model"] = model
    table["localisation_precision_nm"] = precision
    table["labelling_efficiency"] = le
    examples = {}
    if model == "reduced":
        # one parameter combination per sweep; Replica indexes its images,
        # and the metrics measured the first frame of each
        images = next(iter(sweep.acquisition_outputs.values()))
        replica = table["Replica"].astype(int)
        complete = complete_detection(table)
        replicas = pd.concat(
            [replica[complete], replica[~complete]]
        ).to_numpy()[:N_EXAMPLES]
        n_px = int(round(2 * REGION_NM / PIXEL_NM))
        key = example_key(precision, le)
        examples[key] = np.array(
            [
                same_field_of_view(images[r][0], PIXEL_NM, 2 * REGION_NM)
                for r in replicas
            ]
        ).reshape(-1, n_px, n_px)
        examples[key + "_replicas"] = replicas
    return table, examples


def feature_table(tables):
    """Fully labelled and reduced rows for every labelling efficiency."""
    features = pd.concat(tables, ignore_index=True)
    full = features[features.model == "full"]
    reduced = features[features.model == "reduced"]
    # the fully labelled trimers are shared by every labelling efficiency
    return pd.concat(
        [reduced]
        + [
            full.assign(labelling_efficiency=le)
            for le in LABELLING_EFFICIENCIES
        ],
        ignore_index=True,
    )


def plot_feature_distributions(features, prefix):
    """Mean +- s.d. of each PCNA descriptor and brightness against labelling efficiency."""
    reduced = features[features.model == "reduced"]
    precisions = sorted(reduced.localisation_precision_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(precisions)))
    columns = SPECIFIC + TRIANGLE + ["brightness"] + SIMILARITY
    fig, axes = plt.subplots(1, len(columns), figsize=[3.5 * len(columns), 3])
    for ax, feature in zip(axes, columns):
        for precision, colour in zip(precisions, colours):
            summary = (
                reduced[
                    np.isclose(reduced.localisation_precision_nm, precision)
                ]
                .groupby("labelling_efficiency")[feature]
                .agg(["mean", "std"])
            )
            ax.errorbar(
                summary.index,
                summary["mean"],
                yerr=summary["std"],
                fmt="o-",
                color=colour,
                capsize=2,
                markersize=3,
                label=f"{precision:g}",
            )
        ax.set_xlabel("Labelling efficiency")
        ax.set_ylabel(feature.replace("_", " "))
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].legend(title="Precision (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(
        prefix + "feature_distributions.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)


def plot_auc_vs_labelling(
    scores,
    prefix,
    feature_sets=("shape", "brightness", "specific", "similarity", "all"),
):
    """AUC (95% interval) against labelling efficiency, one line per precision."""
    precisions = sorted(scores.localisation_precision_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(precisions)))
    fig, axes = plt.subplots(
        1, len(feature_sets), figsize=[3.5 * len(feature_sets), 3], sharey=True
    )
    for ax, name in zip(axes, feature_sets):
        data = scores[scores.feature_set == name]
        for precision, colour in zip(precisions, colours):
            sub = data[np.isclose(data.localisation_precision_nm, precision)]
            sub = sub.sort_values("labelling_efficiency")
            ax.errorbar(
                sub.labelling_efficiency,
                sub.auc,
                yerr=error_bars(sub.auc, sub.auc_low, sub.auc_high),
                fmt="o-",
                color=colour,
                capsize=2,
                markersize=3,
                label=f"{precision:g}",
            )
        le = np.linspace(min(LABELLING_EFFICIENCIES), 1, 50)
        ax.plot(le, 1 - le**3 / 2, "k--", lw=1, label="perfect site count")
        ax.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
        ax.set_ylim(0.3, 1.05)
        ax.set_xlabel("Labelling efficiency")
        ax.set_title(name, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("ROC AUC (full vs reduced labelling)")
    axes[-1].legend(title="Precision (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "auc_vs_labelling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_examples(features, examples, prefix):
    """Per labelling efficiency: example trimers, one row per localisation precision, those
    with a complete detection first (green frame), each with its corners and peaks.
    """
    reduced = features[features.model == "reduced"]
    extent = [-REGION_NM, REGION_NM, -REGION_NM, REGION_NM]
    for le in LABELLING_EFFICIENCIES:
        fig, axes = plt.subplots(
            len(PRECISIONS_NM),
            N_EXAMPLES,
            figsize=[1.6 * N_EXAMPLES, 2.1 * len(PRECISIONS_NM)],
            squeeze=False,
        )
        for row, precision in zip(axes, PRECISIONS_NM):
            key = example_key(precision, le)
            images = examples.get(key, [])
            replicas = examples.get(key + "_replicas", [])
            data = reduced[
                np.isclose(reduced.localisation_precision_nm, precision)
                & np.isclose(reduced.labelling_efficiency, le)
            ]
            measured = data.set_index(data.Replica.astype(int))
            for k, ax in enumerate(row):
                ax.set_xticks([])
                ax.set_yticks([])
                if k >= len(images):
                    ax.axis("off")
                    continue
                ax.imshow(
                    images[k], cmap="magma", extent=extent, origin="lower"
                )
                if k < len(replicas) and replicas[k] in measured.index:
                    found = measured.loc[replicas[k]]
                    complete = (
                        found.n_corners == N_SEGMENTS
                        and found.n_peaks == N_SEGMENTS
                    )
                    ax.set_title(
                        f"{found.n_corners:.0f} corners, {found.n_peaks:.0f} peaks\n"
                        f"r = {found.pearson_reference:.2f}",
                        fontsize=7,
                        color="C2" if complete else "0.4",
                    )
                    if complete:
                        ax.spines[:].set_color("C2")
                        ax.spines[:].set_linewidth(2.5)
            if len(images):
                scale_bar(row[0])
            row[0].axis("on")
            row[0].set_ylabel(
                f"{precision_label(precision)}\n"
                f"{int(complete_detection(data).sum())} / {len(data)} complete",
                fontsize=9,
            )
        fig.suptitle(
            f"PCNA, labelling efficiency {le:g}: localisation precision per row; "
            f"green, complete structure detected ({N_SEGMENTS} corners and "
            f"{N_SEGMENTS} peaks); scale bar 5 nm",
            x=0.02,
            ha="left",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(
            prefix + f"examples_LE{le:g}.png", dpi=300, bbox_inches="tight"
        )
        plt.close(fig)


def plot_reference(reference, prefix):
    """The reference image of the similarity descriptors."""
    extent = [-REGION_NM, REGION_NM, -REGION_NM, REGION_NM]
    fig, ax = plt.subplots(figsize=[2.5, 2.7])
    ax.imshow(reference, cmap="magma", extent=extent, origin="lower")
    scale_bar(ax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"Reference: LE {REFERENCE_LE:g}, {REFERENCE_PRECISION_NM:g} nm",
        fontsize=9,
    )
    fig.savefig(prefix + "reference.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def model_title(model, le, n_sites):
    if model == "full":
        return f"Fully labelled\n{n_sites}/{N_SEGMENTS} sites labelled"
    control = " (control)" if le == 1 else ""
    return (
        f"Reduced, LE {le:g}{control}\n{n_sites}/{N_SEGMENTS} sites labelled; "
        f"P(all {N_SEGMENTS}) = {le**N_SEGMENTS:.2f}"
    )


def plot_models(examples, prefix):
    """The competing models, one example trimer each (column): its virtual sample (labelled
    sites filled, the unlabelled sites of the fully labelled trimer open) and its image at
    every localisation precision (rows)."""
    models = competing_models()
    limits = [-REGION_NM, REGION_NM]
    all_sites = examples[model_key("full", 1.0) + "_sites"]
    n_rows = 1 + len(PRECISIONS_NM)
    fig, axes = plt.subplots(
        n_rows,
        len(models),
        figsize=[2.0 * len(models), 1.9 * n_rows],
        squeeze=False,
    )
    for column, (model, le) in zip(axes.T, models):
        key = model_key(model, le)
        sites = examples[key + "_sites"]
        column[0].scatter(
            *all_sites.T,
            s=50,
            facecolors="none",
            edgecolors="0.6",
            linewidths=1,
        )
        if len(sites):
            column[0].scatter(*sites.T, s=50, color="C3")
        column[0].set(xlim=limits, ylim=limits, aspect="equal")
        column[0].set_title(model_title(model, le, len(sites)), fontsize=7)
        for ax, image in zip(column[1:], examples[key + "_images"]):
            ax.imshow(
                image, cmap="magma", extent=limits + limits, origin="lower"
            )
        for ax in column:
            ax.set_xticks([])
            ax.set_yticks([])
    scale_bar(axes[0, 0], colour="k")
    axes[0, 0].set_ylabel("Virtual sample", fontsize=9)
    for ax, precision in zip(axes[1:, 0], PRECISIONS_NM):
        scale_bar(ax)
        ax.set_ylabel(precision_label(precision), fontsize=9)
    fig.suptitle(
        "PCNA, competing models: one example trimer each, as a virtual sample "
        "(filled, labelled sites; open, unlabelled) and imaged at every localisation "
        "precision; scale bar 5 nm",
        x=0.02,
        ha="left",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(prefix + "models.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_models_grid(examples, scores, prefix):
    """The pair of images compared in every cell of the distinguishability map (fully
    labelled | reduced), laid out as the map, with the AUC of all features."""
    precisions = sorted(
        PRECISIONS_NM, reverse=True
    )  # as the map: finest precision at the bottom
    n_le = len(LABELLING_EFFICIENCIES)
    ratios = ([1, 1, 0.3] * n_le)[
        :-1
    ]  # a pair of panels per labelling efficiency, then a gap
    fig, axes = plt.subplots(
        len(precisions),
        len(ratios),
        squeeze=False,
        figsize=[0.9 * sum(ratios) + 0.8, 0.95 * len(precisions) + 0.8],
        gridspec_kw=dict(width_ratios=ratios, wspace=0.05, hspace=0.1),
    )
    limits = [-REGION_NM, REGION_NM]
    full_images = examples[model_key("full", 1.0) + "_images"]
    auc = scores[scores.feature_set == "all"]
    for i, (row, precision) in enumerate(zip(axes, precisions)):
        p = PRECISIONS_NM.index(precision)
        for k, le in enumerate(LABELLING_EFFICIENCIES):
            full, reduced = row[3 * k], row[3 * k + 1]
            full.imshow(
                full_images[p],
                cmap="magma",
                extent=limits + limits,
                origin="lower",
            )
            reduced.imshow(
                examples[model_key("reduced", le) + "_images"][p],
                cmap="magma",
                extent=limits + limits,
                origin="lower",
            )
            cell = auc[
                np.isclose(auc.localisation_precision_nm, precision)
                & np.isclose(auc.labelling_efficiency, le)
            ]
            if len(cell):
                reduced.text(
                    0.95,
                    0.93,
                    f"AUC {cell.auc.iloc[0]:.2f}",
                    color="w",
                    fontsize=6,
                    ha="right",
                    va="top",
                    transform=reduced.transAxes,
                )
            if i == 0:
                full.set_title("full", fontsize=7)
                reduced.set_title(f"LE {le:g}", fontsize=7)
        for ax in row:
            ax.set_xticks([])
            ax.set_yticks([])
        for ax in row[2::3]:
            ax.axis("off")
        row[0].set_ylabel(
            precision_label(precision).replace(" (", "\n("), fontsize=7
        )
    scale_bar(axes[-1, 0])
    fig.supylabel("Localisation precision (nm)", fontsize=9)
    fig.supxlabel("Labelling efficiency of the reduced model", fontsize=9)
    fig.suptitle(
        "PCNA: trimers compared in each cell of the distinguishability map "
        "(fully labelled | reduced; AUC, all features); scale bar 5 nm",
        x=0.02,
        ha="left",
        fontsize=10,
    )
    fig.savefig(prefix + "models_grid.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_figure(examples, recovery, scores, prefix):
    """The main figure: (a) one fully labelled trimer at every localisation precision,
    (b) the share of fully labelled trimers whose triangle is seen, (c) the AUC of full
    labelling against FIGURE_LE at every precision, from the triangle score. Precisions are
    labelled in nm only, without the imaging methods."""
    precisions = sorted(PRECISIONS_NM)
    fig = plt.figure(figsize=[11, 6.8])
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.5], hspace=0.35)
    top = grid[0].subgridspec(1, len(precisions), wspace=0.05)
    limits = [-REGION_NM, REGION_NM]
    images = examples[model_key("full", 1.0) + "_images"]
    for k, precision in enumerate(precisions):
        ax = fig.add_subplot(top[k])
        ax.imshow(
            images[PRECISIONS_NM.index(precision)],
            cmap="magma",
            extent=limits + limits,
            origin="lower",
        )
        ax.set_title(f"{precision:g} nm", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        if k == 0:
            scale_bar(ax)
            ax.set_title("a", loc="left", fontweight="bold", x=-0.12)
    bottom = grid[1].subgridspec(1, 2, wspace=0.3)
    ticks = [f"{p:g}" for p in precisions]

    ax = fig.add_subplot(bottom[0])
    recovery = recovery.sort_values("localisation_precision_nm")
    percent = 100 * recovery.fraction_seen
    ax.errorbar(
        recovery.localisation_precision_nm,
        percent,
        yerr=error_bars(
            percent,
            100 * recovery.fraction_low,
            100 * recovery.fraction_high,
        ),
        fmt="o-",
        color="C0",
        capsize=3,
    )
    ax.axhline(
        100 * TRIANGLE_FALSE_ALARM, color="0.5", linestyle=":", linewidth=1
    )
    ax.text(
        precisions[-1],
        100 * TRIANGLE_FALSE_ALARM + 2,
        "chance",
        color="0.5",
        fontsize=7,
        ha="right",
        va="bottom",
    )
    ax.set_ylim(-3, 103)
    ax.set_ylabel("Triangle seen (% of trimers)")
    ax.set_title("b", loc="left", fontweight="bold", x=-0.18)
    ax.set_title(
        f"Fully labelled trimers (n = {recovery.n_trimers.max()} per precision)",
        fontsize=9,
    )

    ax_c = fig.add_subplot(bottom[1])
    data = scores[
        (scores.feature_set == "triangle")
        & np.isclose(scores.labelling_efficiency, FIGURE_LE)
    ].sort_values("localisation_precision_nm")
    ax_c.errorbar(
        data.localisation_precision_nm,
        data.auc,
        yerr=error_bars(data.auc, data.auc_low, data.auc_high),
        fmt="o-",
        color="C0",
        capsize=3,
        label="triangle score",
    )
    ax_c.plot(
        data.localisation_precision_nm,
        data.null_auc,
        "x",
        color="k",
        label="null (full vs full)",
    )
    ax_c.axhline(
        1 - FIGURE_LE**N_SEGMENTS / 2,
        color="k",
        linestyle="--",
        linewidth=1,
        label="perfect site count",
    )
    ax_c.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_c.set_ylim(min(0.3, np.nanmin(data.null_auc) - 0.05), 1.05)
    ax_c.set_ylabel("ROC AUC")
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.18)
    ax_c.set_title(f"Labelling efficiency 1 vs {FIGURE_LE:g}", fontsize=9)
    ax_c.legend(frameon=False, fontsize=7, loc="lower left")
    for axis in (ax, ax_c):
        axis.set_xticks(precisions, ticks, fontsize=8)
        axis.set_xlabel("Localisation precision (nm)")
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(prefix + "pcna_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
        examples = {}
        if os.path.exists(prefix + "examples.npz"):
            examples = dict(np.load(prefix + "examples.npz"))
    else:
        prefix = os.path.join(
            str(output_path),
            ExperimentParametrisation().date_as_string
            + "fig2_pcna_labelling_",
        )
        images_dir = prefix + "images" if SAVE_IMAGES else None
        blocks = [
            (model, precision, le)
            for model in MODELS
            for precision in PRECISIONS_NM
            for le in MODELS[model]
        ]
        # two more seeds, last, for the reference and the examples of the models: the
        # seeds of the sweeps are unchanged
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(
            len(blocks) + 2
        )
        print("reference image")
        reference = simulate_reference(int(seeds[len(blocks)]), images_dir)
        print("examples of the competing models")
        model_examples = simulate_model_examples(
            int(seeds[len(blocks) + 1]), images_dir
        )
        jobs = [
            block + (int(seed), images_dir, reference)
            for block, seed in zip(blocks, seeds)
        ]
        print(
            f"{len(jobs)} sweeps, {N_REALISATIONS} realisations per model and condition, "
            f"{N_WORKERS} worker(s)"
        )
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                outputs = list(
                    tqdm(
                        pool.map(simulate, jobs),
                        total=len(jobs),
                        desc="sweeps",
                    )
                )
        else:
            outputs = [simulate(job) for job in tqdm(jobs, desc="sweeps")]
        tables = [table for table, _ in outputs]
        examples = {
            key: images
            for _, found in outputs
            for key, images in found.items()
        }
        features = feature_table(tables)
        features.to_csv(prefix + "features.csv", index=False)
        examples["reference"] = reference
        examples.update(model_examples)
        np.savez_compressed(prefix + "examples.npz", **examples)
        settings = {
            k: v
            for k, v in globals().items()
            if k.isupper() and isinstance(v, (int, float, str, list, dict))
        }
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    if "reference" in examples:
        plot_reference(examples.pop("reference"), prefix)
    plot_feature_distributions(features, prefix)
    features = features.fillna(FILL_VALUES)
    plot_examples(features, examples, prefix)
    generic = [c for _, names in generic_columns() for c in names]
    feature_sets = dict(
        FEATURE_SETS, all=generic + SPECIFIC + TRIANGLE + SIMILARITY
    )
    scores = []
    for name, columns in feature_sets.items():
        print(f"features: {name}")
        table = distinguishability_from_features(
            features,
            "model",
            ("full", "reduced"),
            columns,
            by=CONDITIONS,
            controls=CONTROLS if name == "all" else None,
            replicate_column="Replica",
            n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED,
            n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    scores = pd.concat(scores, ignore_index=True)
    # AUC of a perfect count of labelled sites: the ceiling of any feature
    scores["auc_ideal"] = 1 - scores.labelling_efficiency**N_SEGMENTS / 2
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    fig = plot_distinguishability_heatmaps(
        scores,
        x="labelling_efficiency",
        y="localisation_precision_nm",
        facet="feature_set",
        facet_order=list(feature_sets),
    )
    n_le = scores.labelling_efficiency.nunique()
    n_precision = scores.localisation_precision_nm.nunique()
    panels = fig.axes[: len(feature_sets)]
    for ax in panels:
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Localisation precision (nm)")
    fig.suptitle(
        "PCNA: fully labelled vs reduced labelling (ROC AUC)",
        x=0.02,
        ha="left",
    )
    fig.set_size_inches(
        1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * n_precision
    )
    fig.savefig(
        prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    plot_auc_vs_labelling(scores, prefix)
    recovery = triangle_recovery(features, site_radius(examples))
    recovery.to_csv(prefix + "triangle_recovery.csv", index=False)
    print(recovery.round(3).to_string(index=False))
    if (
        model_key("full", 1.0) + "_images" in examples
    ):  # not in examples from older runs
        plot_models(examples, prefix)
        plot_models_grid(examples, scores, prefix)
        plot_figure(examples, recovery, scores, prefix)
    print(
        scores.pivot_table(
            index=["feature_set", "localisation_precision_nm"],
            columns="labelling_efficiency",
            values="auc",
            sort=False,
        )
        .round(2)
        .to_string()
    )


if __name__ == "__main__":
    main()
