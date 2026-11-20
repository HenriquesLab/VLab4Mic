"""Fig. 2: can an under-labelled PCNA trimer be told apart from a fully labelled one, and at which resolution?

PCNA (4ZTD) is labelled at SER-186 of each chain (site-specific, no
linkage), i.e. three sites about 6 nm apart on a ring of ~3.3 nm radius
(as in figS4_pcna_sites.py), and imaged by SMLM, seen from above. The
reference is the fully labelled trimer (labelling efficiency 1, three
labelled sites); it is compared with trimers labelled at a lower
labelling efficiency, where each site carries a label with that
probability. The trimers are simulated with run_parameter_sweep over
labelling efficiency x localisation precision: every repetition is one
trimer with its own labelling, in-plane rotation and localisations. The
localisations are rendered with a Gaussian of RENDERING_KERNEL_NM at every
precision, so the image changes only through the spread of the
localisations.

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
  - radius_nm: radius of a geometric circle fit to these points.

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
features.csv (one row per image), distinguishability.csv (one row per
localisation precision, labelling efficiency and feature set),
settings.json, feature_distributions.png, distinguishability_heatmap.png
and auc_vs_labelling.png.

Run:   python fig2_pcna_labelling_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_pcna_labelling_distinguishability.py --features <...features.csv>

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and
condition; e.g. 5 for a quick run), VLAB4MIC_N_PERMUTATIONS (200),
VLAB4MIC_N_WORKERS (1; sweeps and scoring in parallel processes),
VLAB4MIC_SAVE_IMAGES (0; 1 saves the simulated images as TIFF, one folder
per sweep).
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from vlab4mic.analysis import particle_measures
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
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 61

STRUCTURE = "4ZTD"
MODALITY = "SMLM"
PIXEL_NM = 0.5
FOV_NM = 30  # common field of view of the generic descriptors

LABELLING_EFFICIENCIES = [0.3, 0.5, 0.8, 1.0]  # 1: negative control
# localisation precision (nm); 3 ~ DNA-PAINT, 6.5 ~ dSTORM (figS4_pcna_sites.py)
PRECISIONS_NM = [1, 2, 3, 4.5, 6.5]
N_LOCALISATIONS = 10  # per labelled site
RENDERING_KERNEL_NM = 1.0  # same for every precision: the localisations set the resolution

# measurement of the trimer
REGION_NM = 10  # the trimer is at the image centre: ignore pixels further out
N_SEGMENTS = 3  # PCNA trimer: three labelled sites
MEAN_RADIUS_PCNA_NM = 3.3  # SER-186 sites, ~6 nm apart
R_MIN_NM = 1.5  # ring points: between R_MIN_NM and R_MAX_NM of the centre
R_MAX_NM = 5
MIN_SEGMENT_AREA_NM2 = 2  # bright pixels a segment needs to count as detected
MIN_SEPARATION_NM = 3  # between separate peaks

GENERIC = ["shape", "brightness"]
SPECIFIC = ["n_corners", "n_peaks", "max_gap_deg", "radius_nm"]
FEATURE_SETS = {"shape": "shape", "brightness": "brightness",
                "corners": ["n_corners"], "peaks": ["n_peaks"], "largest gap": ["max_gap_deg"],
                "specific": SPECIFIC, "all": None}  # all: filled in from the table
CONTROLS = {feature: feature for feature in GENERIC + SPECIFIC}  # each alone, with "all"
# images without a measurable trimer (e.g. no labels): nothing detected
FILL_VALUES = dict(n_corners=0, n_peaks=0, max_gap_deg=360, radius_nm=0)
CONDITIONS = ["localisation_precision_nm", "labelling_efficiency"]

COMMON = dict(  # arguments of every sweep
    structures=[STRUCTURE], modalities=[MODALITY],
    probe_templates=["Linker"], probe_target_type=["Atom_residue"],
    probe_target_value=[{"atoms": ["CA"], "residues": ["SER"], "position": 186}],
    probe_distance_to_epitope=[0],
    sample_dimensions=[[100, 100, 50]],
    random_orientations=False, random_rotations=True,  # seen from above, rotated in plane
    pixelsize_nm=[PIXEL_NM], psf_voxel_nm=[PIXEL_NM],
    sweep_repetitions=N_REALISATIONS, use_reference=False, default_metrics=[],
    save_analysis_results=False, analysis_plots=False, save_sweep_images=False,
    run_analysis=False, return_generator=True,
)
MODELS = {"full": [1.0], "reduced": LABELLING_EFFICIENCIES}


def measure_pcna(image, pixelsize_nm):
    """Corners, peaks, largest angular gap and radius of the trimer at the image centre."""
    missing = dict(n_corners=np.nan, n_peaks=np.nan, max_gap_deg=np.nan, radius_nm=np.nan)
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        image = image.sum(axis=0)  # frames, x, y
    image = same_field_of_view(image, pixelsize_nm, 2 * REGION_NM)
    rows_px, cols_px = np.indices(image.shape)
    xy = np.c_[cols_px.ravel() + 0.5, rows_px.ravel() + 0.5] * pixelsize_nm  # pixel centres, nm
    image_centre = np.array(image.shape[::-1]) * pixelsize_nm / 2
    region = np.linalg.norm(xy - image_centre, axis=1) <= REGION_NM
    image = np.where(region.reshape(image.shape), image, image.min())
    xy, values = xy[region], image.ravel()[region]
    if np.ptp(values) == 0:
        return missing
    points = xy[values > particle_measures.compute_threshold(values, method="otsu")]
    if len(points) < 3:
        return missing
    _, radius = particle_measures.fit_circle_geometric(points)
    min_pixels = max(1, int(round(MIN_SEGMENT_AREA_NM2 / pixelsize_nm**2)))
    corners = particle_measures.count_ring_corners(
        points, MEAN_RADIUS_PCNA_NM, R_MIN_NM, R_MAX_NM, N_SEGMENTS,
        n_steps=360, min_per_section=min_pixels,
    )
    ring = points[corners["kept"]]
    max_gap = 360.0
    if len(ring) > 1:
        max_gap = particle_measures.angular_gaps(ring, corners["centre"]).max()
    n_peaks = particle_measures.count_resolved_sites(image, pixelsize_nm, MIN_SEPARATION_NM)
    return dict(n_corners=corners["n_corners"], n_peaks=min(n_peaks, N_SEGMENTS),
                max_gap_deg=max_gap, radius_nm=radius)


def generic_columns():
    return [(name, feature_names(name)) for name in GENERIC]


def measure_image(image, pixelsize_nm):
    """Generic and PCNA descriptors of one image, as {column: value}."""
    values = {}
    prepared = same_field_of_view(image, pixelsize_nm, FOV_NM)
    for name, columns in generic_columns():
        values.update(zip(columns, np.atleast_1d(GENERIC_FUNCTIONS[name](prepared))))
    values.update(measure_pcna(image, pixelsize_nm))
    return values


def pcna_metrics():
    """measure_image as custom metrics of a sweep, one per column; each image is measured once."""
    last = {"image": None, "values": None}

    def measure(image, pixelsize_nm):
        # the sweep calls every metric in turn with the same image object
        if last["image"] is not image:
            last["image"] = image
            last["values"] = measure_image(image, float(pixelsize_nm))
        return last["values"]

    columns = [c for _, names in generic_columns() for c in names] + SPECIFIC
    metrics = []
    for name in columns:
        def metric(simulated_image=None, simulated_image_pixelsize_nm=None, _name=name, **unused):
            return float(measure(simulated_image, simulated_image_pixelsize_nm)[_name])
        metric.__name__ = name
        metrics.append(metric)
    return metrics


def simulate(job):
    """One model, localisation precision and labelling efficiency: its feature table."""
    model, precision, le, seed, images_dir = job
    sweep = run_parameter_sweep(labelling_efficiency=[le], custom_metrics=pcna_metrics(),
                                random_seed=seed, **COMMON)
    # the sweep adds the modality from these parameters before imaging; the
    # SMLM image renders the localisations, so its resolution is the
    # localisation precision (not the PSF width of lateral_resolution_nm)
    emitters = sweep.experiment.local_modalities_parameters[MODALITY]["emitters"]
    emitters.update(lateral_precision=precision, axial_precision=precision,
                    nlocalisations=N_LOCALISATIONS, rendering_kernel_nm=RENDERING_KERNEL_NM)
    sweep.run_analysis(save=False, plots=False)
    if images_dir is not None:
        output_directory = os.path.join(
            images_dir, f"{model}_precision_{precision}_labelling_efficiency_{le}")
        os.makedirs(output_directory, exist_ok=True)
        sweep.save_images(output_directory=output_directory)
    table = sweep.analysis["dataframes"].copy()
    table["model"] = model
    table["localisation_precision_nm"] = precision
    table["labelling_efficiency"] = le
    return table


def feature_table(tables):
    """Fully labelled and reduced rows for every labelling efficiency."""
    features = pd.concat(tables, ignore_index=True)
    full = features[features.model == "full"]
    reduced = features[features.model == "reduced"]
    # the fully labelled trimers are shared by every labelling efficiency
    return pd.concat([reduced] + [full.assign(labelling_efficiency=le)
                                  for le in LABELLING_EFFICIENCIES], ignore_index=True)


def plot_feature_distributions(features, prefix):
    """Mean +- s.d. of each PCNA descriptor and brightness against labelling efficiency."""
    reduced = features[features.model == "reduced"]
    precisions = sorted(reduced.localisation_precision_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(precisions)))
    columns = SPECIFIC + ["brightness"]
    fig, axes = plt.subplots(1, len(columns), figsize=[3.5 * len(columns), 3])
    for ax, feature in zip(axes, columns):
        for precision, colour in zip(precisions, colours):
            summary = (reduced[np.isclose(reduced.localisation_precision_nm, precision)]
                       .groupby("labelling_efficiency")[feature].agg(["mean", "std"]))
            ax.errorbar(summary.index, summary["mean"], yerr=summary["std"], fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{precision:g}")
        ax.set_xlabel("Labelling efficiency")
        ax.set_ylabel(feature.replace("_", " "))
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].legend(title="Precision (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "feature_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_auc_vs_labelling(scores, prefix, feature_sets=("shape", "brightness", "specific", "all")):
    """AUC (95% interval) against labelling efficiency, one line per precision."""
    precisions = sorted(scores.localisation_precision_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(precisions)))
    fig, axes = plt.subplots(1, len(feature_sets), figsize=[3.5 * len(feature_sets), 3],
                             sharey=True)
    for ax, name in zip(axes, feature_sets):
        data = scores[scores.feature_set == name]
        for precision, colour in zip(precisions, colours):
            sub = data[np.isclose(data.localisation_precision_nm, precision)]
            sub = sub.sort_values("labelling_efficiency")
            ax.errorbar(sub.labelling_efficiency, sub.auc,
                        yerr=[sub.auc - sub.auc_low, sub.auc_high - sub.auc], fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{precision:g}")
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


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "fig2_pcna_labelling_")
        images_dir = prefix + "images" if SAVE_IMAGES else None
        blocks = [(model, precision, le) for model in MODELS for precision in PRECISIONS_NM
                  for le in MODELS[model]]
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [block + (int(seed), images_dir) for block, seed in zip(blocks, seeds)]
        print(f"{len(jobs)} sweeps, {N_REALISATIONS} realisations per model and condition, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                tables = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="sweeps"))
        else:
            tables = [simulate(job) for job in tqdm(jobs, desc="sweeps")]
        features = feature_table(tables)
        features.to_csv(prefix + "features.csv", index=False)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    plot_feature_distributions(features, prefix)
    features = features.fillna(FILL_VALUES)
    generic = [c for _, names in generic_columns() for c in names]
    feature_sets = dict(FEATURE_SETS, all=generic + SPECIFIC)
    scores = []
    for name, columns in feature_sets.items():
        print(f"features: {name}")
        table = distinguishability_from_features(
            features, "model", ("full", "reduced"), columns, by=CONDITIONS,
            controls=CONTROLS if name == "all" else None,
            replicate_column="Replica", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    scores = pd.concat(scores, ignore_index=True)
    # AUC of a perfect count of labelled sites: the ceiling of any feature
    scores["auc_ideal"] = 1 - scores.labelling_efficiency**N_SEGMENTS / 2
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    fig = plot_distinguishability_heatmaps(scores, x="labelling_efficiency",
                                           y="localisation_precision_nm", facet="feature_set",
                                           facet_order=list(feature_sets))
    n_le = scores.labelling_efficiency.nunique()
    n_precision = scores.localisation_precision_nm.nunique()
    panels = fig.axes[:len(feature_sets)]
    for ax in panels:
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Localisation precision (nm)")
    fig.suptitle("PCNA: fully labelled vs reduced labelling (ROC AUC)", x=0.02, ha="left")
    fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * n_precision)
    fig.savefig(prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    plot_auc_vs_labelling(scores, prefix)
    print(scores.pivot_table(index=["feature_set", "localisation_precision_nm"],
                             columns="labelling_efficiency", values="auc",
                             sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
