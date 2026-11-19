"""Fig. 2: can an incomplete nuclear pore be told apart from an intact one, and at which labelling efficiency?

The nuclear pore complex (7R5K) is labelled at the C-terminus of Nup96 with
a single direct probe (anti-GFP nanobody, no secondary) and imaged by SMLM
and STED, seen from above. Incomplete pores are not made by hand: they come
from the structural integrity model of VLab4Mic, which removes a connected
block of subunits around a random fracture point, leaving the fraction
structural_integrity of the pore. The pores are simulated with
run_parameter_sweep over labelling efficiency x structural integrity:
every repetition is one pore with its own labelling, defects, in-plane
rotation and noise.

Every image is measured from its pixels (the sweep gives images only): the
pixels within REGION_NM of the image centre above an Otsu threshold are
taken as points of the ring, and
- radius_nm: radius of a geometric circle fit to these points;
- n_segments: occupied segments of the 8-fold ring (count_ring_corners:
  points in a radial window around the ring, a segment needs at least
  MIN_SEGMENT_AREA_NM2 of bright pixels);
- max_gap_deg: the largest gap in azimuth between ring points
  (angular_gaps). Missing subunits of an incomplete pore are contiguous,
  so they leave one wide gap, while a low labelling efficiency leaves
  several small, scattered gaps.

The feature table (one row per image) is scored per imaging method,
labelling efficiency and structural integrity with
distinguishability_from_features: out-of-fold ROC AUC of a classifier
(intact vs incomplete) with a 95% bootstrap interval, a null control (two
random halves of the intact pores, expected 0.5) and a permutation p-value,
for each feature set (each metric alone, and all together). The intact
pores are simulated separately and shared by every structural integrity: a
structural integrity of 1 compares two independent sets of intact pores, a
negative control (expected 0.5).

Outputs, in the VLab4Mic output directory (<date>fig2_npc_integrity_...):
features.csv (one row per image), distinguishability.csv (one row per
imaging method, labelling efficiency, structural integrity and feature
set), settings.json, feature_distributions.png and
<modality>_distinguishability_heatmap.png.

Run:   python fig2_npc_integrity_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_npc_integrity_distinguishability.py --features <...features.csv>

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
    distinguishability_from_features,
    plot_distinguishability_heatmaps,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path
from vlab4mic.sweep_generator import run_parameter_sweep

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24

STRUCTURE = "7R5K"
FOV_NM = 600  # lateral size of the virtual sample
MODALITIES = ["SMLM", "STED"]

LABELLING_EFFICIENCIES = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0]
STRUCTURAL_INTEGRITIES = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # 1: intact (negative control)
# subunits: epitopes within 300 A; remaining pore connected at 600 A (7R5K.yaml)
SI_CLUSTERS = dict(structural_integrity_small_cluster=[300],
                   structural_integrity_large_cluster=[600])

# STED as in figS2_npc_three_labels.py, at a fixed exposure for every
# condition: fewer labels give a dimmer pore
STED_PARAMS = dict(lateral_resolution_nm=[15], axial_resolution_nm=[20], psf_voxel_nm=[5],
                   depth_of_field_nm=[100], pixelsize_nm=[15], exp_time=[0.0002])
STED_BINOMIAL_P = 0.7
STED_READ_NOISE = 0.2  # standard deviation of the Gaussian (read) noise
MODALITY_PARAMS = {"SMLM": {}, "STED": STED_PARAMS}

# measurement of the pore
REGION_NM = 100  # the pore is at the image centre: ignore pixels further out
N_SEGMENTS = 8
RADIAL_WINDOW_NM = 25  # ring points: within +- this distance of the fitted radius
MIN_SEGMENT_AREA_NM2 = 50  # bright pixels a segment needs to count as detected

FEATURES = ["n_segments", "radius_nm", "max_gap_deg"]
FEATURE_SETS = {"segments": ["n_segments"], "radius": ["radius_nm"],
                "largest gap": ["max_gap_deg"], "all": FEATURES}
CONTROLS = {feature: feature for feature in FEATURES}  # each alone
# images without a measurable ring (e.g. too few labels): nothing detected
FILL_VALUES = dict(n_segments=0, radius_nm=0, max_gap_deg=360)
CONDITIONS = ["modality_name", "labelling_efficiency", "structural_integrity"]

COMMON = dict(  # arguments of every sweep
    structures=[STRUCTURE],
    probe_templates=["GFP_w_nanobody"], probe_target_type=["Sequence"],
    probe_target_value=["ELAVGSL"], probe_DoL=[2],  # Nup96 C-terminus
    sample_dimensions=[[FOV_NM, FOV_NM, 100]],
    random_orientations=False, random_rotations=True,  # seen from above, rotated in plane
    sweep_repetitions=N_REALISATIONS, use_reference=False, default_metrics=[],
    save_analysis_results=False, analysis_plots=False, save_sweep_images=False,
    run_analysis=False, return_generator=True,
)
MODELS = {"intact": [1.0], "incomplete": STRUCTURAL_INTEGRITIES}


def measure_npc(image, pixelsize_nm):
    """Ring radius, detected segments and largest angular gap of the pore at the image centre."""
    missing = dict(n_segments=np.nan, radius_nm=np.nan, max_gap_deg=np.nan)
    image = np.asarray(image, dtype=float)
    rows_px, cols_px = np.indices(image.shape)
    xy = np.c_[cols_px.ravel() + 0.5, rows_px.ravel() + 0.5] * pixelsize_nm  # pixel centres, nm
    values = image.ravel()
    image_centre = np.array(image.shape[::-1]) * pixelsize_nm / 2
    region = np.linalg.norm(xy - image_centre, axis=1) <= REGION_NM
    xy, values = xy[region], values[region]
    if np.ptp(values) == 0:
        return missing
    points = xy[values > particle_measures.compute_threshold(values, method="otsu")]
    if len(points) < 3:
        return missing
    _, radius = particle_measures.fit_circle_geometric(points)
    min_pixels = max(1, int(round(MIN_SEGMENT_AREA_NM2 / pixelsize_nm**2)))
    corners = particle_measures.count_ring_corners(
        points, radius, radius - RADIAL_WINDOW_NM, radius + RADIAL_WINDOW_NM, N_SEGMENTS,
        n_steps=360, min_per_section=min_pixels,
    )
    ring = points[corners["kept"]]
    max_gap = 360.0
    if len(ring) > 1:
        max_gap = particle_measures.angular_gaps(ring, corners["centre"]).max()
    return dict(n_segments=corners["n_corners"], radius_nm=radius, max_gap_deg=max_gap)


def npc_metrics():
    """measure_npc as custom metrics of a sweep, one per feature; each image is measured once."""
    last = {"image": None, "values": None}

    def measure(image, pixelsize_nm):
        # the sweep calls every metric in turn with the same image object
        if last["image"] is not image:
            last["image"] = image
            last["values"] = measure_npc(image, float(pixelsize_nm))
        return last["values"]

    metrics = []
    for name in FEATURES:
        def metric(simulated_image=None, simulated_image_pixelsize_nm=None, _name=name, **unused):
            return float(measure(simulated_image, simulated_image_pixelsize_nm)[_name])
        metric.__name__ = name
        metrics.append(metric)
    return metrics


def simulate(job):
    """One model, imaging method and labelling efficiency: its feature table (one row per image)."""
    model, modality, le, seed, images_dir = job
    sweep = run_parameter_sweep(
        modalities=[modality], labelling_efficiency=[le],
        structural_integrity=MODELS[model], **SI_CLUSTERS,
        custom_metrics=npc_metrics(), random_seed=seed,
        **COMMON, **MODALITY_PARAMS[modality],
    )
    if modality == "STED":
        # the sweep adds the modality from these parameters before imaging
        noise_model = sweep.experiment.local_modalities_parameters["STED"]["detector"]["noise_model"]
        noise_model["binomial"]["p"] = STED_BINOMIAL_P
        noise_model["gaussian"]["sigma"] = STED_READ_NOISE
    sweep.run_analysis(save=False, plots=False)
    if images_dir is not None:
        output_directory = os.path.join(images_dir, f"{model}_{modality}_labelling_efficiency_{le}")
        os.makedirs(output_directory, exist_ok=True)
        sweep.save_images(output_directory=output_directory)
    table = sweep.analysis["dataframes"].copy()
    table["model"] = model
    table["modality_name"] = modality
    if "labelling_efficiency" not in table:
        table["labelling_efficiency"] = le
    if "structural_integrity" not in table:
        table["structural_integrity"] = MODELS[model][0]
    return table


def feature_table(tables):
    """Intact and incomplete rows for every structural integrity."""
    features = pd.concat(tables, ignore_index=True)
    intact = features[features.model == "intact"]
    incomplete = features[features.model == "incomplete"]
    # the intact pores are shared by every structural integrity
    return pd.concat([incomplete] + [intact.assign(structural_integrity=si)
                                     for si in STRUCTURAL_INTEGRITIES], ignore_index=True)


def plot_feature_distributions(features, prefix):
    """Mean +- s.d. of each feature against structural integrity, per labelling efficiency."""
    incomplete = features[features.model == "incomplete"]
    les = sorted(incomplete.labelling_efficiency.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(les)))
    fig, axes = plt.subplots(len(MODALITIES), len(FEATURES),
                             figsize=[4 * len(FEATURES), 3 * len(MODALITIES)], squeeze=False)
    for row, modality in zip(axes, MODALITIES):
        data = incomplete[incomplete.modality_name == modality]
        for ax, feature in zip(row, FEATURES):
            for le, colour in zip(les, colours):
                summary = (data[np.isclose(data.labelling_efficiency, le)]
                           .groupby("structural_integrity")[feature].agg(["mean", "std"]))
                ax.errorbar(summary.index, summary["mean"], yerr=summary["std"], fmt="o-",
                            color=colour, capsize=2, markersize=3, label=f"{le:g}")
            ax.set_xlabel("Structural integrity")
            ax.set_ylabel(feature.replace("_", " "))
            ax.set_title(modality, loc="left")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, -1].legend(title="Labelling efficiency", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "feature_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "fig2_npc_integrity_")
        images_dir = prefix + "images" if SAVE_IMAGES else None
        blocks = [(model, modality, le) for model in MODELS for modality in MODALITIES
                  for le in LABELLING_EFFICIENCIES]
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
    scores = []
    for name, columns in FEATURE_SETS.items():
        print(f"features: {name}")
        table = distinguishability_from_features(
            features, "model", ("intact", "incomplete"), columns, by=CONDITIONS,
            controls=CONTROLS if name == "all" else None,
            replicate_column="Replica", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    scores = pd.concat(scores, ignore_index=True)
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    for modality in MODALITIES:
        data = scores[scores.modality_name == modality]
        if data.empty:
            continue
        fig = plot_distinguishability_heatmaps(data, x="labelling_efficiency",
                                               y="structural_integrity", facet="feature_set",
                                               facet_order=list(FEATURE_SETS))
        n_le = data.labelling_efficiency.nunique()
        n_si = data.structural_integrity.nunique()
        panels = fig.axes[:len(FEATURE_SETS)]
        for ax in panels:
            ax.set_xlabel("Labelling efficiency")
        panels[0].set_ylabel("Structural integrity")
        fig.suptitle(f"{modality}: intact vs incomplete nuclear pore (ROC AUC)", x=0.02, ha="left")
        fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * n_si)
        fig.savefig(prefix + f"{modality}_distinguishability_heatmap.png", dpi=300,
                    bbox_inches="tight")
        plt.close(fig)
    print(scores.pivot_table(index=["feature_set", "modality_name", "structural_integrity"],
                             columns="labelling_efficiency", values="auc",
                             sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
