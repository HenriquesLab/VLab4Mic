"""Fig. 2 and Supp. Fig. S1: can domed and flat clathrin lattices be told apart by their aspect ratio and brightness?

Each lattice is simulated with run_parameter_sweep. Every repetition is an
independent realisation: one particle with its own random labelling,
structural integrity, orientation, position and noise. The features of
every image are computed by the sweep (feature_metrics) and collected, for
both lattices, in one feature table with one row per image. The table is
scored per condition with distinguishability_from_features: out-of-fold
ROC AUC of a classifier with a 95% bootstrap interval, a null control (two
random halves of the domes, expected 0.5) and a permutation p-value.

Features: aspect ratio and brightness. The aspect ratio is the ratio of
the principal axes of the particle image (smoothed, background subtracted,
thresholded): it measures elongation and is blind to size, orientation and
brightness. Brightness is the total intensity above background: it follows
the number of labelled sites. Each image is first brought to a common field
of view. Each feature is also scored alone (aspect_ratio_auc,
brightness_auc), to show what each contributes.

Geometry: clathrin pits and plaques lie flat on the coverslip at the
ventral surface of cultured cells, so the lattices are placed flat with a
random tilt up to MAX_TILT_DEG and a random in-plane rotation; probe
normals use the "local_plane" mode, correct on flat surfaces. Every
particle is also placed at a random position (POSITION_JITTER_NM), so that
particles do not all sit at the same position relative to the pixels.
Uniform random orientations are a supplementary control.

Conditions: a labelling efficiency x structural integrity grid (its point at
1, 1 is the reference), a linkage error sweep, the lateral and axial PSF
width of BLUR_MODALITY, and uniform orientations.

Outputs, in the VLab4Mic output directory (<date>fig2_clathrin_...):
features.csv (one row per image), distinguishability.csv (one row per
condition and imaging method), settings.json, and the figures.

Run:   python fig2_clathrin_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_clathrin_distinguishability.py --features <...features.csv>

Environment variables: VLAB4MIC_N_REALISATIONS (100 per lattice and
condition; e.g. 5 for a quick run), VLAB4MIC_N_PERMUTATIONS (200),
VLAB4MIC_N_WORKERS (1; sweeps and scoring in parallel processes, about
4 GB of memory each), VLAB4MIC_SAVE_IMAGES (0; 1 saves the simulated images
as TIFF, one folder per sweep).
"""

import json
import os
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tqdm.auto import tqdm

from vlab4mic.analysis.distinguishability import (
    distinguishability_from_features,
    feature_metrics,
    plot_distinguishability_heatmaps,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path
from vlab4mic.sweep_generator import run_parameter_sweep

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24

MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]
FEATURES = ["aspect_ratio", "brightness"]
CONTROLS = {"aspect_ratio": "aspect_ratio", "brightness": "brightness"}  # each alone
FOV_NM = 600  # lateral size of the virtual sample and common field of view
MAX_TILT_DEG = 10  # random tilt of the lattices on the coverslip
POSITION_JITTER_NM = 100  # random particle position, +-100 nm in x and y

# reference condition, and the values swept around it
LABELLING_EFFICIENCY = 1.0
LABELLING_EFFICIENCIES = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0]
STRUCTURAL_INTEGRITY = 1.0
STRUCTURAL_INTEGRITIES = [0.5, 0.7, 1.0]
# subunits: epitopes within 100 A; remaining lattice connected at 300 A
SI_CLUSTERS = dict(structural_integrity_small_cluster=[100],
                   structural_integrity_large_cluster=[300])
LINKAGE_ERROR_NM = 10  # distance to epitope of the Antibody probe template
LINKAGE_ERRORS_NM = [2, 5, 20]
BLUR_MODALITY = "STED"
LATERAL_BLUR_NM = [10, 20, 40, 80, 130]  # PSF s.d., axial at its default
AXIAL_BLUR_NM = [50, 100, 200, 300, 600]  # PSF s.d., lateral at its default

MODELS = {
    "dome": "https://zenodo.org/records/20719843/files/Dome_model_v2.0.cif",
    "flat": "https://zenodo.org/records/20719843/files/Flat_lattice_model_v2.0.cif",
}
COMMON = dict(  # arguments of every sweep
    probe_templates=["Antibody"], probe_target_type=["Sequence"],
    probe_target_value=["ATETQ"], probe_DoL=[4],
    sample_dimensions=[[FOV_NM, FOV_NM, 100]],
    structure_normal_orientation="local_plane",
    position_jitter_nm=POSITION_JITTER_NM,
    sweep_repetitions=N_REALISATIONS, use_reference=False, default_metrics=[],
    save_analysis_results=False, analysis_plots=False, return_generator=True,
)
VENTRAL = dict(random_orientations=False, random_rotations=True, orientation_tilt_max=MAX_TILT_DEG)
UNIFORM = dict(random_orientations=True, random_rotations=True)

# sweep blocks: name, geometry, swept parameters, imaging methods. The grid
# has one block per labelling efficiency, so blocks can run in parallel.
SWEEPS = [("grid", VENTRAL, dict(labelling_efficiency=[le], structural_integrity=STRUCTURAL_INTEGRITIES,
                                  **SI_CLUSTERS), MODALITIES)
          for le in LABELLING_EFFICIENCIES]
SWEEPS += [
    ("linkage_error", VENTRAL, dict(probe_distance_to_epitope=LINKAGE_ERRORS_NM), MODALITIES),
    ("uniform_orientation", UNIFORM, dict(labelling_efficiency=[LABELLING_EFFICIENCY]), MODALITIES),
    ("lateral_blur", VENTRAL, dict(lateral_resolution_nm=LATERAL_BLUR_NM), [BLUR_MODALITY]),
    ("axial_blur", VENTRAL, dict(axial_resolution_nm=AXIAL_BLUR_NM), [BLUR_MODALITY]),
]
CONDITIONS = ["sweep", "modality_name", "labelling_efficiency", "structural_integrity",
              "linkage_error_nm", "lateral_resolution_nm", "axial_resolution_nm"]


def download(url):
    path = Path(tempfile.gettempdir()) / Path(url).name
    if not path.exists():
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path.write_bytes(response.content)
    return str(path)


def resolution(experiment, modality):
    # effective lateral and axial resolution (nm): localisation precision
    # for SMLM, PSF standard deviation otherwise
    emitters = experiment.imaging_modalities[modality].get("emitters") or {}
    if emitters.get("lateral_precision") is not None:
        return emitters["lateral_precision"], emitters["axial_precision"]
    psf = experiment.imaging_modalities[modality]["psf_params"]
    lateral, _, axial = [sd * v for sd, v in zip(psf["std_devs"], psf["voxelsize"])]
    return lateral, axial


def simulate(job):
    """One lattice, one sweep block: its feature table (one row per image)."""
    model, path, (name, geometry, swept, modalities), seed, images_dir = job
    output_directory = None
    if images_dir is not None:
        label = "_".join([name, model] + [f"{k}_{v[0]}" for k, v in swept.items()
                                          if k == "labelling_efficiency"])
        output_directory = os.path.join(images_dir, label)
        os.makedirs(output_directory, exist_ok=True)
    sweep = run_parameter_sweep(
        structures=[path], structure_is_path=True, modalities=modalities,
        custom_metrics=feature_metrics(FEATURES, fov_nm=FOV_NM),
        random_seed=seed, save_sweep_images=images_dir is not None,
        output_directory=output_directory, **COMMON, **geometry, **swept,
    )
    table = sweep.analysis["dataframes"].copy()
    table["model"] = model
    table["sweep"] = name
    table = table.rename(columns={"probe_distance_to_epitope": "linkage_error_nm"})
    # parameters not swept in this block are at their reference value
    for column, value in (("labelling_efficiency", LABELLING_EFFICIENCY),
                          ("structural_integrity", STRUCTURAL_INTEGRITY),
                          ("linkage_error_nm", LINKAGE_ERROR_NM)):
        if column not in table:
            table[column] = value
    for k, column in enumerate(("lateral_resolution_nm", "axial_resolution_nm")):
        default = table["modality_name"].map(
            {m: resolution(sweep.experiment, m)[k] for m in modalities})
        table[column] = table[column].fillna(default) if column in table else default
    return table


def make_figures(scores, prefix):
    reference = scores[scores.sweep == "reference"].set_index("modality_name").loc[MODALITIES]
    uniform = scores[scores.sweep == "uniform_orientation"].set_index("modality_name").loc[MODALITIES]
    grid = scores[scores.sweep.isin(["grid", "reference"])]
    colours = {m: f"C{k}" for k, m in enumerate(MODALITIES)}

    def errorbars(ax, x, data, **kwargs):
        ax.errorbar(x, data.auc, yerr=[data.auc - data.auc_low, data.auc_high - data.auc],
                    capsize=2, **kwargs)

    def style(ax, xlabel, title):
        ax.axhline(0.5, color="0.5", linestyle="--", linewidth=1)
        ax.set_ylim(0.3, 1.05)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("ROC AUC")
        ax.set_title(title, loc="left")
        ax.spines[["top", "right"]].set_visible(False)

    fig, axes = plt.subplots(2, 3, figsize=[16, 9])
    axes = axes.ravel()
    x = np.arange(len(MODALITIES))
    errorbars(axes[0], x - 0.2, reference, fmt="o", label="aspect ratio + brightness")
    for offset, name, marker, colour in [(0.0, "aspect_ratio", "^", "C2"), (0.2, "brightness", "D", "C1")]:
        alone = reference[[f"{name}_auc", f"{name}_auc_low", f"{name}_auc_high"]].set_axis(
            ["auc", "auc_low", "auc_high"], axis=1)
        errorbars(axes[0], x + offset, alone, fmt=marker, color=colour, label=f"{name.replace('_', ' ')} alone")
    axes[0].scatter(x + 0.35, reference.null_auc, marker="x", color="k", label="null (dome vs dome)")
    axes[0].set_xticks(x, MODALITIES)
    axes[0].legend(frameon=False)
    style(axes[0], "", "a  imaging method (reference)")
    errorbars(axes[1], x - 0.1, reference, fmt="o", label="ventral")
    errorbars(axes[1], x + 0.1, uniform, fmt="s", label="uniform orientation")
    axes[1].set_xticks(x, MODALITIES)
    axes[1].legend(frameon=False)
    style(axes[1], "", "b  ventral vs uniform orientations")
    for ax, column, fixed, value, title in [
        (axes[2], "labelling_efficiency", "structural_integrity", STRUCTURAL_INTEGRITY,
         "c  labelling efficiency (integrity 1)"),
        (axes[3], "structural_integrity", "labelling_efficiency", LABELLING_EFFICIENCY,
         "d  structural integrity (labelling 1)"),
    ]:
        data = grid[np.isclose(grid[fixed], value)]
        for m in MODALITIES:
            sub = data[data.modality_name == m].sort_values(column)
            errorbars(ax, sub[column], sub, fmt="o-", color=colours[m], label=m)
        style(ax, column.replace("_", " "), title)
    axes[2].set_xscale("log")
    axes[2].legend(frameon=False, fontsize=8)
    data = scores[scores.sweep.isin(["linkage_error", "reference"])]
    for m in MODALITIES:
        sub = data[data.modality_name == m].sort_values("linkage_error_nm")
        errorbars(axes[4], sub.linkage_error_nm, sub, fmt="o-", color=colours[m], label=m)
    style(axes[4], "linkage error (nm)", "e  linkage error")
    for sweep, column, label in [("lateral_blur", "lateral_resolution_nm", "lateral"),
                                 ("axial_blur", "axial_resolution_nm", "axial")]:
        sub = scores[scores.sweep == sweep].sort_values(column)
        errorbars(axes[5], sub[column], sub, fmt="o-", label=f"{label} PSF width")
    axes[5].set_xscale("log")
    axes[5].legend(frameon=False)
    style(axes[5], f"{BLUR_MODALITY} PSF width, s.d. (nm)", "f  width of the blur")
    fig.tight_layout()
    fig.savefig(prefix + "distinguishability.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig = plot_distinguishability_heatmaps(grid, x="labelling_efficiency",
                                           y="structural_integrity", facet="modality_name",
                                           facet_order=MODALITIES)
    fig.savefig(prefix + "distinguishability_heatmaps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "fig2_clathrin_")
        images_dir = prefix + "images" if SAVE_IMAGES else None
        paths = {name: download(url) for name, url in MODELS.items()}
        blocks = [(model, path, block) for block in SWEEPS for model, path in paths.items()]
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [block + (int(seed), images_dir) for block, seed in zip(blocks, seeds)]
        print(f"{len(jobs)} sweeps, {N_REALISATIONS} realisations per lattice and condition, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                tables = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="sweeps"))
        else:
            tables = [simulate(job) for job in tqdm(jobs, desc="sweeps")]
        features = pd.concat(tables, ignore_index=True)
        reference = ((features.sweep == "grid")
                     & np.isclose(features.labelling_efficiency, LABELLING_EFFICIENCY)
                     & np.isclose(features.structural_integrity, STRUCTURAL_INTEGRITY))
        features.loc[reference, "sweep"] = "reference"
        features.to_csv(prefix + "features.csv", index=False)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    scores = distinguishability_from_features(
        features, "model", ("dome", "flat"), FEATURES, by=CONDITIONS, controls=CONTROLS,
        replicate_column="Replica", n_permutations=N_PERMUTATIONS,
        random_state=RANDOM_SEED, n_workers=N_WORKERS,
    )
    scores.to_csv(prefix + "distinguishability.csv", index=False)
    make_figures(scores, prefix)
    print(scores[["sweep", "modality_name", "labelling_efficiency", "structural_integrity",
                  "auc", "auc_low", "auc_high", "null_auc", "p_value"]].to_string())


if __name__ == "__main__":
    main()
