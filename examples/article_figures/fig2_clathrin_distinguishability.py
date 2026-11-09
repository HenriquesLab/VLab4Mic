"""Fig. 2 and Supp. Fig. S1: distinguishability of domed and flat clathrin lattices.

Structural distinguishability is the accuracy with which a classifier tells
two candidate structures apart under stated simulation parameters. Each
lattice is simulated with vlab4mic.sweep_generator.run_parameter_sweep:
every repetition is an independent realisation (one particle, with its own
random labelling, structural integrity, orientation and noise). Features
of every simulated image are computed as custom sweep metrics and written,
for both lattices, to one feature table (<date>fig2_clathrin_features.csv,
one row per image). The table is then scored with
vlab4mic.analysis.distinguishability: out-of-fold accuracy and ROC AUC of a
classifier, with 95% bootstrap intervals.

The feature table can be scored again without simulating, e.g. to try
other features or classifiers::

    python fig2_clathrin_distinguishability.py --features <date>fig2_clathrin_features.csv

Features (all computed on each image padded or cropped to the same field of
view, FOV_NM, so they do not depend on the image size of a condition):

- shape_00 ... shape_19: vlab4mic.analysis.distinguishability.shape_features
  (radial profile, rotation-averaged power spectrum, principal axes) of the
  background-subtracted image scaled to unit total intensity, so the
  lattices are not told apart by brightness;
- aspect_ratio: size- and brightness-invariant aspect ratio (control);
- brightness: total intensity above background (control).

Classifier: shape features, standardisation, PCA and logistic regression,
5-fold stratified cross-validation. Every score comes with controls: the
AUC from aspect_ratio alone and from brightness alone, null_auc (two random
halves of the domes, expected 0.5) and p_value (permutation p-value of the
shape AUC).

Geometry (Leterrier, co-author): clathrin pits and plaques are imaged at the
ventral surface of cultured cells, lying flat on the coverslip, so the two
forms differ mainly in depth. The lattices lie on the coverslip with a
random tilt up to MAX_TILT_DEG and a random in-plane rotation; probe
normals use the "local_plane" mode, correct on flat surfaces. Uniform
random orientations are a supplementary control.

Conditions (ventral geometry unless stated):

- labelling efficiency x structural integrity grid, every imaging method:
  AUC against each at 1 of the other (Fig. 2c,d) and one heatmap per
  method; its (1, 1) point is the reference (Fig. 2a,b,g);
- linkage error sweep (distance from fluorophore to epitope), every
  imaging method (Fig. 2e);
- width of the microscope's blur: lateral and axial PSF width of
  BLUR_MODALITY swept separately (Fig. 2f);
- uniform random orientations, every imaging method (Supp. Fig. S1).

Environment variables (for a quick run, e.g. VLAB4MIC_N_REALISATIONS=5
VLAB4MIC_N_PERMUTATIONS=0):

- VLAB4MIC_N_REALISATIONS: realisations per lattice and condition (100);
- VLAB4MIC_N_PERMUTATIONS: permutations for the p-value (200);
- VLAB4MIC_MAX_TILT_DEG: maximum tilt on the coverslip in degrees (10);
- VLAB4MIC_N_WORKERS: sweeps (and then scoring) run in parallel, one
  process each (1). Each worker loads one lattice model (about 4 GB of
  memory). Every sweep has its own seed, so results do not depend on this
  number;
- VLAB4MIC_SAVE_IMAGES: 1 to also save the simulated images as TIFF, one
  folder per sweep in <date>fig2_clathrin_images/ (0).
"""

import contextlib
import json
import os
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tqdm.auto import tqdm

from vlab4mic.analysis.distinguishability import (
    aspect_ratio_features,
    brightness_features,
    distinguishability_from_dataframes,
    shape_features,
)
from vlab4mic.sweep_generator import run_parameter_sweep

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
# parallel worker processes, one sweep (or scoring group) each at a time
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24
# maximum random tilt of the lattices on the coverslip (degrees); no
# published reference for this value
MAX_TILT_DEG = float(os.environ.get("VLAB4MIC_MAX_TILT_DEG", 10))
MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]
# lateral size of the virtual sample, and field of view of every image
# before its features are computed (nm)
FOV_NM = 600

# reference condition, and the values swept around it
LABELLING_EFFICIENCY = 1.0
# down to a few labelled sites per particle: antibody crowding caps the
# number of labels at high efficiency, so the outline is lost only at the
# low end
LABELLING_EFFICIENCIES = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0]
# structural integrity: fraction of the lattice left after removing a
# contiguous patch of subunits. Subunits are epitopes within
# SI_SMALL_CLUSTER_A of each other, and the remaining lattice must stay
# connected at SI_LARGE_CLUSTER_A (both in angstrom). These distances keep
# the requested fraction on both lattices; much smaller or larger values
# fail on the clathrin lattices. Below 0.5 the removal often misses the
# requested fraction.
STRUCTURAL_INTEGRITY = 1.0
STRUCTURAL_INTEGRITIES = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
SI_SMALL_CLUSTER_A = 100
SI_LARGE_CLUSTER_A = 300
LINKAGE_ERROR_NM = 10  # distance to epitope of the Antibody probe template
LINKAGE_ERRORS_NM = [2, 5, 10, 20]  # direct label / nanobody to primary+secondary
BLUR_MODALITY = "STED"  # PSF-based method whose PSF width is swept
LATERAL_BLUR_NM = [10, 20, 40, 80, 130]  # PSF standard deviation, axial at default
AXIAL_BLUR_NM = [50, 100, 200, 300, 600]  # PSF standard deviation, lateral at default

MODELS = {
    "dome": "https://zenodo.org/records/20719843/files/Dome_model_v2.0.cif",
    "flat": "https://zenodo.org/records/20719843/files/Flat_lattice_model_v2.0.cif",
}
# probe as run_parameter_sweep arguments (one value each)
PROBE = dict(
    probe_templates=["Antibody"],
    probe_target_type=["Sequence"],
    probe_target_value=["ATETQ"],
    probe_DoL=[4],
)
GEOMETRIES = {
    "ventral": dict(random_orientations=False, random_rotations=True,
                    orientation_tilt_max=MAX_TILT_DEG),
    "uniform": dict(random_orientations=True, random_rotations=True),
}

# features: columns of the feature table
SHAPE_COLUMNS = [f"shape_{k:02d}" for k in range(20)]
FEATURE_COLUMNS = SHAPE_COLUMNS + ["aspect_ratio", "brightness"]
CONTROL_FEATURES = ["aspect_ratio", "brightness"]
# columns that identify a condition in the feature table
CONDITION_COLUMNS = ["sweep", "geometry", "modality", "labelling_efficiency",
                     "structural_integrity", "linkage_error_nm",
                     "lateral_resolution_nm", "axial_resolution_nm"]


def download(url):
    path = Path(tempfile.gettempdir()) / Path(url).name
    if not path.exists():
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path.write_bytes(response.content)
    return str(path)


@contextlib.contextmanager
def quiet():
    # hide the per-simulation messages, inner progress bars and warnings of
    # the library, so the progress bars stay readable; errors are still
    # raised (and reported by the main process)
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull), \
            contextlib.redirect_stderr(devnull), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# features as custom sweep metrics ----------------------------------------
def same_field_of_view(image, pixelsize_nm):
    # pad (with the median) or crop the image, about its centre, to FOV_NM
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        image = image.sum(axis=0)
    n = max(int(round(FOV_NM / pixelsize_nm)), 1)
    out = np.full((n, n), np.median(image))
    src, dst = [], []
    for size in image.shape:
        if size >= n:
            start = (size - n) // 2
            src.append(slice(start, start + n))
            dst.append(slice(0, n))
        else:
            start = (n - size) // 2
            src.append(slice(0, size))
            dst.append(slice(start, start + size))
    out[tuple(dst)] = image[tuple(src)]
    return out


_last_image = [None, None]  # the image object last measured and its features


def image_features(image, pixelsize_nm):
    # the sweep calls every metric with the same image object in turn:
    # compute all features once per image (the reference to the image is
    # kept, so its id cannot be reused by another image)
    if _last_image[0] is not image:
        prepared = same_field_of_view(image, pixelsize_nm)
        _last_image[0] = image
        _last_image[1] = np.concatenate([
            shape_features(prepared),
            aspect_ratio_features(prepared),
            brightness_features(prepared),
        ])
    return _last_image[1]


def make_feature_metric(index, name):
    def metric(simulated_image=None, simulated_image_pixelsize_nm=None, **kwargs):
        return float(image_features(simulated_image, simulated_image_pixelsize_nm)[index])
    metric.__name__ = name
    return metric


FEATURE_METRICS = [make_feature_metric(k, name) for k, name in enumerate(FEATURE_COLUMNS)]


# simulation: one run_parameter_sweep per lattice and sweep block ----------
def list_sweeps():
    # every sweep block: (label, geometry, run_parameter_sweep arguments,
    # modalities). The grid is split by labelling efficiency so that blocks
    # can run in parallel; the slow blocks come first.
    integrity = dict(structural_integrity=STRUCTURAL_INTEGRITIES,
                     structural_integrity_small_cluster=[SI_SMALL_CLUSTER_A],
                     structural_integrity_large_cluster=[SI_LARGE_CLUSTER_A])
    sweeps = [(f"grid, labelling efficiency {le}", "ventral",
               dict(labelling_efficiency=[le], **integrity), MODALITIES)
              for le in LABELLING_EFFICIENCIES]
    sweeps += [
        ("linkage error", "ventral",
         dict(labelling_efficiency=[LABELLING_EFFICIENCY],
              probe_distance_to_epitope=[d for d in LINKAGE_ERRORS_NM if d != LINKAGE_ERROR_NM]),
         MODALITIES),
        ("uniform orientation", "uniform",
         dict(labelling_efficiency=[LABELLING_EFFICIENCY]), MODALITIES),
        (f"{BLUR_MODALITY} lateral blur", "ventral",
         dict(labelling_efficiency=[LABELLING_EFFICIENCY], lateral_resolution_nm=LATERAL_BLUR_NM),
         [BLUR_MODALITY]),
        (f"{BLUR_MODALITY} axial blur", "ventral",
         dict(labelling_efficiency=[LABELLING_EFFICIENCY], axial_resolution_nm=AXIAL_BLUR_NM),
         [BLUR_MODALITY]),
    ]
    return sweeps


def resolution(experiment, modality):
    # effective lateral and axial resolution (nm): localisation precision
    # for localisation-based methods, PSF standard deviation otherwise
    emitters = experiment.imaging_modalities[modality].get("emitters") or {}
    if emitters.get("lateral_precision") is not None:
        return emitters["lateral_precision"], emitters["axial_precision"]
    psf = experiment.imaging_modalities[modality]["psf_params"]
    lateral, _, axial = [sd * v for sd, v in zip(psf["std_devs"], psf["voxelsize"])]
    return lateral, axial


def run_sweep(model, model_path, label, geometry, arguments, modalities, seed, images_dir):
    # one lattice, one sweep block: per-image feature table
    output_directory = None
    if images_dir is not None:
        output_directory = os.path.join(images_dir, f"{label.replace(', ', '_').replace(' ', '_')}_{model}")
        os.makedirs(output_directory, exist_ok=True)
    with quiet():
        generator = run_parameter_sweep(
            structures=[model_path],
            structure_is_path=True,
            use_reference=False,
            modalities=modalities,
            sweep_repetitions=N_REALISATIONS,
            sample_dimensions=[[FOV_NM, FOV_NM, 100]],
            structure_normal_orientation="local_plane",
            custom_metrics=FEATURE_METRICS,
            default_metrics=[],
            random_seed=seed,
            output_directory=output_directory,
            save_sweep_images=images_dir is not None,
            save_analysis_results=False,
            analysis_plots=False,
            run_analysis=True,
            return_generator=True,
            **PROBE,
            **GEOMETRIES[geometry],
            **arguments,
        )
        table = generator.analysis["dataframes"].copy()
        defaults = {m: resolution(generator.experiment, m) for m in modalities}
    table = table.rename(columns={"probe_distance_to_epitope": "linkage_error_nm"})
    table["modality"] = table["modality_name"]
    table["model"] = model
    table["geometry"] = geometry
    table["sweep"] = label.split(",")[0].replace(f"{BLUR_MODALITY} ", "").replace(" ", "_")
    # parameters not swept in this block are at their reference value
    for column, value in (("labelling_efficiency", LABELLING_EFFICIENCY),
                          ("structural_integrity", STRUCTURAL_INTEGRITY),
                          ("linkage_error_nm", LINKAGE_ERROR_NM)):
        if column not in table:
            table[column] = value
    for k, column in enumerate(("lateral_resolution_nm", "axial_resolution_nm")):
        default = table["modality"].map(lambda m: defaults[m][k])
        table[column] = table[column].fillna(default) if column in table else default
    return table


def simulate(paths, n_workers, images_dir):
    sweeps = list_sweeps()
    seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(sweeps) * len(paths))
    jobs = []
    for k, (label, geometry, arguments, modalities) in enumerate(sweeps):
        for j, (model, path) in enumerate(paths.items()):
            jobs.append((model, path, label, geometry, arguments, modalities,
                         int(seeds[k * len(paths) + j]), images_dir))
    tables = run_in_parallel(run_sweep, jobs, n_workers, "sweeps",
                             labels=[f"{job[2]} ({job[0]})" for job in jobs])
    features = pd.concat(tables, ignore_index=True)
    reference = ((features.sweep == "grid")
                 & np.isclose(features.labelling_efficiency, LABELLING_EFFICIENCY)
                 & np.isclose(features.structural_integrity, STRUCTURAL_INTEGRITY))
    features.loc[reference, "sweep"] = "reference"
    keep = ["model"] + CONDITION_COLUMNS + ["Replica"] + FEATURE_COLUMNS
    return features[keep]


def run_in_parallel(function, jobs, n_workers, description, labels):
    # run function(*job) for every job, in order of the results; one process
    # per job at a time when n_workers > 1
    results = [None] * len(jobs)
    progress = tqdm(total=len(jobs), desc=description, unit=description[:-1])
    if n_workers <= 1:
        for k, (job, label) in enumerate(zip(jobs, labels)):
            progress.set_postfix_str(label)
            results[k] = function(*job)
            progress.update()
    else:
        from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

        # one thread per worker for numpy / scipy, so workers do not compete
        for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.setdefault(variable, "1")
        with ProcessPoolExecutor(min(n_workers, len(jobs))) as pool:
            futures = {pool.submit(function, *job): (k, label)
                       for k, (job, label) in enumerate(zip(jobs, labels))}
            pending = set(futures)
            while pending:
                # wake up regularly so the elapsed time keeps moving
                done, pending = wait(pending, timeout=5, return_when=FIRST_COMPLETED)
                for future in done:
                    k, label = futures[future]
                    results[k] = future.result()
                    progress.set_postfix_str(f"done: {label}")
                    progress.update()
                progress.refresh()
    progress.close()
    return results


# scoring: one distinguishability score per condition ---------------------
def interval_columns(prefix, score):
    return {
        f"{prefix}accuracy": score["accuracy"],
        f"{prefix}accuracy_low": score["accuracy_interval"][0],
        f"{prefix}accuracy_high": score["accuracy_interval"][1],
        f"{prefix}auc": score["auc"],
        f"{prefix}auc_low": score["auc_interval"][0],
        f"{prefix}auc_high": score["auc_interval"][1],
    }


def score_condition(condition, dome, flat):
    # dome, flat: feature tables of one condition
    dome = dome.sort_values("Replica")
    flat = flat.sort_values("Replica")
    with quiet():
        score = distinguishability_from_dataframes(
            dome[SHAPE_COLUMNS].to_numpy(), flat[SHAPE_COLUMNS].to_numpy(),
            random_state=RANDOM_SEED, n_permutations=N_PERMUTATIONS,
        )
        row = dict(condition, n_per_class=min(len(dome), len(flat)), features="shape")
        row.update(interval_columns("", score))
        row.update(null_auc=score["null_auc"], p_value=score["p_value"])
        for control in CONTROL_FEATURES:
            control_score = distinguishability_from_dataframes(
                dome[[control]].to_numpy(), flat[[control]].to_numpy(),
                random_state=RANDOM_SEED, null_control=False,
            )
            row.update(interval_columns(f"{control}_", control_score))
    return row, dict(scores=score["scores"], labels=score["labels"])


def score_features(features, n_workers):
    # one row of scores per condition and modality, and the out-of-fold
    # scores of the reference condition per modality (Fig. 2a)
    jobs, labels = [], []
    for key, group in features.groupby(CONDITION_COLUMNS, sort=False):
        condition = dict(zip(CONDITION_COLUMNS, key))
        jobs.append((condition, group[group.model == "dome"], group[group.model == "flat"]))
        labels.append(f"{condition['sweep']} {condition['modality']}")
    outputs = run_in_parallel(score_condition, jobs, n_workers, "scores", labels)
    results = pd.DataFrame([row for row, _ in outputs])
    reference_scores = {row["modality"]: scores for row, scores in outputs
                        if row["sweep"] == "reference"}
    return results, reference_scores


# figures ----------------------------------------------------------------
def errorbars(ax, x, data, prefix="", **kwargs):
    auc = data[f"{prefix}auc"]
    yerr = [auc - data[f"{prefix}auc_low"], data[f"{prefix}auc_high"] - auc]
    ax.errorbar(x, auc, yerr=yerr, capsize=2, **kwargs)


def style(ax, xlabel):
    ax.axhline(0.5, color="0.5", linestyle="--", linewidth=1)
    ax.set_ylim(0.3, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("ROC AUC")
    ax.spines[["top", "right"]].set_visible(False)


def grid_results(results):
    # the labelling efficiency x structural integrity grid, reference included
    return results[results.sweep.isin(["grid", "reference"])]


def make_figures(results, reference_scores, default_lateral, default_axial, out_prefix):
    colours = {m: f"C{k}" for k, m in enumerate(MODALITIES)}
    grid = grid_results(results)
    fig, axes = plt.subplots(2, 4, figsize=[21, 9])
    axes = axes.ravel()

    # a: out-of-fold scores per imaging method (reference)
    for k, modality in enumerate(MODALITIES):
        s = reference_scores[modality]
        for label, offset, colour in [(0, -0.12, "C0"), (1, 0.12, "C1")]:
            values = s["scores"][s["labels"] == label]
            axes[0].scatter(np.full(values.size, k) + offset, values, s=5, color=colour,
                            alpha=0.6, label=["dome", "flat"][label] if k == 0 else None)
    axes[0].set_xticks(range(len(MODALITIES)), MODALITIES)
    axes[0].set_ylabel("P(flat), out-of-fold")
    axes[0].set_title("a  classifier scores", loc="left")
    axes[0].legend(frameon=False)

    # b: AUC against the axial resolution each method delivers (reference)
    reference = results[results.sweep == "reference"]
    for modality in MODALITIES:
        sub = reference[reference.modality == modality]
        errorbars(axes[1], sub.axial_resolution_nm, sub, fmt="o", color=colours[modality],
                  label=modality)
    axes[1].set_xscale("log")
    style(axes[1], "axial resolution (nm)")
    axes[1].set_title("b  imaging method", loc="left")
    axes[1].legend(frameon=False, fontsize=8)

    # c, d: labelling efficiency at integrity 1 and structural integrity at
    # labelling efficiency 1 (row and column of the grid)
    for ax, column, fixed, fixed_value, xlabel, title in [
        (axes[2], "labelling_efficiency", "structural_integrity", STRUCTURAL_INTEGRITY,
         "labelling efficiency", f"c  labelling efficiency (integrity {STRUCTURAL_INTEGRITY:g})"),
        (axes[3], "structural_integrity", "labelling_efficiency", LABELLING_EFFICIENCY,
         "structural integrity", f"d  structural integrity (labelling {LABELLING_EFFICIENCY:g})"),
    ]:
        data = grid[np.isclose(grid[fixed], fixed_value)]
        for modality in MODALITIES:
            sub = data[data.modality == modality].sort_values(column)
            errorbars(ax, sub[column], sub, fmt="o-", color=colours[modality], label=modality)
        style(ax, xlabel)
        ax.set_title(title, loc="left")
    axes[2].set_xscale("log")
    axes[2].legend(frameon=False, fontsize=8)

    # e: linkage error (sweep plus the reference point it passes through)
    data = results[results.sweep.isin(["linkage_error", "reference"])]
    for modality in MODALITIES:
        sub = data[data.modality == modality].sort_values("linkage_error_nm")
        errorbars(axes[4], sub.linkage_error_nm, sub, fmt="o-", color=colours[modality],
                  label=modality)
    style(axes[4], "linkage error (nm)")
    axes[4].set_title("e  linkage error", loc="left")

    # f: width of the blur (PSF standard deviation) of BLUR_MODALITY
    for sweep, column, marker, label in [
        ("lateral_blur", "lateral_resolution_nm", "o-", "lateral (axial %g nm)" % default_axial),
        ("axial_blur", "axial_resolution_nm", "s-", "axial (lateral %g nm)" % default_lateral),
    ]:
        sub = results[results.sweep == sweep].sort_values(column)
        errorbars(axes[5], sub[column], sub, fmt=marker, label=label)
    axes[5].set_xscale("log")
    style(axes[5], f"{BLUR_MODALITY} PSF width, s.d. (nm)")
    axes[5].set_title("f  width of the blur", loc="left")
    axes[5].legend(frameon=False, fontsize=8)

    # g: controls at the reference condition
    x = np.arange(len(MODALITIES))
    reference = reference.set_index("modality").loc[MODALITIES]
    for offset, prefix, label in [(-0.2, "", "shape (main)"), (0, "aspect_ratio_", "aspect ratio"),
                                  (0.2, "brightness_", "brightness only")]:
        errorbars(axes[6], x + offset, reference, prefix, fmt="o", label=label)
    axes[6].scatter(x, reference.null_auc, marker="x", color="k", label="null (dome vs dome)")
    axes[6].set_xticks(x, MODALITIES)
    style(axes[6], "")
    axes[6].set_title("g  controls", loc="left")
    axes[6].legend(frameon=False, fontsize=8)
    axes[7].axis("off")

    fig.tight_layout()
    fig.savefig(out_prefix + "distinguishability.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # labelling efficiency x structural integrity: one AUC heatmap per method
    les = sorted(grid.labelling_efficiency.unique())
    sis = sorted(grid.structural_integrity.unique())
    fig, axes = plt.subplots(1, len(MODALITIES), figsize=[4.2 * len(MODALITIES), 4.2],
                             sharey=True, constrained_layout=True)
    for ax, modality in zip(axes, MODALITIES):
        table = (grid[grid.modality == modality]
                 .pivot_table(index="structural_integrity", columns="labelling_efficiency",
                              values="auc")
                 .reindex(index=sis, columns=les))
        image = ax.imshow(table.values, origin="lower", cmap="viridis", vmin=0.5, vmax=1,
                          aspect="auto")
        for i in range(len(sis)):
            for j in range(len(les)):
                value = table.values[i, j]
                if np.isfinite(value):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7,
                            color="white" if value < 0.8 else "black")
        ax.set_xticks(range(len(les)), [f"{v:g}" for v in les], rotation=45)
        ax.set_yticks(range(len(sis)), [f"{v:g}" for v in sis])
        ax.set_xlabel("labelling efficiency")
        ax.set_title(modality)
    axes[0].set_ylabel("structural integrity")
    fig.colorbar(image, ax=axes, label="ROC AUC", shrink=0.9)
    fig.savefig(out_prefix + "distinguishability_heatmaps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Supp. Fig. S1: ventral geometry against uniform random orientations
    fig, ax = plt.subplots(figsize=[6, 4])
    for offset, sweep, label in [(-0.1, "reference", "ventral"),
                                 (0.1, "uniform_orientation", "uniform")]:
        sub = results[results.sweep == sweep].set_index("modality").loc[MODALITIES]
        errorbars(ax, x + offset, sub, fmt="o", label=label)
    ax.set_xticks(x, MODALITIES)
    style(ax, "")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_prefix + "distinguishability_orientation_S1.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    n_workers = max(1, N_WORKERS)
    if "--features" in sys.argv:
        # score an existing feature table again, without simulating
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        out_prefix = features_path.replace("features.csv", "")
        print(f"Scoring {features_path}: {len(features)} images, {n_workers} worker(s)")
    else:
        paths = {name: download(url) for name, url in MODELS.items()}
        from vlab4mic.experiments import ExperimentParametrisation, output_path

        out_prefix = os.path.join(
            str(output_path), ExperimentParametrisation().date_as_string + "fig2_clathrin_"
        )
        images_dir = None
        if SAVE_IMAGES:
            images_dir = out_prefix + "images"
            os.makedirs(images_dir, exist_ok=True)
        print(f"{len(list_sweeps()) * len(paths)} sweeps, {N_REALISATIONS} realisations per "
              f"lattice and condition, {n_workers} worker(s)"
              + (f", images saved to {images_dir}" if images_dir else ""))
        features = simulate(paths, n_workers, images_dir)
        features.to_csv(out_prefix + "features.csv", index=False)
        settings = {name: value for name, value in globals().items()
                    if name.isupper() and isinstance(value, (int, float, str, list, dict))}
        with open(out_prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)
        print(f"Feature table: {out_prefix}features.csv ({len(features)} images)")

    results, reference_scores = score_features(features, n_workers)
    results.to_csv(out_prefix + "distinguishability.csv", index=False)
    reference_blur = results[(results.sweep == "reference") & (results.modality == BLUR_MODALITY)]
    default_lateral = float(reference_blur.lateral_resolution_nm.iloc[0])
    default_axial = float(reference_blur.axial_resolution_nm.iloc[0])
    make_figures(results, reference_scores, default_lateral, default_axial, out_prefix)
    print(results.drop(columns=[c for c in results.columns if c.endswith(("_low", "_high"))]))


if __name__ == "__main__":
    main()
