"""Fig. 2 and Supp. Fig. S1: distinguishability of domed and flat clathrin lattices.

Structural distinguishability is the accuracy with which a classifier tells
two candidate structures apart under stated simulation parameters. For each
condition, N independent realisations of each lattice are simulated (one
particle each, with its own random labelling pattern, orientation and
noise), and vlab4mic.analysis.distinguishability reports the out-of-fold
accuracy and ROC AUC of a classifier, with 95% bootstrap intervals.

Classifier: each image is background subtracted and scaled to unit total
intensity, then described by rotation-invariant shape features (radial
profile, rotation-averaged power spectrum, principal axes); standardisation,
PCA and logistic regression, 5-fold stratified cross-validation. Because the
images are normalised, the lattices are not told apart by brightness (the
dome carries more labelled sites than the flat lattice). Every score comes
with four controls:

- brightness: the AUC from total intensity alone, the difference the shape
  features deliberately ignore;
- aspect_ratio: the AUC from the size-invariant aspect ratio alone;
- null_auc: two random halves of the domes, which should score 0.5;
- p_value: permutation p-value of the shape AUC.

Geometry (Leterrier, co-author): clathrin pits and plaques are imaged at the
ventral surface of cultured cells, lying flat on the coverslip, so the two
forms differ mainly in depth. The lattices lie on the coverslip (axis +z)
with a random tilt up to MAX_TILT_DEG and a random in-plane rotation. Uniform
random orientations are a supplementary control.

Conditions (ventral geometry unless stated):

- reference: every imaging method at its default effective resolution,
  labelling efficiency 1 and the antibody's linkage error (Fig. 2b,c);
- labelling efficiency sweep, every imaging method (Fig. 2d);
- linkage error sweep (distance from fluorophore to epitope), every
  imaging method (Fig. 2e);
- width of the microscope's blur: lateral and axial PSF width of a
  PSF-based method (BLUR_MODALITY) swept separately (Fig. 2f);
- uniform random orientations, every imaging method (Supp. Fig. S1).

Environment variables (for a quick run, e.g. VLAB4MIC_N_REALISATIONS=5
VLAB4MIC_N_PERMUTATIONS=0):

- VLAB4MIC_N_REALISATIONS: realisations per lattice and condition (100);
- VLAB4MIC_N_PERMUTATIONS: permutations for the p-value (200);
- VLAB4MIC_MAX_TILT_DEG: maximum tilt on the coverslip in degrees (10);
- VLAB4MIC_N_WORKERS: conditions simulated in parallel, one process each
  (1). Each worker loads both lattice models (a few GB of memory). Every
  condition has its own seed, so results do not depend on this number.
"""

import contextlib
import os
import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tqdm.auto import tqdm

from vlab4mic import experiments
from vlab4mic.analysis.distinguishability import distinguishability_from_replicates

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
# parallel worker processes, one condition each at a time (1: serial)
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
RANDOM_SEED = 24
# maximum random tilt of the lattices on the coverslip (degrees); no
# published reference for this value
MAX_TILT_DEG = float(os.environ.get("VLAB4MIC_MAX_TILT_DEG", 10))
MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

# reference condition, and the values swept one at a time around it
LABELLING_EFFICIENCY = 1.0
LABELLING_EFFICIENCIES = [0.2, 0.4, 0.6, 0.8, 1.0]
LINKAGE_ERROR_NM = 10  # distance to epitope of the Antibody probe template
LINKAGE_ERRORS_NM = [2, 5, 10, 20]  # direct label / nanobody to primary+secondary
BLUR_MODALITY = "STED"  # PSF-based method whose PSF width is swept
LATERAL_BLUR_NM = [10, 20, 40, 80, 130]  # PSF standard deviation, axial at default
AXIAL_BLUR_NM = [50, 100, 200, 300, 600]  # PSF standard deviation, lateral at default

# features scored for every condition: the main classifier and controls
MAIN_FEATURES = "shape"
CONTROL_FEATURES = ["aspect_ratio", "brightness"]

MODELS = {
    "dome": "https://zenodo.org/records/20719843/files/Dome_model_v2.0.cif",
    "flat": "https://zenodo.org/records/20719843/files/Flat_lattice_model_v2.0.cif",
}
PROBE = dict(
    probe_template="Antibody",
    probe_name="anti_CHC",
    probe_target_type="Sequence",
    probe_target_value="ATETQ",
    probe_DoL=4,
)


def download(url):
    path = Path(tempfile.gettempdir()) / Path(url).name
    if not path.exists():
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path.write_bytes(response.content)
    return str(path)


def build_experiment(model_path, seed):
    # each lattice model is parsed once; probe, geometry and PSF are changed
    # on the same experiment. Realisations continue the experiment's random
    # stream, so every realisation of every condition is independent.
    _, _, experiment = experiments.image_vsample(
        structure=model_path,
        structure_is_path=True,
        probe_list=[dict(PROBE)],
        multimodal=MODALITIES,
        number_of_particles=1,
        sample_dimensions=[600, 600, 100],
        structure_global_normal_orientation="local_plane",
        run_simulation=False,
        clear_experiment=True,
        random_seed=seed,
    )
    return experiment


def set_probe(experiment, labelling_efficiency, linkage_error_nm):
    experiment.remove_probes()
    experiment.add_probe(
        **PROBE,
        labelling_efficiency=labelling_efficiency,
        probe_distance_to_epitope=linkage_error_nm,
    )


def set_geometry(experiment, geometry):
    if geometry == "ventral":
        experiment.set_virtualsample_params(
            random_orientations=False, orientation_tilt_max=MAX_TILT_DEG,
            random_rotations=True,
        )
    else:
        experiment.virtualsample_params.pop("orientation_tilt_max", None)
        experiment.set_virtualsample_params(random_orientations=True, random_rotations=True)


def resolution(experiment, modality):
    # effective lateral and axial resolution (nm): localisation precision
    # for localisation-based methods, PSF standard deviation otherwise
    emitters = experiment.imaging_modalities[modality].get("emitters") or {}
    if emitters.get("lateral_precision") is not None:
        return emitters["lateral_precision"], emitters["axial_precision"]
    psf = experiment.imaging_modalities[modality]["psf_params"]
    lateral, _, axial = [sd * v for sd, v in zip(psf["std_devs"], psf["voxelsize"])]
    return lateral, axial


def set_blur(modality, lateral_nm, axial_nm):
    for experiment in experiments_by_model.values():
        experiment.update_modality(
            modality, lateral_resolution_nm=lateral_nm, axial_resolution_nm=axial_nm
        )


def realisations(experiment, geometry, labelling_efficiency, linkage_error_nm, modality):
    set_probe(experiment, labelling_efficiency, linkage_error_nm)
    set_geometry(experiment, geometry)
    experiment.build(modules=["particle"])
    return experiment.run_replicates(N_REALISATIONS, modality=modality)


def interval_columns(prefix, score):
    return {
        f"{prefix}accuracy": score["accuracy"],
        f"{prefix}accuracy_low": score["accuracy_interval"][0],
        f"{prefix}accuracy_high": score["accuracy_interval"][1],
        f"{prefix}auc": score["auc"],
        f"{prefix}auc_low": score["auc_interval"][0],
        f"{prefix}auc_high": score["auc_interval"][1],
    }


def compare(sweep, geometry="ventral", labelling_efficiency=LABELLING_EFFICIENCY,
            linkage_error_nm=LINKAGE_ERROR_NM, modalities=MODALITIES):
    # simulate both lattices under one condition and score every modality
    modality = modalities[0] if len(modalities) == 1 else "All"
    replicates = {
        name: realisations(experiment, geometry, labelling_efficiency,
                           linkage_error_nm, modality)
        for name, experiment in experiments_by_model.items()
    }
    rows, scores = [], {}
    for mod in modalities:
        lateral, axial = resolution(experiments_by_model["dome"], mod)
        row = dict(
            sweep=sweep, modality=mod, geometry=geometry,
            labelling_efficiency=labelling_efficiency,
            linkage_error_nm=linkage_error_nm,
            lateral_resolution_nm=lateral, axial_resolution_nm=axial,
            n_per_class=N_REALISATIONS, features=MAIN_FEATURES,
        )
        score = distinguishability_from_replicates(
            replicates["dome"], replicates["flat"], mod, features=MAIN_FEATURES,
            random_state=RANDOM_SEED, n_permutations=N_PERMUTATIONS,
        )
        row.update(interval_columns("", score))
        row.update(null_auc=score["null_auc"], p_value=score["p_value"])
        for control in CONTROL_FEATURES:
            control_score = distinguishability_from_replicates(
                replicates["dome"], replicates["flat"], mod, features=control,
                random_state=RANDOM_SEED, null_control=False,
            )
            row.update(interval_columns(f"{control}_", control_score))
        rows.append(row)
        scores[mod] = score
    return rows, scores


# conditions and parallel execution -----------------------------------------
experiments_by_model = {}  # set in every process by setup_experiments


def setup_experiments(paths):
    # parse both lattice models once per process (worker initializer)
    experiments_by_model.update({
        name: build_experiment(path, RANDOM_SEED + k) for k, (name, path) in enumerate(paths.items())
    })


def list_conditions(default_lateral, default_axial):
    # every condition: (label, compare arguments, (lateral, axial) blur of
    # BLUR_MODALITY). The slow all-modality conditions come first.
    default = (default_lateral, default_axial)
    conditions = [("reference", dict(sweep="reference"), default),
                  ("uniform orientation", dict(sweep="uniform_orientation", geometry="uniform"),
                   default)]
    conditions += [(f"labelling efficiency {le}",
                    dict(sweep="labelling_efficiency", labelling_efficiency=le), default)
                   for le in LABELLING_EFFICIENCIES if le != LABELLING_EFFICIENCY]
    conditions += [(f"linkage error {d} nm", dict(sweep="linkage_error", linkage_error_nm=d),
                    default) for d in LINKAGE_ERRORS_NM if d != LINKAGE_ERROR_NM]
    conditions += [(f"{BLUR_MODALITY} lateral blur {v} nm",
                    dict(sweep="lateral_blur", modalities=[BLUR_MODALITY]), (v, default_axial))
                   for v in LATERAL_BLUR_NM]
    conditions += [(f"{BLUR_MODALITY} axial blur {v} nm",
                    dict(sweep="axial_blur", modalities=[BLUR_MODALITY]), (default_lateral, v))
                   for v in AXIAL_BLUR_NM]
    return conditions


def condition_seeds(n_conditions):
    # three seeds per condition (dome field, flat field, numpy global), so
    # each condition is reproducible whatever the number of workers and the
    # order in which conditions run
    state = np.random.SeedSequence(RANDOM_SEED).generate_state(3 * n_conditions)
    return [tuple(int(v) for v in state[3 * k:3 * k + 3]) for k in range(n_conditions)]


def reseed(seeds):
    # field seeds (placement, orientation) come from each experiment's seed
    # generator; labelling and photon/detector noise use numpy's global one
    for experiment, seed in zip(experiments_by_model.values(), seeds[:2]):
        experiment.random_seed = seed
        experiment._seed_generator = None
    np.random.seed(seeds[2])


@contextlib.contextmanager
def quiet():
    # hide the per-simulation messages and warnings of the library, so the
    # progress bar stays readable; errors are still raised
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull), \
            warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def setup_experiments_quietly(paths):
    with quiet():
        setup_experiments(paths)


def run_condition(index, arguments, blur, seeds):
    with quiet():
        reseed(seeds)
        set_blur(BLUR_MODALITY, *blur)
        rows, scores = compare(**arguments)
    return index, rows, scores


def run_conditions(conditions, paths, n_workers):
    seeds = condition_seeds(len(conditions))
    jobs = [(k, arguments, blur, seeds[k]) for k, (_, arguments, blur) in enumerate(conditions)]
    results = [None] * len(conditions)
    progress = tqdm(total=len(conditions), desc="conditions", unit="condition")
    if n_workers <= 1:
        for job, (label, _, _) in zip(jobs, conditions):
            progress.set_postfix_str(label)
            index, rows, scores = run_condition(*job)
            results[index] = (rows, scores)
            progress.update()
    else:
        from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

        # one thread per worker for numpy / scipy, so workers do not compete
        for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.setdefault(variable, "1")
        with ProcessPoolExecutor(n_workers, initializer=setup_experiments_quietly,
                                 initargs=(paths,)) as pool:
            futures = {pool.submit(run_condition, *job): label
                       for job, (label, _, _) in zip(jobs, conditions)}
            pending = set(futures)
            while pending:
                # wake up regularly so the elapsed time keeps moving while
                # the long conditions run
                done, pending = wait(pending, timeout=5, return_when=FIRST_COMPLETED)
                for future in done:
                    index, rows, scores = future.result()
                    results[index] = (rows, scores)
                    progress.set_postfix_str(f"done: {futures[future]}")
                    progress.update()
                progress.refresh()
    progress.close()
    rows = [row for condition_rows, _ in results for row in condition_rows]
    reference_scores = results[0][1]
    return pd.DataFrame(rows), reference_scores


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


def make_figures(results, reference_scores, default_lateral, default_axial, out_prefix):
    colours = {m: f"C{k}" for k, m in enumerate(MODALITIES)}
    fig, axes = plt.subplots(2, 3, figsize=[16, 9])
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

    # c, d: labelling efficiency and linkage error (each sweep plus the
    # reference point it passes through), one line per method
    for ax, sweep, column, xlabel, title in [
        (axes[2], "labelling_efficiency", "labelling_efficiency", "labelling efficiency",
         "c  labelling efficiency"),
        (axes[3], "linkage_error", "linkage_error_nm", "linkage error (nm)", "d  linkage error"),
    ]:
        data = results[results.sweep.isin([sweep, "reference"])]
        for modality in MODALITIES:
            sub = data[data.modality == modality].sort_values(column)
            errorbars(ax, sub[column], sub, fmt="o-", color=colours[modality], label=modality)
        style(ax, xlabel)
        ax.set_title(title, loc="left")
    axes[2].legend(frameon=False, fontsize=8)

    # e: width of the blur (PSF standard deviation) of BLUR_MODALITY
    for sweep, column, marker, label in [
        ("lateral_blur", "lateral_resolution_nm", "o-", "lateral (axial %g nm)" % default_axial),
        ("axial_blur", "axial_resolution_nm", "s-", "axial (lateral %g nm)" % default_lateral),
    ]:
        sub = results[results.sweep == sweep].sort_values(column)
        errorbars(axes[4], sub[column], sub, fmt=marker, label=label)
    axes[4].set_xscale("log")
    style(axes[4], f"{BLUR_MODALITY} PSF width, s.d. (nm)")
    axes[4].set_title("e  width of the blur", loc="left")
    axes[4].legend(frameon=False, fontsize=8)

    # f: controls at the reference condition
    x = np.arange(len(MODALITIES))
    reference = reference.set_index("modality").loc[MODALITIES]
    for offset, prefix, label in [(-0.2, "", "shape (main)"), (0, "aspect_ratio_", "aspect ratio"),
                                  (0.2, "brightness_", "brightness only")]:
        errorbars(axes[5], x + offset, reference, prefix, fmt="o", label=label)
    axes[5].scatter(x, reference.null_auc, marker="x", color="k", label="null (dome vs dome)")
    axes[5].set_xticks(x, MODALITIES)
    style(axes[5], "")
    axes[5].set_title("f  controls", loc="left")
    axes[5].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_prefix + "distinguishability.png", dpi=300, bbox_inches="tight")
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
    paths = {name: download(url) for name, url in MODELS.items()}
    # the main process also builds the experiments: defaults, output
    # directory and parameter file, and the conditions when N_WORKERS is 1
    setup_experiments_quietly(paths)
    experiment = experiments_by_model["dome"]
    default_lateral, default_axial = resolution(experiment, BLUR_MODALITY)
    conditions = list_conditions(default_lateral, default_axial)
    n_workers = max(1, min(N_WORKERS, len(conditions)))
    print(f"{len(conditions)} conditions, {N_REALISATIONS} realisations per lattice, "
          f"{n_workers} worker(s)")
    results, reference_scores = run_conditions(conditions, paths, n_workers)

    out = experiment.output_directory
    out_prefix = os.path.join(out, experiment.date_as_string + "fig2_clathrin_")
    results.to_csv(out_prefix + "distinguishability.csv", index=False)
    experiment.save_parameters(out, name="fig2_clathrin_distinguishability")
    make_figures(results, reference_scores, default_lateral, default_axial, out_prefix)
    print(results.drop(columns=[c for c in results.columns if c.endswith(("_low", "_high"))]))


if __name__ == "__main__":
    main()
