"""Supp. Fig. S7: HIV-1 capsid response to labelling efficiency and structural integrity.

Can an intact HIV-1 capsid be told apart from the same capsid with reduced
structural integrity, and how does that depend on labelling efficiency?

Same approach as Fig. 2 (fig2_clathrin_distinguishability.py) and the
Methods ("Repeated realisations and the distinguishability score"): every
repetition of run_parameter_sweep is an independent realisation, one capsid
with its own random labelling, structural integrity, orientation, position
and noise. The features of every image are computed by the sweep
(feature_metrics) and collected in one feature table with one row per
image. For each labelling efficiency and imaging method, the intact capsid
(structural integrity 1) is compared with each reduced integrity by
distinguishability_from_features: features standardised (and reduced by
PCA, 95% of the variance, as there are more than two), logistic regression
under stratified 5-fold cross-validation, out-of-fold ROC AUC with a 95%
bootstrap interval (1,000 resamples), a permutation p-value and a null
control (two random halves of the intact realisations).

Features: "shape" (radial intensity profile, rotation-averaged power
spectrum, radius of gyration, principal axes and their ratio), computed
after background subtraction and scaling to unit total intensity, so they
describe the outline and internal distribution of the labels and are blind
to brightness. A damaged capsid also carries fewer labels, so brightness
alone is scored as a control (brightness_auc), to show how much of the
difference is the number of labels rather than the shape.

The intact realisations at one labelling efficiency are shared by the
comparisons with the three reduced integrities.

Sample: one capsid (3J3Y) per realisation, anti-p24 primary antibody
template, uniformly random orientation and in-plane rotation, random
position within +-POSITION_JITTER_NM, in a FOV_NM x FOV_NM x 100 nm field,
imaged by STED, SMLM and Airyscan with their default parameters (Supp.
Table 3). Structural integrity clusters: 20 A (subunits) and 100 A
(neighbouring subunits).

Outputs, in the VLab4Mic output directory (<date>figS7_hiv_...):
- features.csv: one row per image;
- response.csv: one row per labelling efficiency, integrity and imaging
  method (AUC, interval, null control, p-value, brightness control);
- response_heatmaps.png/.pdf/.svg: AUC (top) and brightness-alone AUC
  (bottom) against labelling efficiency and integrity, one panel per method;
- examples.npz and examples.json: realisation 0 of every condition and
  method (image summed over frames) and its pixel size;
- examples.png/.pdf/.svg: gallery of those images at GALLERY_LE;
- settings.json.

Run:   python figS7_hiv_response.py
Score an existing feature table again, without simulating:
       python figS7_hiv_response.py --features <...features.csv>
(the gallery is redrawn from examples.npz if it is next to the table)

Environment variables: VLAB4MIC_N_REALISATIONS (100 per condition; e.g. 5
for a quick run), VLAB4MIC_N_PERMUTATIONS (200), VLAB4MIC_N_WORKERS (1;
sweeps and scoring in parallel processes). Each sweep parses 3J3Y, which
takes several GB of memory, and one sweep runs per labelling efficiency, so
at most len(LABELLING_EFFICIENCIES) sweeps run at once.
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
RANDOM_SEED = 24

STRUCTURE = "3J3Y"
PROBE = "anti-p24_primary_antibody_HIV"
MODALITIES = ["STED", "SMLM", "AiryScan"]
FEATURES = "shape"
CONTROLS = {"brightness": "brightness"}  # scored alone, as a control
FOV_NM = 600  # lateral size of the virtual sample and common field of view
POSITION_JITTER_NM = 100  # random particle position, +-100 nm in x and y

LABELLING_EFFICIENCIES = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0]
INTACT = 1.0
STRUCTURAL_INTEGRITIES = [0.25, 0.5, 0.75]  # each compared with INTACT
SI_CLUSTERS = dict(structural_integrity_small_cluster=[20],  # A
                   structural_integrity_large_cluster=[100])  # A
GALLERY_LE = [0.05, 0.25, 1.0]  # labelling efficiencies shown in the gallery
GALLERY_CROP_NM = 300  # gallery window, centred on the particle
SCALE_BAR_NM = 100

COMMON = dict(  # arguments of every sweep
    structures=[STRUCTURE], probe_templates=[PROBE], modalities=MODALITIES,
    sample_dimensions=[[FOV_NM, FOV_NM, 100]],
    random_orientations=True, random_rotations=True,
    position_jitter_nm=POSITION_JITTER_NM,
    structural_integrity=[INTACT] + STRUCTURAL_INTEGRITIES, **SI_CLUSTERS,
    sweep_repetitions=N_REALISATIONS, use_reference=False, default_metrics=[],
    custom_metrics=None, save_analysis_results=False, analysis_plots=False,
    save_sweep_images=False, return_generator=True,
)


def simulate(job):
    """
    One labelling efficiency: the feature table of every integrity and
    method, and realisation 0 of every condition as an example image.
    """
    labelling_efficiency, seed = job
    sweep = run_parameter_sweep(
        labelling_efficiency=[labelling_efficiency], random_seed=seed,
        **dict(COMMON, custom_metrics=feature_metrics([FEATURES, *CONTROLS.values()],
                                                      fov_nm=FOV_NM)),
    )
    table = sweep.analysis["dataframes"].copy()
    table["labelling_efficiency"] = labelling_efficiency
    examples = []
    first = table[table.Replica.astype(int) == 0]
    for _, row in first.iterrows():
        image = np.asarray(sweep.acquisition_outputs[row.Combination_id][0], dtype=float)
        if image.ndim == 3:
            image = image.sum(axis=0)
        pixelsize = sweep.acquisition_outputs_parameters[row.Combination_id][6]["pixelsize"]
        examples.append(dict(modality=row.modality_name, labelling_efficiency=labelling_efficiency,
                             structural_integrity=float(row.structural_integrity),
                             pixelsize_nm=float(pixelsize), image=image))
    return table, examples


def example_key(example):
    return (f"{example['modality']}_le{example['labelling_efficiency']:g}"
            f"_si{example['structural_integrity']:g}")


def save_examples(examples, prefix):
    np.savez_compressed(prefix + "examples.npz",
                        **{example_key(e): e["image"] for e in examples})
    with open(prefix + "examples.json", "w") as handle:
        json.dump([{k: v for k, v in e.items() if k != "image"} | {"key": example_key(e)}
                   for e in examples], handle, indent=1)


def load_examples(prefix):
    if not (os.path.exists(prefix + "examples.npz") and os.path.exists(prefix + "examples.json")):
        return None
    images = np.load(prefix + "examples.npz")
    with open(prefix + "examples.json") as handle:
        return [dict(e, image=images[e["key"]]) for e in json.load(handle)]


def score(features):
    """Intact against each reduced integrity, per labelling efficiency and method."""
    tables = []
    for integrity in STRUCTURAL_INTEGRITIES:
        pair = features[np.isclose(features.structural_integrity, INTACT)
                        | np.isclose(features.structural_integrity, integrity)].copy()
        pair["class"] = np.where(np.isclose(pair.structural_integrity, INTACT),
                                 "intact", "damaged")
        scores = distinguishability_from_features(
            pair, "class", ("intact", "damaged"), FEATURES,
            by=["modality_name", "labelling_efficiency"], controls=CONTROLS,
            replicate_column="Replica", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        scores.insert(2, "structural_integrity", integrity)
        tables.append(scores)
    return pd.concat(tables, ignore_index=True)


def save_figure(fig, path_without_extension):
    for extension in ("png", "pdf", "svg"):
        fig.savefig(f"{path_without_extension}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_heatmaps(scores, prefix):
    """AUC (shape features) and brightness-alone AUC, as in Fig. 2b."""
    data = scores.copy()
    data["labelling efficiency (%)"] = (100 * data.labelling_efficiency).round().astype(int)
    data["structural integrity"] = data.structural_integrity
    fig, axes = plt.subplots(2, len(MODALITIES), figsize=(3.4 * len(MODALITIES), 5.6),
                             sharey=True, constrained_layout=True)
    for row, (metric, title) in enumerate([("auc", "shape features"),
                                           ("brightness_auc", "brightness alone")]):
        plot_distinguishability_heatmaps(
            data, x="labelling efficiency (%)", y="structural integrity",
            facet="modality_name", facet_order=MODALITIES, metric=metric, axes=axes[row])
        for ax, modality in zip(axes[row], MODALITIES):
            ax.set_title(f"{modality.replace('AiryScan', 'Airyscan')}, {title}")
    save_figure(fig, prefix + "response_heatmaps")


def crop_on_particle(image, pixelsize_nm, size_nm=GALLERY_CROP_NM):
    """Square window of size_nm centred on the intensity centroid (above the median)."""
    n = max(int(round(size_nm / pixelsize_nm)), 1)
    if min(image.shape) <= n:
        return image
    signal = np.clip(image - np.median(image), 0, None)
    total = signal.sum()
    x, y = np.indices(image.shape)
    cx, cy = ((signal * x).sum() / total, (signal * y).sum() / total) if total > 0 else \
        (image.shape[0] / 2, image.shape[1] / 2)
    x0 = int(np.clip(round(cx - n / 2), 0, image.shape[0] - n))
    y0 = int(np.clip(round(cy - n / 2), 0, image.shape[1] - n))
    return image[x0:x0 + n, y0:y0 + n]


def make_gallery(examples, prefix):
    """Realisation 0 at GALLERY_LE: rows imaging methods, columns integrity."""
    integrities = [INTACT] + sorted(STRUCTURAL_INTEGRITIES, reverse=True)
    examples = [dict(e, image=crop_on_particle(e["image"], e["pixelsize_nm"])) for e in examples]
    lookup = {(e["modality"], round(e["labelling_efficiency"], 4),
               round(e["structural_integrity"], 4)): e for e in examples}
    n_cols = len(GALLERY_LE) * len(integrities)
    fig, axes = plt.subplots(len(MODALITIES), n_cols, figsize=(1.1 * n_cols, 1.25 * len(MODALITIES)),
                             squeeze=False)
    for r, modality in enumerate(MODALITIES):
        for b, le in enumerate(GALLERY_LE):
            block = [lookup.get((modality, round(le, 4), round(si, 4))) for si in integrities]
            # one intensity scale per method and labelling efficiency
            vmax = max((np.percentile(e["image"], 99.9) for e in block if e is not None),
                       default=1.0)
            for c, (si, example) in enumerate(zip(integrities, block)):
                ax = axes[r, b * len(integrities) + c]
                ax.set_axis_off()
                if example is None:
                    continue
                ax.imshow(example["image"], cmap="gray", vmin=0, vmax=vmax or 1.0)
                if r == 0:
                    ax.set_title(f"{int(round(100 * le))}% LE\nintegrity {si:g}", fontsize=7)
                if c == 0 and b == 0:
                    ax.text(-0.08, 0.5, modality.replace("AiryScan", "Airyscan"), fontsize=8,
                            rotation=90, ha="right", va="center", transform=ax.transAxes)
                if modality == "SMLM" and b == 0 and c == 0:
                    # every panel is GALLERY_CROP_NM wide: one scale bar
                    xlim, ylim = ax.get_xlim(), ax.get_ylim()
                    length_px = SCALE_BAR_NM / example["pixelsize_nm"]
                    height, width = example["image"].shape
                    ax.plot([width * 0.95 - length_px, width * 0.95], [height * 0.92] * 2,
                            color="white", linewidth=2, solid_capstyle="butt")
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
    fig.subplots_adjust(wspace=0.04, hspace=0.06)
    save_figure(fig, prefix + "examples")


def main():
    examples = None
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
        examples = load_examples(prefix)
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "figS7_hiv_")
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(LABELLING_EFFICIENCIES))
        jobs = [(le, int(seed)) for le, seed in zip(LABELLING_EFFICIENCIES, seeds)]
        print(f"{len(jobs)} sweeps, {N_REALISATIONS} realisations per condition, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(min(N_WORKERS, len(jobs))) as pool:
                outputs = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="sweeps"))
        else:
            outputs = [simulate(job) for job in tqdm(jobs, desc="sweeps")]
        features = pd.concat([table for table, _ in outputs], ignore_index=True)
        features.to_csv(prefix + "features.csv", index=False)
        examples = [example for _, block in outputs for example in block]
        save_examples(examples, prefix)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    scores = score(features)
    scores.to_csv(prefix + "response.csv", index=False)
    make_heatmaps(scores, prefix)
    if examples:
        make_gallery(examples, prefix)
    print(scores[["modality_name", "labelling_efficiency", "structural_integrity",
                  "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc", "p_value",
                  "brightness_auc"]].to_string())


if __name__ == "__main__":
    main()
