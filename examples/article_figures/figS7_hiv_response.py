"""Supp. Fig. S7: HIV-1 capsid response curves.

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
features.csv (one row per image), response.csv (one row per labelling
efficiency, integrity and imaging method), settings.json and response.png.

Run:   python figS7_hiv_response.py
Score an existing feature table again, without simulating:
       python figS7_hiv_response.py --features <...features.csv>

Environment variables: VLAB4MIC_N_REALISATIONS (100 per condition; e.g. 5
for a quick run), VLAB4MIC_N_PERMUTATIONS (200), VLAB4MIC_N_WORKERS (1;
sweeps and scoring in parallel processes). Parsing 3J3Y takes several GB of
memory per sweep, so keep the number of workers within the memory available.
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

LABELLING_EFFICIENCIES = [0.25, 0.5, 0.75, 1.0]
INTACT = 1.0
STRUCTURAL_INTEGRITIES = [0.25, 0.5, 0.75]  # each compared with INTACT
SI_CLUSTERS = dict(structural_integrity_small_cluster=[20],  # A
                   structural_integrity_large_cluster=[100])  # A

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
CONDITIONS = ["modality_name", "labelling_efficiency", "structural_integrity"]


def simulate(job):
    """One labelling efficiency: feature table of every integrity and method."""
    labelling_efficiency, seed = job
    sweep = run_parameter_sweep(
        labelling_efficiency=[labelling_efficiency], random_seed=seed,
        **dict(COMMON, custom_metrics=feature_metrics([FEATURES, *CONTROLS.values()],
                                                      fov_nm=FOV_NM)),
    )
    table = sweep.analysis["dataframes"].copy()
    table["labelling_efficiency"] = labelling_efficiency
    return table


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


def make_figure(scores, prefix):
    fig, axes = plt.subplots(1, len(MODALITIES), figsize=[5 * len(MODALITIES), 4], sharey=True)
    for ax, modality in zip(axes, MODALITIES):
        sub = scores[scores.modality_name == modality]
        for k, integrity in enumerate(STRUCTURAL_INTEGRITIES):
            s = sub[np.isclose(sub.structural_integrity, integrity)].sort_values("labelling_efficiency")
            ax.errorbar(s.labelling_efficiency, s.auc, yerr=[s.auc - s.auc_low, s.auc_high - s.auc],
                        fmt="o-", capsize=2, color=f"C{k}", label=f"integrity {integrity:g}")
            ax.plot(s.labelling_efficiency, s.brightness_auc, ":", color=f"C{k}", alpha=0.6)
        ax.axhline(0.5, color="0.5", linestyle="--", linewidth=1)
        ax.set_ylim(0.3, 1.05)
        ax.set_title(modality)
        ax.set_xlabel("labelling efficiency")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("ROC AUC, intact vs reduced integrity")
    axes[0].legend(frameon=False, title="shape features (dotted: brightness alone)")
    fig.tight_layout()
    fig.savefig(prefix + "response.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "figS7_hiv_")
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(LABELLING_EFFICIENCIES))
        jobs = [(le, int(seed)) for le, seed in zip(LABELLING_EFFICIENCIES, seeds)]
        print(f"{len(jobs)} sweeps, {N_REALISATIONS} realisations per condition, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(min(N_WORKERS, len(jobs))) as pool:
                tables = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="sweeps"))
        else:
            tables = [simulate(job) for job in tqdm(jobs, desc="sweeps")]
        features = pd.concat(tables, ignore_index=True)
        features.to_csv(prefix + "features.csv", index=False)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    scores = score(features)
    scores.to_csv(prefix + "response.csv", index=False)
    make_figure(scores, prefix)
    print(scores[["modality_name", "labelling_efficiency", "structural_integrity",
                  "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc", "p_value",
                  "brightness_auc"]].to_string())


if __name__ == "__main__":
    main()
