"""Fig. 2b-d and Supp. Fig. S1: distinguishability of domed and flat clathrin lattices.

For each imaging method, N independent realisations of each lattice (one
particle each, own labelling, orientation and noise) are simulated and the
distinguishability (accuracy and ROC AUC with 95% bootstrap intervals) is
computed with vlab4mic.analysis.distinguishability.

- Ventral geometry (Fig. 2b,c): lattices lie on the coverslip (axis +z)
  with a random tilt up to MAX_TILT_DEG and a random in-plane rotation.
- Labelling efficiency sweep (Fig. 2d) at the ventral geometry.
- Uniform random orientations (Supp. Fig. S1), as a control.

Set VLAB4MIC_N_REALISATIONS to change the number of realisations per class
(default 100; e.g. 5 for a quick run).
"""

import os
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

from vlab4mic import experiments
from vlab4mic.analysis.distinguishability import distinguishability_from_replicates

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 24
MAX_TILT_DEG = 10
MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]
LABELLING_EFFICIENCIES = [0.2, 0.4, 0.6, 0.8, 1.0]
LE_SWEEP_MODALITIES = ["STED", "SMLM"]
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
    # each lattice model is parsed once; geometry and labelling efficiency
    # are changed on the same experiment, and every realisation is imaged
    # in all modalities
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


def realisations(experiment, labelling_efficiency, geometry):
    probe_name = list(experiment.probe_parameters)[0]
    experiment.probe_parameters[probe_name]["labelling_efficiency"] = labelling_efficiency
    if geometry == "ventral":
        experiment.set_virtualsample_params(
            random_orientations=False, orientation_tilt_max=MAX_TILT_DEG,
            random_rotations=True,
        )
    else:
        experiment.virtualsample_params.pop("orientation_tilt_max", None)
        experiment.set_virtualsample_params(random_orientations=True, random_rotations=True)
    experiment.build(modules=["particle"])
    return experiment.run_replicates(N_REALISATIONS, modality="All")


def resolution(experiment, modality):
    psf = experiment.imaging_modalities[modality]["psf_params"]
    return [sd * v for sd, v in zip(psf["std_devs"], psf["voxelsize"])][::2]


def compare(labelling_efficiency, geometry, modalities):
    dome = realisations(experiments_by_model["dome"], labelling_efficiency, geometry)
    flat = realisations(experiments_by_model["flat"], labelling_efficiency, geometry)
    rows, scores = [], {}
    for modality in modalities:
        score = distinguishability_from_replicates(dome, flat, modality, random_state=RANDOM_SEED)
        lateral, axial = resolution(experiments_by_model["dome"], modality)
        rows.append(dict(
            modality=modality, geometry=geometry, labelling_efficiency=labelling_efficiency,
            lateral_resolution_nm=lateral, axial_resolution_nm=axial,
            n_per_class=N_REALISATIONS, accuracy=score["accuracy"],
            accuracy_low=score["accuracy_interval"][0], accuracy_high=score["accuracy_interval"][1],
            auc=score["auc"], auc_low=score["auc_interval"][0], auc_high=score["auc_interval"][1],
        ))
        scores[(modality, geometry)] = score
    return rows, scores


paths = {name: download(url) for name, url in MODELS.items()}
experiments_by_model = {
    name: build_experiment(path, RANDOM_SEED + k) for k, (name, path) in enumerate(paths.items())
}
experiment = experiments_by_model["dome"]
rows, scores = [], {}
for geometry in ["ventral", "uniform"]:
    new_rows, new_scores = compare(1.0, geometry, MODALITIES)
    rows += new_rows
    scores.update(new_scores)
for le in LABELLING_EFFICIENCIES[:-1]:
    new_rows, _ = compare(le, "ventral", LE_SWEEP_MODALITIES)
    rows += new_rows

results = pd.DataFrame(rows)
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "fig2_clathrin_")
results.to_csv(prefix + "distinguishability.csv", index=False)
experiment.save_parameters(out, name="fig2_clathrin_distinguishability")

# Fig. 2b: score distributions (ventral); Fig. 2c: AUC vs resolution;
# Fig. 2d: AUC vs labelling efficiency; Fig. S1: uniform orientations
fig, axes = plt.subplots(1, 3, figsize=[18, 5])
for k, modality in enumerate(MODALITIES):
    s = scores[(modality, "ventral")]
    axes[0].scatter(np.full(s["n_a"], k) - 0.1, s["scores"][s["labels"] == 0], s=6, label="dome" if k == 0 else None)
    axes[0].scatter(np.full(s["n_b"], k) + 0.1, s["scores"][s["labels"] == 1], s=6, label="flat" if k == 0 else None)
axes[0].set_xticks(range(len(MODALITIES)), MODALITIES)
axes[0].set_ylabel("P(flat), out-of-fold")
axes[0].legend()
full = results[results.labelling_efficiency == 1.0]
for geometry, marker in [("ventral", "o"), ("uniform", "x")]:
    sub = full[full.geometry == geometry]
    for res, colour in [("axial_resolution_nm", "C0"), ("lateral_resolution_nm", "C1")]:
        axes[1].errorbar(sub[res], sub.auc, yerr=[sub.auc - sub.auc_low, sub.auc_high - sub.auc],
                         fmt=marker, color=colour, label=f"{geometry}, {res.split('_')[0]}")
axes[1].set_xscale("log")
axes[1].set_xlabel("effective resolution (nm)")
axes[1].set_ylabel("AUC")
axes[1].legend(fontsize=7)
for modality in LE_SWEEP_MODALITIES:
    sub = results[(results.modality == modality) & (results.geometry == "ventral")].sort_values("labelling_efficiency")
    axes[2].errorbar(sub.labelling_efficiency, sub.auc, yerr=[sub.auc - sub.auc_low, sub.auc_high - sub.auc], fmt="o-", label=modality)
axes[2].set_xlabel("labelling efficiency")
axes[2].set_ylabel("AUC")
axes[2].legend()
fig.savefig(prefix + "distinguishability.png", dpi=300, bbox_inches="tight")
plt.close()
print(results)
