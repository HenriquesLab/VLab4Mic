"""Supp. Fig. S7: HIV-1 capsid response curves.

Distinguishability (ROC AUC with 95% bootstrap intervals) of the intact
HIV-1 capsid (3J3Y, structural integrity 1) from the same capsid with
reduced structural integrity, against labelling efficiency and structural
integrity, for STED, SMLM and Airyscan. Each condition uses N independent
realisations per class (one capsid each, anti-p24 primary antibody,
uniform random orientations). Structural integrity clusters: small 20,
large 100 (Supp. Note 4).

Set VLAB4MIC_N_REALISATIONS to change the number of realisations.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd

from vlab4mic import experiments
from vlab4mic.analysis.distinguishability import distinguishability_from_replicates

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 24
MODALITIES = ["STED", "SMLM", "AiryScan"]
LABELLING_EFFICIENCIES = [0.25, 0.5, 0.75, 1.0]
STRUCTURAL_INTEGRITIES = [0.25, 0.5, 0.75]


# The capsid is parsed once (the most expensive and memory-hungry step);
# labelling efficiency and structural integrity are changed on the same
# experiment and only the labelled particle is rebuilt. Every realisation
# is imaged in all modalities.
_, _, experiment = experiments.image_vsample(
    structure="3J3Y",
    probe_template="anti-p24_primary_antibody_HIV",
    structural_integrity_small_cluster=20,
    structural_integrity_large_cluster=100,
    multimodal=MODALITIES,
    number_of_particles=1,
    sample_dimensions=[300, 300, 100],
    random_orientations=True,
    run_simulation=False,
    clear_experiment=True,
    random_seed=RANDOM_SEED,
)
probe_name = list(experiment.probe_parameters)[0]


def realisations(labelling_efficiency, structural_integrity):
    experiment.probe_parameters[probe_name]["labelling_efficiency"] = labelling_efficiency
    experiment.set_structural_integrity(structural_integrity=structural_integrity)
    experiment.build(modules=["particle"])
    return experiment.run_replicates(N_REALISATIONS, modality="All")


rows = []
for le in LABELLING_EFFICIENCIES:
    intact = realisations(le, 1.0)
    for si in STRUCTURAL_INTEGRITIES:
        damaged = realisations(le, si)
        for modality in MODALITIES:
            score = distinguishability_from_replicates(intact, damaged, modality, random_state=RANDOM_SEED)
            rows.append(dict(
                modality=modality, labelling_efficiency=le, structural_integrity=si,
                n_per_class=N_REALISATIONS, auc=score["auc"],
                auc_low=score["auc_interval"][0], auc_high=score["auc_interval"][1],
                accuracy=score["accuracy"],
            ))

results = pd.DataFrame(rows)
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "figS7_hiv_")
results.to_csv(prefix + "response.csv", index=False)
experiment.save_parameters(out, name="figS7_hiv_response")

fig, axes = plt.subplots(1, len(MODALITIES), figsize=[5 * len(MODALITIES), 4], sharey=True)
for ax, modality in zip(axes, MODALITIES):
    sub = results[results.modality == modality]
    for si in STRUCTURAL_INTEGRITIES:
        s = sub[sub.structural_integrity == si].sort_values("labelling_efficiency")
        ax.errorbar(s.labelling_efficiency, s.auc, yerr=[s.auc - s.auc_low, s.auc_high - s.auc],
                    fmt="o-", label=f"integrity {si}")
    ax.set_title(modality)
    ax.set_xlabel("labelling efficiency")
axes[0].set_ylabel("AUC, intact vs damaged")
axes[0].legend()
fig.savefig(prefix + "response.png", dpi=300, bbox_inches="tight")
plt.close()
print(results)
