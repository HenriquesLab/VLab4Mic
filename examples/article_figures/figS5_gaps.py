"""Supp. Fig. S5: apparent breaks and resolved corners in the nuclear pore.

The nuclear pore has eight Nup96 corners, 45 degrees apart (cryo-EM). For
each labelling efficiency, N independent realisations (SMLM, one pore seen
from above) give:
- the fraction of particles with an apparent break: an angular gap between
  localisations wider than BREAK_GAP_DEG (more than one corner spacing,
  i.e. at least one unlabelled corner);
- the corner detection rate: occupied 45-degree sectors / 8;
- the fraction of particles with all eight corners resolved.

Set VLAB4MIC_N_REALISATIONS to change the number of realisations.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd

from vlab4mic import experiments
from vlab4mic.analysis.particle_measures import (
    count_occupied_sectors,
    fit_circle,
    has_apparent_break,
)

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 51
LABELLING_EFFICIENCIES = [0.2, 0.4, 0.6, 0.8, 1.0]
BREAK_GAP_DEG = 67.5  # 1.5 corner spacings

rows = []
for i, le in enumerate(LABELLING_EFFICIENCIES):
    _, _, experiment = experiments.image_vsample(
        structure="7R5K",
        probe_template="NPC_Nup96_Cterminal_direct",
        labelling_efficiency=le,
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[400, 400, 200],
        sample_inital_orientation=[0, 0, 1],
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED + i,
    )
    experiment.set_virtualsample_params(random_orientations=False, random_rotations=True)
    for result in experiment.run_replicates(N_REALISATIONS, modality="SMLM"):
        record = list(result["positions"]["SMLM"]["ch0"].values())[0]
        emitters = record["emitters"]
        points = record["localisations"] if record["localisations"] is not None else emitters
        if len(points) < 3:
            rows.append(dict(labelling_efficiency=le, apparent_break=True, corners=0))
            continue
        # centre from the emitters' ring (known structure), not the sparse data
        centre = fit_circle(emitters)[0] if len(emitters) >= 3 else None
        rows.append(dict(
            labelling_efficiency=le,
            apparent_break=has_apparent_break(points, BREAK_GAP_DEG, centre),
            corners=count_occupied_sectors(points, 8, centre),
        ))

results = pd.DataFrame(rows)
summary = results.groupby("labelling_efficiency").agg(
    break_fraction=("apparent_break", "mean"),
    corner_detection_rate=("corners", lambda c: c.mean() / 8),
    all_corners_fraction=("corners", lambda c: (c == 8).mean()),
    n=("corners", "size"),
)
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "figS5_gaps_")
results.to_csv(prefix + "particles.csv", index=False)
summary.to_csv(prefix + "summary.csv")
experiment.save_parameters(out, name="figS5_gaps")

fig, ax = plt.subplots(figsize=[5, 4])
for column in ["break_fraction", "corner_detection_rate", "all_corners_fraction"]:
    ax.plot(summary.index, summary[column], "o-", label=column.replace("_", " "))
ax.set_xlabel("labelling efficiency")
ax.set_ylabel("fraction")
ax.legend()
fig.savefig(prefix + "summary.png", dpi=300, bbox_inches="tight")
plt.close()
print(summary)
