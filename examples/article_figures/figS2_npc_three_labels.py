"""Fig. 2e and Supp. Fig. S2c: nuclear pore ring radius for three labels.

Nup96 (7R5K, C-terminal motif ELAVGSL) labelled with a SNAP-tag, an
anti-GFP nanobody (GFP + nanobody, 6XZF) and an antibody, imaged by SMLM
seen from above (pore axis along z). For each label, N independent
realisations (one pore each) give the ring radius and width from a circle
fit to the localisations, compared with the published radii of the
reference dataset (Thevathasan et al. 2019): 53.7 +/- 2.1, 55.0 +/- 1.9
and 64.3 +/- 2.6 nm.

The antibody label is modelled as in the experiment: GFP fused to the
Nup96 C-terminus (primary) and an anti-GFP antibody (Antibody template)
bound to an epitope on the GFP barrel (secondary). Nup96-GFP is a genetic
knock-in, so every Nup96 carries GFP (labelling efficiency 1); the antibody
binds with the labelling efficiency of the other labels.

Only pores whose localisations cover at least MIN_CORNERS of the eight
corners (around the pore centre) are fitted, as sparse partial arcs give
unstable circle fits; the number of pores used is reported.

Set VLAB4MIC_N_REALISATIONS to change the number of realisations.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vlab4mic import experiments
from vlab4mic.analysis.particle_measures import (
    count_occupied_sectors,
    fit_circle,
    ring_measures,
)

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 44
LABELLING_EFFICIENCY = 0.4  # Supp. Note 4
MIN_CORNERS = 5
TARGET = dict(probe_target_type="Sequence", probe_target_value="ELAVGSL")
LABELS = {
    "SNAP-tag": [dict(probe_template="SNAP-tag", probe_name="SNAP-tag", **TARGET)],
    "anti-GFP nanobody": [dict(probe_template="GFP_w_nanobody", probe_name="GFP_w_nanobody", probe_DoL=2, **TARGET)],
    "anti-GFP antibody": [
        dict(probe_template="GFP", probe_name="Nup96-GFP", as_primary=True,
             labelling_efficiency=1.0, **TARGET),
        dict(probe_template="Antibody", probe_name="anti-GFP",
             probe_target_type="Primary", probe_target_value="Nup96-GFP"),
    ],
}
PUBLISHED_RADIUS_NM = {
    "SNAP-tag": (53.7, 2.1),
    "anti-GFP nanobody": (55.0, 1.9),
    "anti-GFP antibody": (64.3, 2.6),
}

rows = []
for i, (label, probes) in enumerate(LABELS.items()):
    _, _, experiment = experiments.image_vsample(
        structure="7R5K",
        probe_list=[dict({"labelling_efficiency": LABELLING_EFFICIENCY}, **p) for p in probes],
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
        points = record["localisations"]
        if points is None or len(points) < 3:
            rows.append(dict(label=label, used=False, n_localisations=0))
            continue
        # pore centre from the emitters (known structure) to count corners
        centre = fit_circle(record["emitters"])[0]
        corners = count_occupied_sectors(points, 8, centre)
        if corners < MIN_CORNERS:
            rows.append(dict(label=label, used=False, corners=corners,
                             n_localisations=len(points)))
            continue
        measures = ring_measures(points)
        rows.append(dict(label=label, used=True, corners=corners,
                         radius_nm=measures["radius"], width_nm=measures["width"],
                         n_localisations=len(points)))

results = pd.DataFrame(rows)
summary = results[results.used].groupby("label").agg(
    radius_mean=("radius_nm", "mean"), radius_sd=("radius_nm", "std"),
    width_mean=("width_nm", "mean"), n_used=("radius_nm", "size"),
).reindex(list(LABELS))
summary["n_total"] = results.groupby("label").size().reindex(list(LABELS))
summary["published_radius"] = [PUBLISHED_RADIUS_NM[k][0] for k in summary.index]
summary["published_sd"] = [PUBLISHED_RADIUS_NM[k][1] for k in summary.index]
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "figS2_npc_three_labels_")
results.to_csv(prefix + "particles.csv", index=False)
summary.to_csv(prefix + "summary.csv")
experiment.save_parameters(out, name="figS2_npc_three_labels")

fig, ax = plt.subplots(figsize=[5, 5])
ax.errorbar(summary.published_radius, summary.radius_mean, xerr=summary.published_sd,
            yerr=summary.radius_sd, fmt="o")
for label, row in summary.iterrows():
    ax.annotate(label, (row.published_radius, row.radius_mean))
lims = [45, 75]
ax.plot(lims, lims, "--", color="grey")
ax.set_xlabel("published ring radius (nm)")
ax.set_ylabel("predicted ring radius (nm)")
fig.savefig(prefix + "radius.png", dpi=300, bbox_inches="tight")
plt.close()
print(summary)
