"""Supp. Fig. S4: resolving the three labelled sites of PCNA.

PCNA (4ZTD) labelled at SER-186 of each chain (site-specific, no linkage),
about 6 nm apart (Helmerich et al. 2024). For each localisation precision,
N independent realisations (one trimer seen from above, SMLM) give:
- the fraction of particles in which the three sites are resolved: a
  Gaussian mixture model of the localisations selects three components by
  BIC (particle_measures.sites_resolved);
- the recovered inter-site distance (k-means with three clusters on the
  localisations).
dSTORM (6.5 nm) and DNA-PAINT (3 nm) precisions of the published data are
marked.

Set VLAB4MIC_N_REALISATIONS to change the number of realisations.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans

from vlab4mic import experiments
from vlab4mic.analysis.particle_measures import sites_resolved

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 61
PRECISIONS_NM = [1, 2, 3, 4, 5, 6.5, 8]
PUBLISHED = {"DNA-PAINT": 3.0, "dSTORM": 6.5}
PIXEL_NM = 0.5
LOCALISATIONS_PER_SITE = 50
RENDERING_KERNEL_NM = 1.0
probe = dict(
    probe_template="Linker",
    probe_name="pcna_ser186",
    probe_target_type="Atom_residue",
    probe_target_value={"atoms": ["CA"], "residues": ["SER"], "position": 186},
    probe_distance_to_epitope=0,
)

rows = []
for i, precision in enumerate(PRECISIONS_NM):
    _, _, experiment = experiments.image_vsample(
        structure="4ZTD",
        probe_list=[probe],
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[100, 100, 50],
        sample_inital_orientation=[0, 0, 1],
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED + i,
    )
    experiment.set_virtualsample_params(random_orientations=False, random_rotations=True)
    experiment.update_modality(
        modality_name="SMLM", pixelsize_nm=PIXEL_NM, psf_voxel_nm=PIXEL_NM,
        simulate_localistations=True, lateral_precision=precision,
        axial_precision=precision, nlocalisations=LOCALISATIONS_PER_SITE,
        rendering_kernel_nm=RENDERING_KERNEL_NM,
    )
    for result in experiment.run_replicates(N_REALISATIONS, modality="SMLM"):
        record = list(result["positions"]["SMLM"]["ch0"].values())[0]
        locs = record["localisations"]
        resolved = False
        distance = np.nan
        if locs is not None and len(locs) >= 3:
            resolved = sites_resolved(locs, 3)[0]
            centres = KMeans(n_clusters=3, n_init=5, random_state=0).fit(locs[:, :2]).cluster_centers_
            distance = pdist(centres).mean()
        rows.append(dict(precision_nm=precision, resolved=resolved, inter_site_nm=distance))

results = pd.DataFrame(rows)
summary = results.groupby("precision_nm").agg(
    resolved_fraction=("resolved", "mean"),
    inter_site_mean=("inter_site_nm", "mean"),
    inter_site_sd=("inter_site_nm", "std"),
    n=("resolved", "size"),
)
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "figS4_pcna_")
results.to_csv(prefix + "particles.csv", index=False)
summary.to_csv(prefix + "summary.csv")
experiment.save_parameters(out, name="figS4_pcna_sites")

fig, ax = plt.subplots(figsize=[5, 4])
ax.plot(summary.index, summary.resolved_fraction, "o-")
for name, p in PUBLISHED.items():
    ax.axvline(p, ls="--", color="grey")
    ax.annotate(name, (p, 0.95))
ax.set_xlabel("localisation precision (nm)")
ax.set_ylabel("fraction with three resolved sites")
fig.savefig(prefix + "resolved.png", dpi=300, bbox_inches="tight")
plt.close()
print(summary)
