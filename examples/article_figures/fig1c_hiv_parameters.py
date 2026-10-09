"""Fig. 1c: parameter choices on the HIV-1 capsid.

HIV-1 capsid (3J3Y) labelled with the anti-p24 primary antibody, varying:
- labelling efficiency 100, 50 and 20% (structural integrity 100%, SMLM),
- structural integrity 100, 60 and 40% (labelling efficiency 100%, SMLM),
- imaging method STED, SMLM and Airyscan (both at 100%), with the default
  effective resolutions of Supp. Table S3.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from vlab4mic import experiments

random_seed = 1
LABELLING_EFFICIENCIES = [1.0, 0.5, 0.2]
STRUCTURAL_INTEGRITIES = [1.0, 0.6, 0.4]
MODALITIES = ["STED", "SMLM", "AiryScan"]


# the capsid is parsed once; only the labelled particle is rebuilt for
# each condition, and each realisation is imaged in all three modalities
_, _, experiment = experiments.image_vsample(
    structure="3J3Y",
    probe_template="anti-p24_primary_antibody_HIV",
    structural_integrity_small_cluster=20,
    structural_integrity_large_cluster=100,
    number_of_particles=1,
    sample_dimensions=[300, 300, 100],
    multimodal=MODALITIES,
    run_simulation=False,
    clear_experiment=True,
    random_seed=random_seed,
)
probe_name = list(experiment.probe_parameters)[0]


def simulate(labelling_efficiency, structural_integrity, modality):
    experiment.probe_parameters[probe_name]["labelling_efficiency"] = labelling_efficiency
    experiment.set_structural_integrity(structural_integrity=structural_integrity)
    experiment.build(modules=["particle"])
    result = experiment.run_replicates(1, modality="All")[0]
    return np.asarray(result["images"][modality]["ch0"])[0], experiment


panels = []
for le in LABELLING_EFFICIENCIES:
    panels.append((f"labelling efficiency {int(le * 100)}%", simulate(le, 1.0, "SMLM")))
for si in STRUCTURAL_INTEGRITIES:
    panels.append((f"structural integrity {int(si * 100)}%", simulate(1.0, si, "SMLM")))
for mod in MODALITIES:
    panels.append((mod, simulate(1.0, 1.0, mod)))

fig, axes = plt.subplots(3, 3, figsize=[12, 12])
for ax, (title, (image, _)) in zip(axes.ravel(), panels):
    ax.imshow(image, cmap="grey")
    ax.set_title(title)
    ax.set_axis_off()
experiment = panels[-1][1][1]
filename = os.path.join(
    experiment.output_directory,
    experiment.date_as_string + "vlab4mic_fig1c_hiv_parameters.png",
)
fig.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()
experiment.save_parameters(experiment.output_directory, name="fig1c_hiv_parameters")
