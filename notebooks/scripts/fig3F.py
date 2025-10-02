import numpy as np
import os
import matplotlib.pyplot as plt
from vlab4mic import experiments

np.random.seed(44)

# parameters
structure = "7R5K"
probe_template = "NPC_Nup96_Cterminal_direct"
modalities = ["Widefield_Thev2016", "Confocal_Thev2016", "AiryScan_Thev2016", "STED_Thev2016", "STORM_Thev2016"]

# Run simulation
images , noiseless_,  experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    clear_experiment=True,
    random_rotations=True,
    rotation_angles=None,
    number_of_particles = 10,
    multimodal=modalities,
    run_simulation=True
)
plt.rcParams['figure.figsize'] = [20, 10]
nmods = len(modalities)

fig, axs = plt.subplots(1, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][nframe], cmap="magma")
    axs[i].set_title(mod)

filename = os.path.join(experiment.output_directory, 'vlab4mic_fig3F.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()