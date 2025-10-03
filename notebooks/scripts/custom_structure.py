import numpy as np
import os
import matplotlib.pyplot as plt
from vlab4mic import experiments

np.random.seed(44)

# parameters
structure = "YOUR_CUSTOM_STRUCTURE_PATH"
modalities = ["STED", "SMLM",]
# Run simulation
images , noiseless_,  experiment = experiments.image_vsample(
    structure=structure,
    structure_is_path=True,
    clear_experiment=True,
    multimodal=modalities,
    run_simulation=True
)
plt.rcParams['figure.figsize'] = [20, 10]

fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(121, projection="3d")
experiment.particle.gen_axis_plot(axis_object=ax)
ax = fig.add_subplot(122)
ax.imshow(images["SMLM"][0], cmap="grey")
filename = os.path.join(experiment.output_directory, 'vlab4mic_custom_structure.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()