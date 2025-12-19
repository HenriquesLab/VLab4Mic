from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt
from IPython.utils import io
import os
import numpy as np
np.random.seed(44)

modalities = ["STED", "SMLM",]
target_colour="#01579D"


############  NPC 
image_outputs1, image_outputs_noiseless1, experiment1 = image_vsample(
    structure = "7R5K",
    probe_template = "GFP_w_nanobody",
    probe_target_type = "Sequence",
    probe_target_value = "ELAVGSL",
    multimodal=modalities,
    clear_experiment=True
)
############  T4 capsid head
image_outputs2, image_outputs_noiseless2, experiment2 = image_vsample(
    structure = "8GMO",
    probe_template = "NHS_ester",
    labelling_efficiency = 0.01,
    multimodal=modalities,
    random_orientations=True,
    yz_orientations=[90,],
    clear_experiment=True
)
############  HIV-capsid
image_outputs3, image_outputs_noiseless3, experiment3 = image_vsample(
    structure = "3J3Y",
    probe_template = "anti-p24_primary_antibody_HIV",
    multimodal=modalities,
    clear_experiment=True
)
############ Clathrin coated pit, primary and secodnary labelling

primary = dict(
    probe_template = "Antibody",
    probe_name="Primary-clathrin",
    probe_target_type = "Sequence",
    probe_target_value = "EQATETQ",
)
secondary = dict(
    probe_template = "Antibody",
    probe_name="Secondary-clathrin",
    probe_target_type = "Primary",
    probe_target_value = "Primary-clathrin"
)
image_outputs4, image_outputs_noiseless4, experiment4 = image_vsample(
    structure = "1XI5",
    primary_probe = primary,
    secondary_probe = secondary,
    multimodal=modalities,
    clear_experiment=True,
    run_simulation=True
)


############ save images as figures
#######
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(121, projection="3d")
experiment1.particle.gen_axis_plot(with_sources=True, axis_object=ax,target_colour=target_colour, source_plotsize=5, emitter_plotsize=5)
ax = fig.add_subplot(122)
ax.imshow(image_outputs1["SMLM"]["ch0"][0], cmap="grey")
filename = os.path.join(experiment1.output_directory, 'vlab4mic_fig2E_7R5K.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
#######
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(121, projection="3d")
experiment2.particle.gen_axis_plot(with_sources=True, axis_object=ax,target_colour=target_colour, source_plotsize=0.5, emitter_plotsize=10, source_plotalpha=0.1)
ax = fig.add_subplot(122)
ax.imshow(image_outputs2["SMLM"]["ch0"][0], cmap="grey")
filename = os.path.join(experiment2.output_directory, 'vlab4mic_fig2E_8GMO.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
#######
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(121, projection="3d")
experiment3.particle.gen_axis_plot(with_sources=True, axis_object=ax,target_colour=target_colour, source_plotsize=1, emitter_plotsize=1, source_plotalpha=1)
ax = fig.add_subplot(122)
ax.imshow(image_outputs3["SMLM"]["ch0"][0], cmap="grey")
filename = os.path.join(experiment3.output_directory, 'vlab4mic_fig2E_3J3Y.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
#######
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(121, projection="3d")
experiment4.particle.gen_axis_plot(with_sources=True, axis_object=ax,target_colour=target_colour, source_plotsize=2, emitter_plotsize=1, source_plotalpha=1)
ax = fig.add_subplot(122)
ax.imshow(image_outputs4["SMLM"]["ch0"][0], cmap="grey")
filename = os.path.join(experiment4.output_directory, 'vlab4mic_fig2E_1XI5.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()