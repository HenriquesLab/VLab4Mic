from vlab4mic import experiments
import matplotlib.pyplot as plt
from vlab4mic.analysis.metrics import zoom_img
import os
random_seed = 1

modalities = ["SMLM", ]
images, noiselees,  myexp = experiments.image_vsample(
    structure="3J3Y",
    probe_template="HIV_capsid_p24_direct",
    structural_integrity=1,
    structural_integrity_small_cluster=20,
    structural_integrity_large_cluster=100,
    multimodal=modalities,
    clear_experiment=True,
    random_seed=random_seed,
    run_simulation=False
)

myexp.set_structural_integrity(structural_integrity=1)
myexp.build(modules=["particle", "coordinate_field", "imager"])
images, noiseless = myexp.run_simulation()

fig = plt.figure(figsize=[10,10])
ax1 = fig.add_subplot(321)
ax1.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax1.set_axis_off()
ax2 = fig.add_subplot(322, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax2, emitter_plotsize=1, view_init=[20,0,0])

###### 
myexp.set_structural_integrity(structural_integrity=0.6)
myexp.build(modules=["particle", "coordinate_field", "imager"],)
images, noiseless = myexp.run_simulation()

ax3 = fig.add_subplot(323)
ax3.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax3.set_axis_off()
ax4 = fig.add_subplot(324, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax4, emitter_plotsize=1, view_init=[20,0,0])
#### 
myexp.set_structural_integrity(structural_integrity=0.4)
myexp.build(modules=["particle", "coordinate_field", "imager"])
images, noiseless = myexp.run_simulation()
ax5 = fig.add_subplot(325)
ax5.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax5.set_axis_off()
ax6 = fig.add_subplot(326, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax6, emitter_plotsize=1, view_init=[20,0,0])
plt.tight_layout()

filename = myexp.date_as_string + 'vlab4mic_fig3_panelc_integrity.png'
filename2 = os.path.join(myexp.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()


myexp.set_structural_integrity(structural_integrity=1)
myexp.remove_probes()
myexp.clear_labelled_structure()
myexp.add_probe(
    probe_template="HIV_capsid_p24_direct",
    labelling_efficiency=1)
myexp.build(modules=["particle", "coordinate_field", "imager"])
images, noiseless = myexp.run_simulation()

fig = plt.figure(figsize=[10,10])
ax1 = fig.add_subplot(321)
ax1.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax1.set_axis_off()
ax2 = fig.add_subplot(322, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax2, emitter_plotsize=1, view_init=[20,0,0])

###### 
myexp.remove_probes()
myexp.clear_labelled_structure()
myexp.add_probe(
    probe_template="HIV_capsid_p24_direct",
    labelling_efficiency=0.5)
myexp.build(modules=["particle", "coordinate_field", "imager"])
images, noiseless = myexp.run_simulation()

ax3 = fig.add_subplot(323)
ax3.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax3.set_axis_off()
ax4 = fig.add_subplot(324, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax4, emitter_plotsize=1, view_init=[20,0,0])
#### 
myexp.remove_probes()
myexp.clear_labelled_structure()
myexp.add_probe(
    probe_template="HIV_capsid_p24_direct",
    labelling_efficiency=0.2)
myexp.build(modules=["particle", "coordinate_field", "imager"])
images, noiseless = myexp.run_simulation()
ax5 = fig.add_subplot(325)
ax5.imshow(zoom_img(images[modalities[0]]["ch0"][0], 0.6), cmap="grey")
ax5.set_axis_off()
ax6 = fig.add_subplot(326, projection="3d")
myexp.coordinate_field.molecules[0].gen_axis_plot(axis_object=ax6, emitter_plotsize=1, view_init=[20,0,0])
plt.tight_layout()

filename = myexp.date_as_string + 'vlab4mic_fig3_panelc_efficiency.png'
filename2 = os.path.join(myexp.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()