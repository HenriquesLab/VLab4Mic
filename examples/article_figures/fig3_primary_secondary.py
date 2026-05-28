from vlab4mic.experiments import image_vsample
from vlab4mic.analysis import metrics
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.pyplot as plt
import os
import numpy as np
random_seed = 25

modalities = ["STED", "SMLM",]

## Primary probe only
image_outputs1, image_outputs_noiseless1, experiment1 = image_vsample(
    structure = "3J3Y",
    probe_template = "Antibody",
    probe_target_type = "Sequence",
    probe_target_value = "SPRTLNA",
    probe_DoL=4,
    multimodal=modalities,
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)
primary = dict(
    probe_template = "Antibody",
    probe_name="Primary-HIV",
    probe_target_type = "Sequence",
    probe_target_value = "SPRTLNA",
)
secondary_nanobody = dict(
    probe_template = "Nanobody",
    probe_name="Secondary-nanobody",
    probe_target_type = "Primary",
    probe_target_value = "Primary-HIV",
    probe_DoL=2,
)
secondary_antibody = dict(
    probe_template = "Antibody",
    probe_name="Secondary-antibody",
    probe_target_type = "Primary",
    probe_target_value = "Primary-HIV",
    probe_DoL=4,
)
emitter_plotsize = 1
zoom = 0.5
target_colour="#01579D"
labelled_figure = plt.figure(figsize=[10,15])
### Primary 
### Primary and secondary nanobody
experiment1.remove_probes()
experiment1.clear_labelled_structure()
experiment1.add_probe(as_primary=False, **primary)
experiment1.build(modules=["particle", "coordinate_field", "imager"])
image_outputs1, noiseless1 = experiment1.run_simulation()
ax1 = labelled_figure.add_subplot(3, 3, 1, projection="3d")
experiment1.particle.gen_axis_plot(
    with_sources=True, 
    axis_object=ax1,
    target_colour=target_colour,
    source_plotsize=2, 
    emitter_plotsize=10, 
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700],
    view_init=[90,0,0])
ax2 = labelled_figure.add_subplot(3, 3, 2)
capsid_primary_sted1 = metrics.zoom_img(image_outputs1["STED"]["ch0"][0], zoom_in=zoom)
ax2.imshow(capsid_primary_sted1, cmap="grey")
ax2.set_axis_off()

ax3 = labelled_figure.add_subplot(3, 3, 3)
capsid_primary_smlm1 = metrics.zoom_img(image_outputs1["SMLM"]["ch0"][0], zoom_in=zoom)
ax3.imshow(capsid_primary_smlm1, cmap="grey")
ax3.set_axis_off()

length_nm = 100
nm = 1e-09
pixelsize = (experiment1.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment1.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10

scalebar1 = AnchoredSizeBar(ax3.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax3.add_artist(scalebar1)

##############
### Primary and secondary nanobody
experiment1.remove_probes()
experiment1.clear_labelled_structure()
experiment1.add_probe(as_primary=True, **primary)
experiment1.add_probe(as_primary=False, **secondary_nanobody)
experiment1.build(modules=["particle", "coordinate_field", "imager"])
ax4 = labelled_figure.add_subplot(3, 3, 4, projection="3d")
experiment1.particle.gen_axis_plot(
    with_sources=True, 
    axis_object=ax4,
    target_colour=target_colour,
    source_plotsize=2, 
    emitter_plotsize=10, 
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700],
    view_init=[90,0,0])
image_outputs2, noiseless2 = experiment1.run_simulation()
ax5 = labelled_figure.add_subplot(3, 3, 5)
capsid_primary_sted2 = metrics.zoom_img(image_outputs2["STED"]["ch0"][0], zoom_in=zoom)
ax5.imshow(capsid_primary_sted2, cmap="grey")
ax5.set_axis_off()

ax6 = labelled_figure.add_subplot(3, 3, 6)
capsid_primary_smlm2 = metrics.zoom_img(image_outputs2["SMLM"]["ch0"][0], zoom_in=zoom)
ax6.imshow(capsid_primary_smlm2, cmap="grey")
ax6.set_axis_off()

length_nm = 100
nm = 1e-09
pixelsize = (experiment1.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment1.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10

scalebar2 = AnchoredSizeBar(ax6.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax6.add_artist(scalebar2)
##############
### Primary and secondary antibody
experiment1.remove_probes()
experiment1.clear_labelled_structure()
experiment1.add_probe(as_primary=True, **primary)
experiment1.add_probe(as_primary=False, **secondary_antibody)
experiment1.build(modules=["particle", "coordinate_field", "imager"])
ax7 = labelled_figure.add_subplot(3, 3, 7, projection="3d")
experiment1.particle.gen_axis_plot(
    with_sources=True, 
    axis_object=ax7,
    target_colour=target_colour,
    source_plotsize=2, 
    emitter_plotsize=10, 
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700],
    view_init=[90,0,0])
image_outputs3, noiseless3 = experiment1.run_simulation()
ax8 = labelled_figure.add_subplot(3, 3, 8)
capsid_primary_sted3 = metrics.zoom_img(image_outputs3["STED"]["ch0"][0], zoom_in=zoom)
ax8.imshow(capsid_primary_sted3, cmap="grey")
ax8.set_axis_off()

ax9 = labelled_figure.add_subplot(3, 3, 9)
capsid_primary_smlm3 = metrics.zoom_img(image_outputs3["SMLM"]["ch0"][0], zoom_in=zoom)
ax9.imshow(capsid_primary_smlm3, cmap="grey")
ax9.set_axis_off()

length_nm = 100
nm = 1e-09
pixelsize = (experiment1.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment1.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10

scalebar3 = AnchoredSizeBar(ax9.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax9.add_artist(scalebar3)
plt.tight_layout()

# save figure
name = experiment1.date_as_string + 'Fig3_primaryandsecondary_panel.png'
filename = os.path.join(experiment1.output_directory, name)
labelled_figure.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()