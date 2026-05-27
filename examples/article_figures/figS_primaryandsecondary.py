from vlab4mic.experiments import image_vsample
from vlab4mic.analysis import metrics
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.pyplot as plt
import os
import numpy as np
random_seed = 24

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
    run_simulation=True
)

## Primary and Secondary probes

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

image_outputs2, image_outputs_noiseless2, experiment2 = image_vsample(
    structure = "3J3Y",
    primary_probe = primary,
    secondary_probe = secondary_nanobody,
    multimodal=modalities,
    clear_experiment=True,
    run_simulation=True,
    random_seed=random_seed
)


## create a plot 

emitter_plotsize = 1
zoom = 0.5
target_colour="#01579D"
labelled_figure = plt.figure(figsize=[10,10])
### Primary 
ax1 = labelled_figure.add_subplot(3, 2, 1, projection="3d")
experiment1.particle.gen_axis_plot(
    with_sources=True, 
    axis_object=ax1,
    target_colour=target_colour,
    source_plotsize=2, 
    emitter_plotsize=5, 
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700])
ax2 = labelled_figure.add_subplot(3, 2, 2)
capsid_primary_smlm = metrics.zoom_img(image_outputs1["SMLM"]["ch0"][0], zoom_in=zoom)
ax2.imshow(capsid_primary_smlm, cmap="grey")
ax2.set_axis_off()

length_nm = 100
nm = 1e-09
pixelsize = (experiment1.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment1.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10

scalebar1 = AnchoredSizeBar(ax2.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax2.add_artist(scalebar1)

### Primary and secondary nanobody
ax3 = labelled_figure.add_subplot(3, 2, 3, projection="3d", sharex=ax1, sharey=ax1)
experiment2.particle.gen_axis_plot(
    with_sources=True,
    axis_object=ax3,
    target_colour=target_colour,
    source_plotsize=2,
    emitter_plotsize=5, 
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700])
ax4 = labelled_figure.add_subplot(3, 2, 4)
capsid_primary_secondary_smlm = metrics.zoom_img(image_outputs2["SMLM"]["ch0"][0], zoom_in=zoom)
ax4.imshow(capsid_primary_secondary_smlm, cmap="grey")
ax4.set_axis_off()
scalebar2 = AnchoredSizeBar(ax4.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax4.add_artist(scalebar2)
### Primary and secondary antibody
experiment2.remove_probes()
experiment2.clear_labelled_structure()
experiment2.add_probe(as_primary=True, **primary)
experiment2.add_probe(as_primary=False, **secondary_antibody)
experiment2.build(modules=["particle", "coordinate_field", "imager"])
image_outputs3, noiseless3 = experiment2.run_simulation()
ax5 = labelled_figure.add_subplot(3, 2, 5, projection="3d")
experiment2.particle.gen_axis_plot(
    with_sources=True,
    axis_object=ax5,
    target_colour=target_colour,
    source_plotsize=2,
    emitter_plotsize=5,
    source_plotalpha=1,
    xlim=[-200,800],
    ylim=[0,1200],
    zlim=[0,700])
ax6 = labelled_figure.add_subplot(3, 2, 6)
capsid_primary_secondary_smlm_2 = metrics.zoom_img(image_outputs3["SMLM"]["ch0"][0], zoom_in=zoom)
ax6.imshow(capsid_primary_secondary_smlm_2, cmap="grey")
ax6.set_axis_off()
scalebar3 = AnchoredSizeBar(ax6.transData,
                           length_px, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax6.add_artist(scalebar3)

plt.tight_layout()

# save figure
name = experiment1.date_as_string + 'FigureS_primaryandsecondary.png'
filename = os.path.join(experiment1.output_directory, name)
labelled_figure.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
