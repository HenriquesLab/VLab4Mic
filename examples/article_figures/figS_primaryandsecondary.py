from vlab4mic.experiments import image_vsample
from vlab4mic.analysis import metrics
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.pyplot as plt
import os

modalities = ["STED", "SMLM",]

## Primary probe only
image_outputs1, image_outputs_noiseless1, experiment1 = image_vsample(
    structure = "3J3Y",
    probe_template = "Antibody",
    probe_target_type = "Sequence",
    probe_target_value = "SPRTLNA",
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
secondary = dict(
    probe_template = "Nanobody",
    probe_name="Secondary-HIV",
    probe_target_type = "Primary",
    probe_target_value = "Primary-HIV"
)
image_outputs2, image_outputs_noiseless2, experiment2 = image_vsample(
    structure = "3J3Y",
    primary_probe = primary,
    secondary_probe = secondary,
    multimodal=modalities,
    clear_experiment=True,
    run_simulation=True
)

## create a plot 

emitter_plotsize = 1
zoom = 0.5
target_colour="#01579D"
labelled_figure = plt.figure(figsize=[10,10])
###
ax1 = labelled_figure.add_subplot(2, 2, 1, projection="3d")
experiment1.particle.gen_axis_plot(with_sources=True, axis_object=ax1,target_colour=target_colour, source_plotsize=2, emitter_plotsize=1, source_plotalpha=1)
ax2 = labelled_figure.add_subplot(2, 2, 2)
capsid_primary_smlm = metrics.zoom_img(image_outputs1["SMLM"]["ch0"][0], zoom_in=zoom)
ax2.imshow(capsid_primary_smlm, cmap="grey")
ax2.set_axis_off()
scalebar1 = AnchoredSizeBar(ax2.transData,
                           20, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax2.add_artist(scalebar1)
### 
ax3 = labelled_figure.add_subplot(2, 2, 3, projection="3d", sharex=ax1, sharey=ax1)
experiment2.particle.gen_axis_plot(with_sources=True, axis_object=ax3,target_colour=target_colour, source_plotsize=2, emitter_plotsize=1, source_plotalpha=1)
ax4 = labelled_figure.add_subplot(2, 2, 4)
capsid_primary_secondary_smlm = metrics.zoom_img(image_outputs2["SMLM"]["ch0"][0], zoom_in=zoom)
ax4.imshow(capsid_primary_secondary_smlm, cmap="grey")
ax4.set_axis_off()
scalebar2 = AnchoredSizeBar(ax4.transData,
                           20, '100 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax4.add_artist(scalebar2)
plt.tight_layout()

# save figure
filename = os.path.join(experiment1.output_directory, 'FigureS_primaryandsecondary.png')
labelled_figure.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
