from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt
from IPython.utils import io
from matplotlib.transforms import Bbox
import os
import numpy as np
from vlab4mic.analysis.metrics import zoom_img
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar 
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
figx = 15
figy = 10
fig = plt.figure(figsize=[figy,figx])
ax_general = plt.axes(projection='3d')
#ax_general = fig.add_subplot(111, projection="3d")
hide_axes=True
experiment1.particle.transform_translate(np.array([0,0,0]))
experiment1.particle.gen_axis_plot(axesoff=hide_axes, with_sources=True, axis_object=ax_general,target_colour=target_colour, source_plotsize=5, emitter_plotsize=5)

experiment4.particle.transform_translate(np.array([0,1400,0]))
experiment4.particle.gen_axis_plot(axesoff=hide_axes,with_sources=True, axis_object=ax_general,target_colour=target_colour, source_plotsize=2, emitter_plotsize=0.5, source_plotalpha=1)

experiment2.particle.transform_translate(np.array([0,2600,0]))

experiment2.particle.gen_axis_plot(axesoff=hide_axes, with_sources=True, axis_object=ax_general,target_colour=target_colour, source_plotsize=0.5, emitter_plotsize=10, source_plotalpha=0.05)

experiment3.particle.transform_translate(np.array([0,3800,0]))
experiment3.particle.gen_axis_plot(axesoff=hide_axes, with_sources=True, axis_object=ax_general,target_colour=target_colour, source_plotsize=1, emitter_plotsize=0.5, source_plotalpha=1)
ax_general.set_zlim3d(bottom=-300, top=1000)
ax_general.set_ylim3d(bottom=-500, top=4500)
ax_general.set_xlim3d(left=0, right=1000)

ax_general.set_box_aspect(
        [
            ub - lb
            for lb, ub in (
                getattr(ax_general, f"get_{a}lim")() for a in "xyz"
            )
        ],
    )
x0,x1 = figx*0.08, figx*0.59
y0,y1 = figy*0.64, figy*0.82

bbox = Bbox([[x0,y0],[x1,y1]])
filename1 = os.path.join(experiment1.output_directory, 'vlab4mic_fig2_E.png')

#bbox = bbox.transformed(ax_general.transData).transformed(fig.dpi_scale_trans.inverted())
fig.savefig(filename1, dpi=300, bbox_inches=bbox)
fig = plt.figure(figsize=[figy,figx])
ax = fig.add_subplot(141)
img_zoom = 0.6
ax.imshow(zoom_img(image_outputs1["SMLM"]["ch0"][0], img_zoom), cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(142)
ax.imshow(zoom_img(image_outputs4["SMLM"]["ch0"][0],img_zoom), cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(143)
ax.imshow(zoom_img(image_outputs2["SMLM"]["ch0"][0],img_zoom), cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(144)
ax.imshow(zoom_img(image_outputs3["SMLM"]["ch0"][0],img_zoom), cmap="grey")
ax.set_axis_off()
length_nm = 250
mods_scalebar = dict()
nm = 1e-9
pixelsize = (experiment3.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment3.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 20 # to get a ratio of 1:10 aspect 
#print(f"{mod} pixelsize: {pixelsize} nm, scalebar length: {length_px} px for {length_nm} nm")
mods_scalebar
scalebar = AnchoredSizeBar(ax.transData,
                           length_px, 
                           #str(length_px)+' nm', 
                           "",
                           "lower right",
                           borderpad=0,
                           color='white',
                           frameon=False,
                           size_vertical=hight_px)
ax.add_artist(scalebar)

plt.tight_layout()

x0,x1 = figx*0, figx*0.67
y0,y1 = figy*0.6, figy*0.89

bbox = Bbox([[x0,y0],[x1,y1]])
filename2 = os.path.join(experiment1.output_directory, 'vlab4mic_fig2_F.png')
fig.savefig(filename2,dpi=300, bbox_inches=bbox)