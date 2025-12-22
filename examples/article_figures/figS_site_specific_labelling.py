from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.axes_rgb import RGBAxes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar 
import numpy as np
from vlab4mic.analysis.metrics import zoom_img
np.random.seed(44)

# parameters for simulation
structure = "4ZTD"
modalities = ["SMLM1", "SMLM2"]
probe_1 = dict(
    probe_template = "Linker",
    probe_name = "pcnb_A",
    probe_target_type = "Atom_residue",
    probe_target_value = {
            "atoms": ["CA", ],
            "residues": ["SER", ],
            "position": 186,
    },
    probe_distance_to_epitope = 0,
    probe_fluorophore = "fluorophore1",
    fluorophore_parameters = {
        "emission":{
            "type": "constant",
            "photon_yield": 100000
        }
    }
)
# Create experiment 
images, noisless_images, my_experiment = image_vsample(
    structure=structure,
    probe_list=[probe_1,],
    multimodal=modalities,
    SMLM1={"exp_time":1, 
          "channels":["ch0",], 
          "lateral_resolution_nm":10, # PSF sigma
          "pixelsize_nm":2, # each pixel in the final image represents 2nm
          "psf_voxel_nm":1, # sampling rate for PSF model
          "depth_of_field_nm":500,
          "modality_template":"SMLM"},
    SMLM2={"exp_time":1, 
          "channels":["ch0",], 
          "lateral_resolution_nm":5, # PSF sigma
          "pixelsize_nm":2, # each pixel in the final image represents 2nm
          "psf_voxel_nm":1, # sampling rate for PSF model
          "depth_of_field_nm":500,
          "modality_template":"SMLM"},
    sample_dimensions=[500,500,100],
    run_simulation=True,
    clear_experiment=True,
)

# prepare output image for multichannel plot
zoom_factor = 0.88 # focus on the center of the image
modalitiy_name1 = modalities[0]
image_zoomed1 = zoom_img(
    images[modalitiy_name1]["ch0"][0], 
    zoom_in=zoom_factor,
)
modalitiy_name2 = modalities[1]
image_zoomed2 = zoom_img(
    images[modalitiy_name2]["ch0"][0], 
    zoom_in=zoom_factor,
)
# create figure
plt.rcParams['figure.figsize'] = [15, 5]

fig = plt.figure(figsize=[10,10])
ax1 = fig.add_subplot(1,3,1, projection="3d")
my_experiment.structure.show_target_labels(
    view_init=[90,0,0],
    show_axis=False,
    with_assembly_atoms=True,
    assembly_fraction=1,
    reference_point=False,
    axis_object=ax1,
    target_size=30,
    axesoff=False
    )
ax1.set_xlabel("Angstroms")
ax1.set_zticks([])
ax1.set_yticklabels([])
ax1.set_xticklabels([0,20,40,60,80,100])
# 
ax2 = fig.add_subplot(1,3,2)
ax2.imshow(image_zoomed1, cmap="grey")
ax2.set_axis_off()
scalebar = AnchoredSizeBar(ax2.transData,
                           10, '20 nm', 'lower right', # pixelsize is 2 nm/pixel
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=1)
ax2.add_artist(scalebar)
#
ax3 = fig.add_subplot(1,3,3)
ax3.imshow(image_zoomed2, cmap="grey")
ax3.set_axis_off()
scalebar = AnchoredSizeBar(ax3.transData,
                           10, '20 nm', 'lower right', # pixelsize is 2 nm/pixel
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=1)
ax3.add_artist(scalebar)

# save figure
filename = os.path.join(my_experiment.output_directory, 'FigureS_site_specific_labelling.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()