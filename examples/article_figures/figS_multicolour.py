from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.axes_rgb import RGBAxes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar 
import numpy as np
from vlab4mic.analysis.metrics import zoom_img
np.random.seed(44)

# parameters for simulation
structure = "3JA9"
modalities = ["SMLM",]
probe_1 = dict(
    probe_template = "Linker",
    probe_name = "pcnb_A",
    probe_target_type = "Atom_residue",
    probe_target_value = {
            "atoms": ["CA", ],
            "residues": ["SER", ],
            "position": 186,
            "chains": ["A",] # specific for chain A
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
probe_2 = dict(
    probe_template = "Linker",
    probe_name = "pcnb_B",
    probe_target_type = "Atom_residue",
    probe_target_value = {
            "atoms": ["CA", ],
            "residues": ["SER", ],
            "position": 186,
            "chains": ["B",] # specific for chain B
    },
    probe_distance_to_epitope = 0,
    probe_fluorophore = "fluorophore2",
    fluorophore_parameters = {
        "emission":{
            "type": "constant",
            "photon_yield": 100000
        }
    }
)
probe_3 = dict(
    probe_template = "Linker",
    probe_name = "pcnb_C",
    probe_target_type = "Atom_residue",
    probe_target_value = {
            "atoms": ["CA", ],
            "residues": ["SER", ],
            "position": 186,
            "chains": ["C",] # specific for chain C
    },
    probe_distance_to_epitope = 0,
    probe_fluorophore = "fluorophore3",
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
    probe_list=[probe_1, probe_2, probe_3],
    multimodal=modalities,
    SMLM={"exp_time":1, 
          "channels":["ch0","ch1","ch2"], 
          "lateral_resolution_nm":10, # PSF sigma
          "pixelsize_nm":2, # each pixel in the final image represents 2nm
          "psf_voxel_nm":1, # sampling rate for PSF model
          "depth_of_field_nm":200},
    sample_dimensions=[500,500,100],
    run_simulation=False,
    clear_experiment=True,
    random_orientations=True,
    yz_orientations = [150,],
)
# expand sample for better visualisation and run simulation
my_experiment.coordinate_field.expand_isotropically(factor=5.0)
my_experiment.build(modules=["imager",])
images, noisless_images = my_experiment.run_simulation()
# prepare output image for multichannel plot
zoom_factor = 0.8 # focus on the center of the image
modalitiy_name = modalities[0]
R = zoom_img(images[modalitiy_name]["ch0"][0], zoom_in=zoom_factor)
R = R / np.max(R)
G = zoom_img(images[modalitiy_name]["ch1"][0], zoom_in=zoom_factor)
G = G / np.max(G)
B = zoom_img(images[modalitiy_name]["ch2"][0], zoom_in=zoom_factor)
B = B / np.max(B)
# plot 
fig = plt.figure()
ax = RGBAxes(fig, [0.1, 0.1, 0.8, 0.8], pad=0.0)
ax.imshow_rgb(R, G, B)
ax.RGB.set_xticklabels([])
ax.RGB.set_yticklabels([])
scalebar = AnchoredSizeBar(ax.RGB.transData,
                           10, '20 nm', 'lower right', # we know each pixel is 2nm
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2)
ax.RGB.add_artist(scalebar)
ax.RGB.text(x=22, y=3, s="Merge", c="white", fontsize=15)
ax.R.text(x=10, y=10, s="Chain A", c="white", fontsize=15)
ax.G.text(x=10, y=10, s="Chain B", c="white", fontsize=15)
ax.B.text(x=10, y=10, s="Chain C", c="white", fontsize=15)

# save figure
filename = os.path.join(my_experiment.output_directory, 'VLab4Mic_multicolour_example.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()