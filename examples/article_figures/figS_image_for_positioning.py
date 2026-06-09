import matplotlib.pyplot as plt
from vlab4mic import experiments
import numpy as np
from IPython.utils import io
import requests
import tifffile as tif

import io
import os
random_seed=24

img, noiseless, my_experiment = experiments.image_vsample(
    structure="7R5K",
    probe_template = "GFP_w_nanobody",
    probe_target_type="Sequence",
    probe_target_value="ELAVGSL", # Nup96 C-terminal
    particle_positions=[[0.5,0.5,0.5],],
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)

fig = plt.figure(figsize=[20,10])
ax = fig.add_subplot(1, 3, 1, projection="3d")
my_experiment.particle.gen_axis_plot(axis_object=ax)
ax = fig.add_subplot(1, 3, 2, projection="3d")
my_experiment.coordinate_field.show_field(axis_object=ax, axesoff=True, emitters_plotsize=5, initial_pos=False, view_init=[20,0,0])
my_experiment.clear_virtual_sample()
my_experiment.set_virtualsample_params(
    particle_positions=[
        [0.2,0.9,0.0],
        [0.3,0.3,0.6],
        [0.5,0.6,0.1],
        [0.8,0.2,0.8],
        [0.8,0.7,0.8]
        ]
)
my_experiment.build(modules=["coordinate_field",])
ax = fig.add_subplot(1, 3, 3, projection="3d")
my_experiment.coordinate_field.show_field(axis_object=ax, axesoff=True, emitters_plotsize=5, initial_pos=False, view_init=[20,0,0])

filename_ = my_experiment.date_as_string + 'figS_vsample_construction.pdf'
filename = os.path.join(my_experiment.output_directory, filename_)
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()


rng_ = np.random.default_rng(seed=random_seed)
img_mask = rng_.random([100,100])
p = 0.99
img_mask[img_mask >= p] = 1
img_mask[img_mask < p] = 0
npixels = np.sum(img_mask)
my_experiment.use_image_for_positioning(img_mask, mode="mask", pixelsize=15, npositions=npixels, min_distance=200)

fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(1, 2, 1)
ax.imshow(img_mask, cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(1, 2, 2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], initial_pos=False, axis_object=ax, emitters_plotsize=5, axesoff=True)

filename_ = my_experiment.date_as_string + 'figS_vsample_position_with_mask.pdf'
filename = os.path.join(my_experiment.output_directory, filename_)
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()

def read_image_from_url(url):
    try:
        resp = requests.get(url)
    except:
        resp = requests.get(url)
    # Check that request succeeded
    return tif.imread(io.BytesIO(resp.content))

sted_experimental_img_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/STED/GFP_Gallery-STED_181026_1.tif"
sted_experimental_img =  read_image_from_url(sted_experimental_img_path)

centery = int(800) #vertical
centerx = int(700)
width = 50
experimental_image_patch = sted_experimental_img[centery-width:centery+width:,centerx-width:centerx+width].copy()

background = 32768
pixelsize = 15
# pre processing
sigma = 2
min_distance = 6*pixelsize
threshold = 4
my_experiment.use_image_for_positioning(
    img=experimental_image_patch,
    mode="localmaxima",
    background=background,
    pixelsize=pixelsize,
    sigma=sigma, 
    min_distance=min_distance,
    threshold=threshold)

fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(1, 2, 1)
ax.imshow(experimental_image_patch, cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(1, 2, 2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], initial_pos=False, axis_object=ax, emitters_plotsize=5, axesoff=True)


filename_ = my_experiment.date_as_string + 'figS_vsample_position_from_experimentimage.pdf'
filename = os.path.join(my_experiment.output_directory, filename_)
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()