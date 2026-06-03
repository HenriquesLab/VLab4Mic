from vlab4mic import experiments
from vlab4mic.analysis import metrics
import pandas as pd
import seaborn as sns
import requests
import tifffile as tif
import io
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
random_seed= 24

# parameters for simulation
structure = "7R5K"
probe_template = "GFP_w_nanobody"
probe_target_type = "Sequence"
probe_target_value = "ELAVGSL" # Sequence in the C-terminal of Nup96
modalities = ["STED", ]

# fetch experimental image from URL 
STED_example = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/STED/GFP_Gallery-STED_181026_1.tif"
print("Downloading experimental image from Thevathasan et al. 2019")
print("#### If you find this image useful, please consider citing: ####")
print("Jervis Vermal Thevathasan, Maurice Kahnwald, Konstanty Cieśliński, Philipp Hoess, Sudheer Kumar Peneti, Manuel Reitberger, Daniel Heid, Krishna Chaitanya Kasuba, Sarah Janice Hoerner, Yiming Li, Yu-Le Wu, Markus Mund, Ulf Matti, Pedro Matos Pereira, Ricardo Henriques, Bianca Nijmeijer, Moritz Kueblbeck, Vilma Jimenez Sabinina, Jan Ellenberg & Jonas Ries (2019). Nuclear pores as versatile reference standards for quantitative superresolution microscopy . BioStudies, S-BIAD8. Retrieved from https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD8")
print("########")
resp = requests.get(STED_example)
img_tif = tif.imread(io.BytesIO(resp.content))


# Prepare experiment object
_1, _2, experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    probe_target_type=probe_target_type,
    probe_target_value=probe_target_value,
    multimodal=modalities,
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)

# Get a crop of the experimental image to reduce the field of view
experimental_img = np.array(img_tif, dtype=int)
halfpatch = 200
center = int(experimental_img.shape[0]/2)
experimental_img_processed = experimental_img[center-halfpatch:center+halfpatch:,center-halfpatch:center+halfpatch]
# parameters from experimental image metadata
background = 32768 # background level
pixelsize = 15 # in nm
# parameters for spot detection
sigma = 2 # in pixels
min_distance = 6*pixelsize # in physical size, in nm
threshold = 4

# Use experimental image for positioning
experiment.use_image_for_positioning(
    img=experimental_img_processed,
    mode="localmaxima",
    sigma=sigma,
    background=background,
    threshold=threshold,
    pixelsize=pixelsize,
    min_distance=min_distance)

# Run simulation
## Modality parameter of pixelsize_nm is set to match Thevathasan et al. 2019.
## All other parameters are estimations to emulate Thevathasan et al. 2019 STED modality
experiment.update_modality(
    modality_name="STED",
    lateral_resolution_nm=15,
    axial_resolution_nm=15,
    psf_voxel_nm=5,
    depth_of_field_nm=100,
    pixelsize_nm=15)
experiment.set_modality_acq(
    modality_name="STED",
    exp_time=0.0002)
images, noiselsess = experiment.run_simulation()

# Analyse images: fit circles
# expected parameters of objects
expected_min_radius = 40  # in nanometers
expected_MAX_radius = 70  # in nanometers
min_radius_px = expected_min_radius / pixelsize
MAX_radius_px = expected_MAX_radius / pixelsize
minRadius_round = np.ceil(min_radius_px).astype('int64')
maxRadius_round = np.ceil(MAX_radius_px).astype('int64')
minDist = minRadius_round

if images["STED"]["ch0"][0].min() < 0:
    images["STED"]["ch0"][0] += -images["STED"]["ch0"][0].min()

#### Simulated Data
HCparams = dict(dp=1, minDist=maxRadius_round, 
        param1=10, param2=7, minRadius=minRadius_round, maxRadius=maxRadius_round)
circles_sim, img_blurred_sim, c_params_sim = metrics.get_circles(images["STED"]["ch0"][0].astype(np.uint8), **HCparams)
radii_simulated= []
for (x, y, r) in circles_sim[0]:
    radii_simulated.append((r*pixelsize))

#### Experimental Data
circles_exp, img_blurred_exp, c_params_exp = metrics.get_circles(experimental_img_processed.astype(np.uint8), **HCparams)
radii_experimental = []
for (x, y, r) in circles_exp[0]:
    radii_experimental.append((r*pixelsize))

# Plot experimental vs simulated image

plt.rcParams['figure.figsize'] = [15, 5]

fig, axs = plt.subplots(1, 3)
axs[0].imshow(experimental_img_processed, cmap="grey")
axs[0].set_title("Modified from Thevathasan et al. 2019. Nature Methods.")
axs[1].imshow(images["STED"]["ch0"][0], cmap="grey")
axs[1].set_title("VLab4Mic simulated image")

df = pd.DataFrame({
    'value': np.concatenate([radii_experimental, radii_simulated]),
    'Condition': ['Experimental'] * len(radii_experimental) + ['Simulated'] * len(radii_simulated)
})
sns.histplot(data=df, x="value", hue="Condition", binrange=[45,70], bins=15, kde=True, ax=axs[2])
plt.xlabel("Radius of circle fit (nm)")

length_nm = 1000
nm = 1e-09
pixelsize = (experiment.imaging_modalities["STED"]["detector"]["scale"] / nm) * experiment.imaging_modalities["STED"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10


fontprops = fm.FontProperties(size=20)
scalebar = AnchoredSizeBar(axs[1].transData,
                           length_px, '1 µm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=5,
                           fontproperties=fontprops)
axs[1].add_artist(scalebar)
# remove axis ticks
axs[0].set_xticks([])
axs[0].set_yticks([])
axs[1].set_xticks([])
axs[1].set_yticks([])
plt.tight_layout()

# save figure
name = experiment.date_as_string + 'NPC_quantification_STED.png'
filename = os.path.join(experiment.output_directory, name)
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()