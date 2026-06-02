from vlab4mic import experiments
from vlab4mic.experiments import generate_virtual_sample
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

smlm_experimental_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/dSTORM2D/SNAP_Gallery-STORM-2D_181214_1_sml.csv"

localisations_all = pd.read_csv(smlm_experimental_path)
smlm_experimental_locs = localisations_all[localisations_all.numberInGroup == 1]

# SMLM
loc_prec = np.max(smlm_experimental_locs.locprecnm)
xrange = smlm_experimental_locs.xnm.max() - smlm_experimental_locs.xnm.min()
yrange = smlm_experimental_locs.ynm.max() - smlm_experimental_locs.ynm.min()
####################
def render_from_localisations(localisations, crop_size_nm = 3000, centerx = None, centery = None,precision_nm=None, random_seed=None, pixelsize=10, voxelsize=5):
    if precision_nm is None:
        precision_nm = np.mean(localisations.locprecnm)
    minx = np.min(localisations.xnm)
    maxx = np.max(localisations.xnm)
    miny = np.min(localisations.ynm)
    maxy = np.max(localisations.ynm)
    xrange_nm = maxx - minx
    yrange_nm = maxy - miny
    newlocs_x = localisations.xnm - minx
    newlocs_y = localisations.ynm - miny
    if centerx is not None:
        centerx_nm = centerx
    else:
        centerx_nm = xrange_nm / 2
        print(centerx_nm, xrange_nm )
    if centery is not None:
        centery_nm = centery
    else:
        centery_nm = yrange_nm / 2
        print(centery_nm, yrange_nm)
    crop_minx = centerx_nm - crop_size_nm / 2
    crop_maxx = centerx_nm + crop_size_nm / 2
    crop_miny = centery_nm - crop_size_nm / 2
    crop_maxy = centery_nm + crop_size_nm / 2
    x_cropped = []
    y_cropped = []
    for x,y in zip(newlocs_x, newlocs_y):
        if (x >= crop_minx) and (x <= crop_maxx) and (y >= crop_miny) and (y <= crop_maxy):
            x_cropped.append(x)
            y_cropped.append(y)
    normalised_x = (np.array(x_cropped) - np.array(x_cropped).min()) / np.array(x_cropped).max()
    normalised_x = normalised_x / normalised_x.max()
    normalised_y = (np.array(y_cropped) - np.array(y_cropped).min()) / np.array(y_cropped).max()
    normalised_y = normalised_y / normalised_y.max()
    sample, experiment = generate_virtual_sample(
        structure=None,
        particle_positions=(normalised_x, normalised_y),
        random_seed=random_seed
    )
    experiment.clear_structure()
    experiment.remove_probes()
    experiment.clear_labelled_structure()
    experiment.clear_virtual_sample()
    list_of_positions = [ [x,y,z] for x,y,z in zip(normalised_x, normalised_y, np.zeros_like(normalised_x))]
    experiment.set_virtualsample_params(
        sample_dimensions=(crop_size_nm, crop_size_nm, 100),
        particle_positions=list_of_positions
    )
    experiment.add_modality("SMLM", modality_template="SMLM")
    experiment.build(modules=["coordinate_field", "imager",])
    experiment.update_modality(modality_name="SMLM",lateral_resolution_nm = precision_nm, simulate_localistations=False, pixelsize_nm=pixelsize, psf_voxel_nm=voxelsize)
    experiment.set_modality_acq(modality_name="SMLM", exp_time=1)
    #images = {}
    images, noiseless = experiment.run_simulation(modality="SMLM")
    return list_of_positions, images, experiment

FOV_nm = 4000
smlm_pixelsize = 5
list_of_positions, rendered_smlm, smlm_experiment = render_from_localisations(
    smlm_experimental_locs,
    crop_size_nm = FOV_nm,
    centery = yrange*0.7,
    centerx = xrange*0.7,
    precision_nm=loc_prec,
    pixelsize=smlm_pixelsize,
    voxelsize=smlm_pixelsize,
    random_seed=random_seed
    )

# parameters for simulation
structure = "7R5K"
probe_template = "GFP_w_nanobody"
probe_target_type = "Sequence"
probe_target_value = "ELAVGSL" # Sequence in the C-terminal of Nup96
modalities = ["SMLM", ]

# Prepare experiment object
_1, _2, experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    probe_target_type=probe_target_type,
    probe_target_value=probe_target_value,
    probe_DoL=2,
    labelling_efficiency=0.8,
    multimodal=modalities,
    SMLM={"exp_time": 0.0002},
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)

experimental_img = np.array(rendered_smlm["SMLM"]["ch0"][0], dtype=int)
background = 0
sigma = 5
min_distance = 20*smlm_pixelsize
threshold = 0.01

# Use experimental image for positioning
experiment.use_image_for_positioning(
    img=experimental_img,
    mode="localmaxima",
    sigma=sigma,
    background=background,
    threshold=threshold,
    pixelsize=smlm_pixelsize,
    min_distance=min_distance)

# Run simulation
experiment.update_modality(modality_name="SMLM", simulate_localisations=False, pixelsize_nm=smlm_pixelsize, lateral_resolution_nm=5, psf_voxel_nm=smlm_pixelsize)
experiment.set_modality_acq(modality_name="SMLM", exp_time=0.0002)
images, noiselsess = experiment.run_simulation()


# Analyse images: fit circles
# expected parameters of objects
expected_min_radius = 40  # in nanometers
expected_MAX_radius = 70  # in nanometers
min_radius_px = expected_min_radius / smlm_pixelsize
MAX_radius_px = expected_MAX_radius / smlm_pixelsize
minRadius_round = np.ceil(min_radius_px).astype('int64')
maxRadius_round = np.ceil(MAX_radius_px).astype('int64')
minDist = minRadius_round
HCparams = dict(dp=1, minDist=maxRadius_round, 
        param1=1, param2=7, minRadius=minRadius_round, maxRadius=maxRadius_round)
if images["SMLM"]["ch0"][0].min() < 0:
    images["SMLM"]["ch0"][0] += -images["SMLM"]["ch0"][0].min()

#### Simulated Data
circles_sim, img_blurred_sim, c_params_sim = metrics.get_circles(images["SMLM"]["ch0"][0].astype(np.uint8), **HCparams)
radii_simulated= []
for (x, y, r) in circles_sim[0]:
    radii_simulated.append((r*smlm_pixelsize))

#### Experimental Data
circles_exp, img_blurred_exp, c_params_exp = metrics.get_circles(experimental_img.astype(np.uint8), **HCparams)
radii_experimental = []
for (x, y, r) in circles_exp[0]:
    radii_experimental.append((r*smlm_pixelsize))


# Plot experimental vs simulated image

plt.rcParams['figure.figsize'] = [15, 5]

fig, axs = plt.subplots(1, 3)
axs[0].imshow(experimental_img, cmap="grey")
axs[0].set_title("Experimental example SMLM")
axs[1].imshow(images["SMLM"]["ch0"][0], cmap="grey")
axs[1].set_title("VLab4Mic simulated image SMLM")

df = pd.DataFrame({
    'value': np.concatenate([radii_experimental, radii_simulated]),
    'Condition': ['Experimental'] * len(radii_experimental) + ['Simulated'] * len(radii_simulated)
})
sns.histplot(data=df, x="value", hue="Condition", binrange=[45,70], bins=15, kde=True, ax=axs[2])
plt.xlabel("Radius of circle fit (nm)")

length_nm = 1000
nm = 1e-09
pixelsize = (experiment.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * experiment.imaging_modalities["SMLM"]["detector"]["pixelsize"]
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
name = experiment.date_as_string + 'NPC_quantification_SMLM.png'
filename = os.path.join(experiment.output_directory, name)
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()