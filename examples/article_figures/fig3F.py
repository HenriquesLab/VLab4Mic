import numpy as np
import os
import requests
import matplotlib.pyplot as plt
from vlab4mic import experiments
from vlab4mic.experiments import generate_virtual_sample
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import tifffile as tif
import io
import pandas as pd
np.random.seed(44)
sim_vs_exp = dict() 

def mimic_experimental_image(experiment, experimental_image, halfpatch, centerx, centery, as_int=True, **kwargs):
    if as_int:
        experimental_image_int = np.array(experimental_image, dtype=int)
    else:
        experimental_image_int = experimental_image
    if halfpatch is not None:
        experimental_image_patch = experimental_image_int[centery-halfpatch:centery+halfpatch:,centerx-halfpatch:centerx+halfpatch]
    else:
        experimental_image_patch = experimental_image_int
    experiment.use_image_for_positioning(
        img=experimental_image_patch,
        background=experimental_image_patch.min(),
        **kwargs
    )
    return experimental_image_patch, experiment.virtualsample_params["relative_positions"]

def read_image_from_url(url):
    try:
        resp = requests.get(url)
    except:
        resp = requests.get(url)
    # Check that request succeeded
    return tif.imread(io.BytesIO(resp.content))

def render_from_localisations(localisations, crop_size_nm = 3000, centerx = None, centery = None,precision_nm=None):
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
        particle_positions=(normalised_x, normalised_y)
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
    experiment.build(modules=["coordinate_field", "imager",])
    experiment.update_modality(modality_name="SMLM",lateral_resolution_nm = precision_nm)
    experiment.set_modality_acq(modality_name="SMLM", exp_time=0.01)
    #images = {}
    images, noiseless = experiment.run_simulation(modality="SMLM")
    return list_of_positions, images, experiment


# widefield
wf_experimental_img_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/Widefield/GFP_Gallery-Widefield_190124_1.tif"
wf_experimental_img =  read_image_from_url(wf_experimental_img_path)
# confocal
confocal_experimental_img_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/Confocal/GFP_Gallery-Confocal_180917_1.tif"
confocal_experimental_img =  read_image_from_url(confocal_experimental_img_path)
# AiryScan
airy_experimental_img_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/Airy/GFP_Gallery-Airy_181122_1.tif"
airy_experimental_img =  read_image_from_url(airy_experimental_img_path)
# STED
sted_experimental_img_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/STED/GFP_Gallery-STED_181026_1.tif"
sted_experimental_img =  read_image_from_url(sted_experimental_img_path)
# SMLM
smlm_experimental_path = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/dSTORM2D/SNAP_Gallery-STORM-2D_181214_1_sml.csv"
localisations_all = pd.read_csv(smlm_experimental_path)
smlm_experimental_locs = localisations_all[localisations_all.numberInGroup == 1]

# parameters
structure = "7R5K"
probe_template = "NPC_Nup96_Cterminal_direct"
modalities = ["Widefield", "Confocal", "AiryScan","STED", "SMLM"]
# Run simulation
# Prepare experiment object
_1, _2, my_experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    labelling_efficiency=1,
    multimodal=modalities,
    random_rotation=True,
    clear_experiment=True,
    run_simulation=False,
)
FOV_nm = 1500
widefield_pixelsize = 100
confocal_pixelsize = 70
airyscan_pixelsize = 40
sted_pixelsize = 15
smlm_pixelsize = 5


sample_from_image_params = dict(
    Widefield = dict(
        centery = int(167), #vertical
        centerx = int(389),
        pixelsize = widefield_pixelsize, # in nm
        sigma = 0.5, # in pixels
        min_distance = 1*widefield_pixelsize,
        threshold = 0,
        halfpatch = FOV_nm//(2*widefield_pixelsize),
    ),
    Confocal = dict(
        centery = int(294), #vertical
        centerx = int(179), 
        pixelsize = confocal_pixelsize, # in nm
        sigma = 0.5, # in pixels
        min_distance = 1*confocal_pixelsize,
        threshold = 400,
        halfpatch = FOV_nm//(2*confocal_pixelsize),
    ),
    AiryScan = dict(
        centery = int(70), #vertical
        centerx = int(310),
        pixelsize = airyscan_pixelsize, # in nm
        sigma = 1, # in pixels
        min_distance = 1*airyscan_pixelsize,
        threshold = 100,
        halfpatch = FOV_nm//(2*airyscan_pixelsize),
    ),
    STED = dict(
        centery = int(630), #vertical
        centerx = int(1000), 
        pixelsize = sted_pixelsize, # in nm
        sigma = 3, # in pixels
        min_distance = 6*sted_pixelsize,
        threshold = 3,
        halfpatch = FOV_nm//(2*sted_pixelsize),
    ),
    SMLM = dict(
        centery = None,
        centerx = None,
        pixelsize = smlm_pixelsize, # in nm
        sigma = 6, # in pixels
        min_distance = 24*smlm_pixelsize,
        threshold = 10,
        halfpatch = None,
    ),    
)
# Widefield
wf_experimental_img_patch, relative_positions = mimic_experimental_image(
    experiment=my_experiment,
    experimental_image=wf_experimental_img,
    mode="localmaxima",
    exclude_border=False,
    **sample_from_image_params["Widefield"]
)
my_experiment.imager.modalities["Widefield"]["detector"]["noise_model"]["binomial"]["p"] = 0.8
my_experiment.imager.modalities["Widefield"]["detector"]["noise_model"]["baselevel"]["bl"] = 4000
my_experiment.imager.modalities["Widefield"]["detector"]["noise_model"]["gaussian"]["sigma"] = 4000**(1/2)
# Run simulation
my_experiment.update_modality(
    modality_name="Widefield",
    lateral_resolution_nm = np.mean(smlm_experimental_locs.PSFxnm))
my_experiment.set_modality_acq(modality_name="Widefield", exp_time=0.01)
wf_images, noiselsess = my_experiment.run_simulation(modality="Widefield")
sim_vs_exp["Widefield"] = (wf_images["Widefield"]["ch0"][0], wf_experimental_img_patch)
# Confocal
confocal_experimental_img_patch, relative_positions = mimic_experimental_image(
    experiment=my_experiment,
    experimental_image=confocal_experimental_img,
    mode="localmaxima",
    exclude_border=False,
    **sample_from_image_params["Confocal"]
)
my_experiment.update_modality(modality_name="Confocal", lateral_resolution_nm = 100)
my_experiment.set_modality_acq(modality_name="Confocal", exp_time=0.001)
confocal_images, noiselsess = my_experiment.run_simulation(modality="Confocal")
sim_vs_exp["Confocal"] = (confocal_images["Confocal"]["ch0"][0], confocal_experimental_img_patch)
# AiryScan
airy_experimental_img_patch, relative_positions = mimic_experimental_image(
    experiment=my_experiment,
    experimental_image=airy_experimental_img,
    mode="localmaxima",
    exclude_border=True,
    **sample_from_image_params["AiryScan"]
)
my_experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["binomial"]["p"] = 0.8
my_experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["baselevel"]["bl"] = 600
my_experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["gaussian"]["sigma"] = 25
# Run simulation
my_experiment.update_modality(modality_name="AiryScan", lateral_resolution_nm = 60)
my_experiment.set_modality_acq(modality_name="AiryScan", exp_time=0.01)
airy_images, noiselsess = my_experiment.run_simulation(modality="AiryScan")
sim_vs_exp["AiryScan"] = (airy_images["AiryScan"]["ch0"][0], airy_experimental_img_patch)
# STED
STED_experimental_img_patch, relative_positions = mimic_experimental_image(
    experiment=my_experiment,
    experimental_image=sted_experimental_img,
    mode="localmaxima",
    exclude_border=True,
    **sample_from_image_params["STED"]
)
my_experiment.imager.modalities["STED"]["detector"]["noise_model"]["binomial"]["p"] = 0.4
my_experiment.imager.modalities["STED"]["detector"]["noise_model"]["baselevel"]["bl"] = 2
my_experiment.imager.modalities["STED"]["detector"]["noise_model"]["gaussian"]["sigma"] = 0
# Run simulation
my_experiment.update_modality(modality_name="STED", lateral_resolution_nm = 15)
my_experiment.set_modality_acq(modality_name="STED", exp_time=0.001)
sted_images, noiselsess = my_experiment.run_simulation(modality="STED")
sim_vs_exp["STED"] = (sted_images["STED"]["ch0"][0], STED_experimental_img_patch)
# SMLM
loc_prec = np.max(smlm_experimental_locs.locprecnm)
xrange = smlm_experimental_locs.xnm.max() - smlm_experimental_locs.xnm.min()
yrange = smlm_experimental_locs.ynm.max() - smlm_experimental_locs.ynm.min()
list_of_positions, rendered_smlm, smlm_experiment = render_from_localisations(
    smlm_experimental_locs,
    crop_size_nm = FOV_nm,
    centery = yrange*0.6,
    centerx = xrange*0.5,
    precision_nm=loc_prec
    )
SMLM_experimental_img_patch, relative_positions = mimic_experimental_image(
    experiment=my_experiment,
    experimental_image=rendered_smlm["SMLM"]["ch0"][0],
    as_int=False,
    mode="localmaxima",
    exclude_border=False,
    **sample_from_image_params["SMLM"]
)
my_experiment.imager.modalities["SMLM"]["detector"]["noise_model"]["binomial"]["p"] = 0.4
my_experiment.imager.modalities["SMLM"]["detector"]["noise_model"]["baselevel"]["bl"] = 1
my_experiment.imager.modalities["SMLM"]["detector"]["noise_model"]["gaussian"]["sigma"] = 0
# Run simulation
my_experiment.update_modality(modality_name="SMLM", lateral_resolution_nm = loc_prec)
my_experiment.set_modality_acq(modality_name="SMLM", exp_time=0.01)
smlm_images, noiselsess = my_experiment.run_simulation(modality="SMLM")
sim_vs_exp["SMLM"] = (smlm_images["SMLM"]["ch0"][0], SMLM_experimental_img_patch)
### 
length_nm = 500
nm = 1e-09
mods_scalebar = dict()
for mod in modalities:
    pixelsize = (my_experiment.imaging_modalities[mod]["detector"]["scale"] / nm) * my_experiment.imaging_modalities[mod]["detector"]["pixelsize"]
    pixelsize = np.ceil(pixelsize)
    length_px = length_nm / pixelsize
    hight_px = length_px / 10
    mods_scalebar[mod] = {
        "length_px": length_px,
        "hight": hight_px,
    }



plt.rcParams['figure.figsize'] = [24, 10]
nmods = len(modalities)

fig, axs = plt.subplots(2, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[0,i].imshow(sim_vs_exp[mod][1], cmap="grey")
    axs[0,i].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    axs[0,i].set_title(f"{mod}", fontsize=20)
    axs[1,i].imshow(sim_vs_exp[mod][0], cmap="magma")
    axs[1,i].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    scalebar = AnchoredSizeBar(axs[1,i].transData,
                           mods_scalebar[mod]["length_px"],
                           #str(length_nm)+' nm', 
                           '',
                           'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=mods_scalebar[mod]["hight"])
    axs[1,i].add_artist(scalebar)
axs[0,0].set_ylabel("Thevathasan, et al. 2019", fontsize=24)
axs[1,0].set_ylabel("VLab4Mic simulation", fontsize=24)
plt.tight_layout(h_pad=0, w_pad=3)

filename = os.path.join(my_experiment.output_directory, 'vlab4mic_fig3F_url.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()