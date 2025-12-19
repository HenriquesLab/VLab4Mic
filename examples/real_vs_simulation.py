from vlab4mic import experiments
import requests
import tifffile as tif
import io
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# parameters for simulation
structure = "7R5K"
probe_template = "GFP_w_nanobody"
probe_target_type = "Sequence"
probe_target_value = "ELAVGSL" # Sequence in the C-terminal of Nup96
modalities = ["STED_Thev2016", ]

# fetch experimental image from URL 
STED_example = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/STED/GFP_Gallery-STED_181026_1.tif"
print("Downloading experimental image from Thevathasan et al. 2019")
print("#### If you find this image useful, please consider citing: ####")
print("Jervis Vermal Thevathasan, Maurice Kahnwald, Konstanty Cieśliński, Philipp Hoess, Sudheer Kumar Peneti, Manuel Reitberger, Daniel Heid, Krishna Chaitanya Kasuba, Sarah Janice Hoerner, Yiming Li, Yu-Le Wu, Markus Mund, Ulf Matti, Pedro Matos Pereira, Ricardo Henriques, Bianca Nijmeijer, Moritz Kueblbeck, Vilma Jimenez Sabinina, Jan Ellenberg & Jonas Ries (2019). Nuclear pores as versatile reference standards for quantitative superresolution microscopy . BioStudies, S-BIAD8. Retrieved from https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD8")
print("########")
resp = requests.get(STED_example)
img_tif = tif.imread(io.BytesIO(resp.content))

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

# Prepare experiment object
_1, _2, experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    probe_target_type=probe_target_type,
    probe_target_value=probe_target_value,
    multimodal=modalities,
    STED_Thev2016={"exp_time": 0.0004},
    clear_experiment=True,
    run_simulation=False,
)

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
images, noiselsess = experiment.run_simulation()

# Plot experimental vs simulated image


plt.rcParams['figure.figsize'] = [20, 10]

fig, axs = plt.subplots(1, 2)
axs[0].imshow(experimental_img_processed, cmap="grey")
axs[0].set_title("Modified from Thevathasan et al. 2019. Nature Methods.")
axs[1].imshow(images["STED_Thev2016"]["ch0"][0], cmap="grey")
axs[1].set_title("VLab4Mic simulated image")

fontprops = fm.FontProperties(size=20)
scalebar = AnchoredSizeBar(axs[1].transData,
                           67, '1 µm', 'lower right', 
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
# save figure
filename = os.path.join(experiment.output_directory, 'vlab4mic_real_vs_sim.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()