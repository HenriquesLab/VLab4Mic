from vlab4mic import experiments
from vlab4mic.analysis import metrics, particle_measures
import pandas as pd
import seaborn as sns
import requests
import tifffile as tif
import io
import os
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
random_seed= 24

# parameters for simulation
structure = "7R5K"
probe_template = "GFP_w_nanobody"
probe_target_type = "Sequence"
probe_target_value = "ELAVGSL" # Sequence in the C-terminal of Nup96
modalities = ["AiryScan", ]

# fetch experimental image from URL
AiryScan_example = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/Airy/GFP_Gallery-Airy_181122_1.tif"
print("Downloading experimental image from Thevathasan et al. 2019")
print("#### If you find this image useful, please consider citing: ####")
print("Jervis Vermal Thevathasan, Maurice Kahnwald, Konstanty Cieśliński, Philipp Hoess, Sudheer Kumar Peneti, Manuel Reitberger, Daniel Heid, Krishna Chaitanya Kasuba, Sarah Janice Hoerner, Yiming Li, Yu-Le Wu, Markus Mund, Ulf Matti, Pedro Matos Pereira, Ricardo Henriques, Bianca Nijmeijer, Moritz Kueblbeck, Vilma Jimenez Sabinina, Jan Ellenberg & Jonas Ries (2019). Nuclear pores as versatile reference standards for quantitative superresolution microscopy . BioStudies, S-BIAD8. Retrieved from https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD8")
print("########")
resp = requests.get(AiryScan_example)
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
output_prefix = os.path.join(experiment.output_directory, experiment.date_as_string + 'NPC_quantification_AiryScan')

# Get a crop of the experimental image to reduce the field of view
experimental_img = np.array(img_tif, dtype=int)
halfpatch = 75 # same field of view as STED (6 µm)
center = int(experimental_img.shape[0]/2)
experimental_img_processed = experimental_img[center-halfpatch:center+halfpatch:,center-halfpatch:center+halfpatch]
# parameters from experimental image metadata
background = experimental_img_processed.min() # background level
pixelsize = 40 # in nm
min_distance = 1*pixelsize # in physical size, in nm
min_distance_px = int(np.ceil(min_distance / pixelsize))

# Calibration of the spot detection: number of local maxima for each sigma and threshold
calibration_sigmas = [0.5, 1, 1.5, 2]
calibration_thresholds = np.arange(0, 425, 25)
calibration = []
for s in calibration_sigmas:
    for t in calibration_thresholds:
        xy, _ = metrics.local_maxima_positions(
            experimental_img_processed, min_distance=min_distance_px,
            threshold=t, background=background, sigma=s)
        calibration.append(dict(sigma=s, threshold=t, n_maxima=len(xy)))
calibration = pd.DataFrame(calibration)
print("Number of local maxima (rows: threshold, columns: sigma)")
print(calibration.pivot(index="threshold", columns="sigma", values="n_maxima"))
calibration.to_csv(output_prefix + '_localmaxima_calibration.csv', index=False)

# parameters for spot detection
# for sigma = 1 the number of maxima is flat for thresholds between 100 and 200
# (noise maxima removed, NPC maxima kept): threshold in the middle of the plateau
sigma = 1 # in pixels
threshold = 150

fig, ax = plt.subplots(figsize=(5, 4))
sns.lineplot(data=calibration, x="threshold", y="n_maxima", hue="sigma", marker="o", palette="viridis", ax=ax)
chosen = calibration[(calibration.sigma == sigma) & (calibration.threshold == threshold)]
ax.scatter(chosen.threshold, chosen.n_maxima, s=150, facecolors="none", edgecolors="red", zorder=3, label="chosen")
ax.set_xlabel("Threshold (counts above background)")
ax.set_ylabel("Number of local maxima")
ax.legend(title="sigma (px)")
fig.tight_layout()
fig.savefig(output_prefix + '_localmaxima_calibration.png', dpi=300, bbox_inches='tight')
plt.close(fig)

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
## All other parameters are estimations to emulate Thevathasan et al. 2019 AiryScan modality
experiment.update_modality(
    modality_name="AiryScan",
    lateral_resolution_nm=60,
    pixelsize_nm=40)
experiment.set_modality_acq(
    modality_name="AiryScan",
    exp_time=0.01)
experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["binomial"]["p"] = 0.8
experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["baselevel"]["bl"] = 600
experiment.imager.modalities["AiryScan"]["detector"]["noise_model"]["gaussian"]["sigma"] = 25
images, noiselsess = experiment.run_simulation()

# Analyse images: fit circles
# expected parameters of objects
expected_min_radius = 40  # in nanometers
expected_MAX_radius = 70  # in nanometers

if images["AiryScan"]["ch0"][0].min() < 0:
    images["AiryScan"]["ch0"][0] += -images["AiryScan"]["ch0"][0].min()
simulated_img = images["AiryScan"]["ch0"][0]

def to_uint8(img):
    # scale the full intensity range to 0-255; casting 16-bit data to uint8
    # directly wraps values modulo 256
    img = np.asarray(img, dtype=float)
    img = img - img.min()
    return np.uint8(img / img.max() * 255) if img.max() > 0 else np.uint8(img)

circle_fit_failures = []

def fit_radii(img, pixelsize, image_label, method):
    """Radii (nm) of the circles found in img; [] if none, and the settings are
    recorded in circle_fit_failures."""
    minRadius_round = np.ceil(expected_min_radius / pixelsize).astype('int64')
    maxRadius_round = np.ceil(expected_MAX_radius / pixelsize).astype('int64')
    HCparams = dict(dp=1, minDist=maxRadius_round,
            param1=10, param2=7, minRadius=minRadius_round, maxRadius=maxRadius_round)
    circles, img_blurred, c_params = metrics.get_circles(to_uint8(img), **HCparams)
    try:
        return [r*pixelsize for (x, y, r) in circles[0]]
    except TypeError:  # get_circles returns None when no circle is found
        circle_fit_failures.append(dict(image=image_label, method=method, pixelsize_nm=pixelsize, **HCparams))
        return []

#### Circle fit at the native pixel size, and on images upsampled to the STED pixel size
upsampled_pixelsize = 15 # in nm
native_method = f"native {pixelsize} nm/px"
upsampled_method = f"upsampled {upsampled_pixelsize} nm/px"
radii = {}
for image_label, img in [("Experimental", experimental_img_processed), ("Simulated", simulated_img)]:
    radii[(image_label, native_method)] = fit_radii(img, pixelsize, image_label, native_method)
    img_upsampled = ndimage.zoom(np.asarray(img, dtype=float), pixelsize / upsampled_pixelsize, order=3)
    radii[(image_label, upsampled_method)] = fit_radii(img_upsampled, upsampled_pixelsize, image_label, upsampled_method)

if circle_fit_failures:
    circle_fit_failures = pd.DataFrame(circle_fit_failures)
    print("No circles found with these parameters:")
    print(circle_fit_failures.to_string(index=False))
    circle_fit_failures.to_csv(output_prefix + '_circle_fit_failures.csv', index=False)

#### Apparent blob size: FWHM of the projections of a window around each local maximum
# window: ring radius + 2 PSF standard deviations
blob_half_window_nm = 200

def blob_fwhms(img, peaks, pixelsize, half_nm=blob_half_window_nm):
    """FWHM (nm) of each blob: mean of the FWHMs of the projections onto x and y of the
    window of half_nm around the peak. Blobs at the image edge or not measurable are skipped."""
    img = np.asarray(img, dtype=float)
    half = int(round(half_nm / pixelsize))
    widths = []
    for row, col in peaks:
        if row - half < 0 or col - half < 0 or row + half + 1 > img.shape[0] or col + half + 1 > img.shape[1]:
            continue
        win = img[row-half:row+half+1, col-half:col+half+1]
        blob_widths = []
        for profile in (win.sum(axis=0), win.sum(axis=1)):
            x = (np.arange(profile.size) + 0.5) * pixelsize
            found = particle_measures._fwhm_1d(profile, x, 1, min_snr=0, mode="single")
            if found is None or found[0] is None:
                break
            blob_widths.append(found[0])
        if len(blob_widths) == 2:
            widths.append(float(np.mean(blob_widths)))
    return widths

peaks_exp, img_pre_exp = metrics.local_maxima_positions(
    experimental_img_processed, min_distance=min_distance_px,
    threshold=threshold, background=background, sigma=sigma)
# simulated intensities are in other units: same threshold relative to the maximum
peaks_sim, img_pre_sim = metrics.local_maxima_positions(
    simulated_img, min_distance=min_distance_px,
    threshold=None, background=simulated_img.min(), sigma=sigma,
    threshold_rel=threshold / img_pre_exp.max())
fwhm_experimental = blob_fwhms(experimental_img_processed, peaks_exp, pixelsize)
fwhm_simulated = blob_fwhms(simulated_img, peaks_sim, pixelsize)

# Measurements table and summary
measurements = [dict(image=image_label, measure="circle radius", method=method, value_nm=r)
                for (image_label, method), values in radii.items() for r in values]
measurements += [dict(image="Experimental", measure="blob FWHM", method="projection", value_nm=w) for w in fwhm_experimental]
measurements += [dict(image="Simulated", measure="blob FWHM", method="projection", value_nm=w) for w in fwhm_simulated]
measurements = pd.DataFrame(measurements, columns=["image", "measure", "method", "value_nm"])
measurements.to_csv(output_prefix + '_measurements.csv', index=False)
print(f"Local maxima: experimental {len(peaks_exp)}, simulated {len(peaks_sim)}")
print(measurements.groupby(["measure", "method", "image"]).value_nm.describe()[["count", "25%", "50%", "75%"]])

# Plot experimental vs simulated image

plt.rcParams['figure.figsize'] = [20, 5]

fig, axs = plt.subplots(1, 4)
axs[0].imshow(experimental_img_processed, cmap="grey")
axs[0].scatter(peaks_exp[:, 1], peaks_exp[:, 0], s=12, facecolors="none", edgecolors="red", linewidths=0.6)
axs[0].set_title("Modified from Thevathasan et al. 2019. Nature Methods.")
axs[1].imshow(simulated_img, cmap="grey")
axs[1].set_title("VLab4Mic simulated image")

df = measurements[measurements.measure == "circle radius"].copy()
df['Condition'] = df.image + " (" + df.method + ")"
if len(df):
    sns.histplot(data=df, x="value_nm", hue="Condition", binrange=[45,70], bins=15, kde=True, ax=axs[2])
if isinstance(circle_fit_failures, pd.DataFrame):
    failed = "\n".join(circle_fit_failures.image + " (" + circle_fit_failures.method + ")")
    axs[2].text(0.02, 0.98, "No circles found:\n" + failed, transform=axs[2].transAxes,
                va="top", ha="left", fontsize=8)
axs[2].set_xlabel("Radius of circle fit (nm)")

df = measurements[measurements.measure == "blob FWHM"].rename(columns={"image": "Condition"})
sns.histplot(data=df, x="value_nm", hue="Condition", bins=20, kde=True, ax=axs[3])
axs[3].set_xlabel("Apparent blob size, FWHM (nm)")

length_nm = 1000
nm = 1e-09
pixelsize = (experiment.imaging_modalities["AiryScan"]["detector"]["scale"] / nm) * experiment.imaging_modalities["AiryScan"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = length_px / 10


fontprops = fm.FontProperties(size=20)
scalebar = AnchoredSizeBar(axs[1].transData,
                           length_px, '1 µm', 'lower right',
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=2,
                           fontproperties=fontprops)
axs[1].add_artist(scalebar)
# remove axis ticks
axs[0].set_xticks([])
axs[0].set_yticks([])
axs[1].set_xticks([])
axs[1].set_yticks([])
plt.tight_layout()

# save figure
filename = output_prefix + '.png'
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
