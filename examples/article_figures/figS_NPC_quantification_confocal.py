from vlab4mic import experiments
from vlab4mic.analysis import metrics
import pandas as pd
import seaborn as sns
import requests
import tifffile as tif
import io
import os
import numpy as np
from scipy import ndimage
from scipy.stats import norm
from scipy.optimize import curve_fit
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
random_seed= 24

# parameters for simulation
structure = "7R5K"
probe_template = "GFP_w_nanobody"
probe_target_type = "Sequence"
probe_target_value = "ELAVGSL" # Sequence in the C-terminal of Nup96
modalities = ["Confocal", ]

# fetch experimental image from URL
Confocal_example = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/S-BIAD0-99/S-BIAD8/Files/Library/Gallery%20Fig%201/Confocal/GFP_Gallery-Confocal_180917_1.tif"
print("Downloading experimental image from Thevathasan et al. 2019")
print("#### If you find this image useful, please consider citing: ####")
print("Jervis Vermal Thevathasan, Maurice Kahnwald, Konstanty Cieśliński, Philipp Hoess, Sudheer Kumar Peneti, Manuel Reitberger, Daniel Heid, Krishna Chaitanya Kasuba, Sarah Janice Hoerner, Yiming Li, Yu-Le Wu, Markus Mund, Ulf Matti, Pedro Matos Pereira, Ricardo Henriques, Bianca Nijmeijer, Moritz Kueblbeck, Vilma Jimenez Sabinina, Jan Ellenberg & Jonas Ries (2019). Nuclear pores as versatile reference standards for quantitative superresolution microscopy . BioStudies, S-BIAD8. Retrieved from https://www.ebi.ac.uk/biostudies/BioImages/studies/S-BIAD8")
print("########")
resp = requests.get(Confocal_example)
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
results_directory = os.path.join(experiment.output_directory, experiment.date_as_string + 'NPC_quantification_Confocal')
os.makedirs(results_directory, exist_ok=True)
output_prefix = os.path.join(results_directory, experiment.date_as_string + 'NPC_quantification_Confocal')

# Get a crop of the experimental image to reduce the field of view
experimental_img = np.array(img_tif, dtype=int)
halfpatch = 43 # same field of view as STED (6 µm)
center = int(experimental_img.shape[0]/2)
experimental_img_processed = experimental_img[center-halfpatch:center+halfpatch:,center-halfpatch:center+halfpatch]
# parameters from experimental image metadata
background = experimental_img_processed.min() # background level
pixelsize = 70 # in nm
min_distance = 1*pixelsize # in physical size, in nm
min_distance_px = int(np.ceil(min_distance / pixelsize))

# parameters for spot detection: values of fig2_npc_panel.py (maxima at the image edges excluded)
sigma = 0.5 # in pixels
threshold = 400

# Analyse images: fit circles
# expected parameters of objects
expected_min_radius = 40  # in nanometers
expected_MAX_radius = 70  # in nanometers

def to_uint8(img):
    # scale the full intensity range to 0-255; casting 16-bit data to uint8
    # directly wraps values modulo 256
    img = np.asarray(img, dtype=float)
    img = img - img.min()
    return np.uint8(img / img.max() * 255) if img.max() > 0 else np.uint8(img)

circle_fit_failures = []

def fit_circles(img, pixelsize, image_label, method):
    """Circles (x, y, r in pixels of img) found in img; None if none, and the settings
    are recorded in circle_fit_failures."""
    minRadius_round = np.ceil(expected_min_radius / pixelsize).astype('int64')
    maxRadius_round = np.ceil(expected_MAX_radius / pixelsize).astype('int64')
    HCparams = dict(dp=1, minDist=maxRadius_round,
            param1=10, param2=7, minRadius=minRadius_round, maxRadius=maxRadius_round)
    failure = dict(image=image_label, method=method, pixelsize_nm=pixelsize, **HCparams)
    if maxRadius_round <= minRadius_round:
        # expected radii within one pixel: cv2.HoughCircles would replace maxRadius by
        # minRadius + 2 and return circles outside the expected range
        circle_fit_failures.append(dict(failure, reason="radius range within one pixel"))
        return None
    circles, _, _ = metrics.get_circles(to_uint8(img), **HCparams)
    try:
        return np.asarray(circles[0])
    except TypeError:  # get_circles returns None when no circle is found
        circle_fit_failures.append(dict(failure, reason="no circles found"))
        return None

def upsample(img):
    return ndimage.zoom(np.asarray(img, dtype=float), pixelsize / upsampled_pixelsize, order=3)

#### Circle fit at the native pixel size, and on images upsampled to the STED pixel size
upsampled_pixelsize = 15 # in nm
native_method = f"native {pixelsize} nm/px"
upsampled_method = f"upsampled {upsampled_pixelsize} nm/px"
radii = {}

def fit_both(img, image_label):
    """Radii (nm) of the circles of the native and upsampled fits, stored in radii."""
    for method, fit_img, fit_pixelsize in [(native_method, img, pixelsize),
                                           (upsampled_method, upsample(img), upsampled_pixelsize)]:
        circles = fit_circles(fit_img, fit_pixelsize, image_label, method)
        radii[(image_label, method)] = [] if circles is None else list(circles[:, 2] * fit_pixelsize)

#### Local maxima of the experimental image: positions of the NPCs
peaks_exp, _ = metrics.local_maxima_positions(
    experimental_img_processed, min_distance=min_distance_px,
    threshold=threshold, background=background, sigma=sigma)
print(f"{len(peaks_exp)} NPCs placed at the local maxima of the experimental image")

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
## All other parameters are estimations to emulate Thevathasan et al. 2019 Confocal modality
## Point-scanning detector: exp_time is the pixel dwell time. 3 µs gives a signal-to-noise
## ratio of the NPCs (29) close to the experimental one (33)
experiment.update_modality(
    modality_name="Confocal",
    lateral_resolution_nm=100,
    pixelsize_nm=70,
    scanning=True)
experiment.set_modality_acq(
    modality_name="Confocal",
    exp_time=3e-6)
images, _ = experiment.run_simulation()

simulated_img = images["Confocal"]["ch0"][0]

#### Circle fit of the experimental and simulated images
fit_both(experimental_img_processed, "Experimental")
fit_both(simulated_img, "Simulated")

if circle_fit_failures:
    circle_fit_failures = pd.DataFrame(circle_fit_failures)
    print("No circles found with these parameters:")
    print(circle_fit_failures.to_string(index=False))
    circle_fit_failures.to_csv(output_prefix + '_circle_fit_failures.csv', index=False)

#### Apparent blob size: FWHM of a 2D Gaussian fitted around each local maximum
# window: ring radius + 2 PSF standard deviations
blob_half_window_nm = 250

def gaussian_2d(xy, amplitude, x0, y0, s, offset):
    x, y = xy
    return (amplitude * np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * s**2)) + offset).ravel()

def blob_gaussian_fwhms(img, peaks, pixelsize, half_nm=blob_half_window_nm):
    """FWHM (nm) of each blob from a fit of a 2D Gaussian (isotropic, with offset) to the
    window of half_nm around the peak. Blobs at the image edge, fits that do not converge,
    and fits centred more than 2 pixels from the peak or wider than the window are skipped."""
    img = np.asarray(img, dtype=float)
    half = int(round(half_nm / pixelsize))
    y, x = np.mgrid[-half:half+1, -half:half+1]
    widths = []
    for row, col in peaks:
        if row - half < 0 or col - half < 0 or row + half + 1 > img.shape[0] or col + half + 1 > img.shape[1]:
            continue
        win = img[row-half:row+half+1, col-half:col+half+1]
        p0 = [win.max() - win.min(), 0, 0, 1.5, win.min()]
        try:
            (amplitude, x0, y0, s, offset), _ = curve_fit(gaussian_2d, (x, y), win.ravel(), p0=p0, maxfev=2000)
        except RuntimeError:
            continue
        if amplitude > 0 and abs(x0) < 2 and abs(y0) < 2 and 0.3 < abs(s) < half:
            widths.append(2 * np.sqrt(2 * np.log(2)) * abs(s) * pixelsize)
    return widths

# simulated image: blobs at the positions where the NPCs were placed
gaussian_fwhm_experimental = blob_gaussian_fwhms(experimental_img_processed, peaks_exp, pixelsize)
gaussian_fwhm_simulated = blob_gaussian_fwhms(simulated_img, peaks_exp, pixelsize)

def bimodal(values):
    """Two-component Gaussian mixture of values, and the BIC difference to a single Gaussian."""
    values = np.asarray(values).reshape(-1, 1)
    gmm = {k: GaussianMixture(n_components=k, random_state=random_seed, n_init=5).fit(values) for k in (1, 2)}
    delta_bic = gmm[2].bic(values) - gmm[1].bic(values)
    order = np.argsort(gmm[2].means_.ravel())
    return pd.DataFrame(dict(
        component=[1, 2],
        mean_nm=gmm[2].means_.ravel()[order],
        sd_nm=np.sqrt(gmm[2].covariances_.ravel())[order],
        weight=gmm[2].weights_[order],
        delta_bic_2_minus_1=delta_bic))

# Bimodal fit of the experimental and simulated blob size, compared to a single Gaussian
gaussian_fwhms = {"Experimental": gaussian_fwhm_experimental, "Simulated": gaussian_fwhm_simulated}
bimodal_fits = {image_label: bimodal(values) for image_label, values in gaussian_fwhms.items()}
bimodal_fit = pd.concat(bimodal_fits, names=["image"]).reset_index(level=0)
bimodal_fit.to_csv(output_prefix + '_fwhm_bimodal_fit.csv', index=False)
print("Bimodal fit of the blob FWHM")
print(bimodal_fit.round(2).to_string(index=False))

# Measurements table and summary
measurements = [dict(image=image_label, measure="circle radius", method=method, value_nm=r)
                for (image_label, method), values in radii.items() for r in values]
for image_label, values in [("Experimental", gaussian_fwhm_experimental), ("Simulated", gaussian_fwhm_simulated)]:
    measurements += [dict(image=image_label, measure="blob FWHM", method="gaussian fit", value_nm=w) for w in values]
measurements = pd.DataFrame(measurements, columns=["image", "measure", "method", "value_nm"])
measurements.to_csv(output_prefix + '_measurements.csv', index=False)
print(measurements.groupby(["measure", "method", "image"]).value_nm.describe()[["count", "25%", "50%", "75%"]])

# Plot experimental vs simulated image

plt.rcParams['figure.figsize'] = [20, 5]

fig, axs = plt.subplots(1, 4)
axs[0].imshow(experimental_img_processed, cmap="grey")
axs[0].set_title("Modified from Thevathasan et al. 2019. Nature Methods.")
axs[1].imshow(simulated_img, cmap="grey")
axs[1].set_title("VLab4Mic simulated image")

df = measurements[measurements.measure == "circle radius"].copy()
df['Condition'] = df.image + " (" + df.method + ")"
if len(df):
    sns.histplot(data=df, x="value_nm", hue="Condition", binrange=[45,70], bins=15, kde=True, ax=axs[2])
    sns.move_legend(axs[2], "center left")
axs[2].text(0.02, 0.98, f"Expected radius: {expected_min_radius}–{expected_MAX_radius} nm",
            transform=axs[2].transAxes, va="top", ha="left")
axs[2].set_xlabel("Radius of circle fit (nm)")

df = measurements[measurements.measure == "blob FWHM"].rename(columns={"image": "Condition"})
bin_edges = np.histogram_bin_edges(df.value_nm, bins=20)
condition_colours = dict(zip(gaussian_fwhms, sns.color_palette(n_colors=len(gaussian_fwhms))))
sns.histplot(data=df, x="value_nm", hue="Condition", bins=bin_edges, palette=condition_colours, ax=axs[3])
x = np.linspace(bin_edges[0], bin_edges[-1], 300)
fit_text = []
for image_label, fit in bimodal_fits.items():
    colour = sns.set_hls_values(condition_colours[image_label], l=0.3)  # darker than the bars
    scale = len(gaussian_fwhms[image_label]) * (bin_edges[1] - bin_edges[0])
    components = [c.weight * norm.pdf(x, c.mean_nm, c.sd_nm) * scale for c in fit.itertuples()]
    for component in components:
        axs[3].plot(x, component, color=colour, linestyle="--", linewidth=1)
    axs[3].plot(x, np.sum(components, axis=0), color=colour, linewidth=1.5)
    fit_text.append(f"{image_label}, bimodal fit:\n"
                    + "\n".join(f"μ = {c.mean_nm:.0f} nm, σ = {c.sd_nm:.0f} nm ({c.weight:.0%})" for c in fit.itertuples())
                    + f"\nΔBIC (2 − 1 components) = {fit.delta_bic_2_minus_1.iloc[0]:+.1f}")
axs[3].text(0.98, 0.72, "\n\n".join(fit_text), transform=axs[3].transAxes, va="top", ha="right", fontsize=8)
axs[3].set_title("Blob size: 2D Gaussian fit")
axs[3].set_xlabel("Apparent blob size, FWHM (nm)")

length_px = 1000 / pixelsize  # 1 µm
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

# Intensity histograms of the two images of the figure. Each image minus its background
# (5th percentile of the image smoothed with a 1 px Gaussian), drawn up to its 99.9th
# percentile: a few brighter pixels (e.g. repeated blinking in SMLM) would stretch the axis.
# Last panel: both images divided by their 99.9th percentile, to compare the shapes
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
normalised = {}
for ax, img, title in [(axs[0], experimental_img_processed, "Experimental (crop)"),
                       (axs[1], simulated_img, "VLab4Mic simulated")]:
    img = np.asarray(img, dtype=float)
    img_background = np.percentile(ndimage.gaussian_filter(img, 1), 5)
    values = np.ravel(img) - img_background
    upper = np.percentile(values, 99.9)
    ax.hist(values, bins=100, range=(values.min(), upper))
    ax.set_yscale("log")
    ax.set_title(f"{title}\n{np.sum(values > upper)} pixels above {upper:.0f} not shown (max {values.max():.0f})")
    ax.set_xlabel("Intensity above background (counts)")
    ax.set_ylabel("Number of pixels")
    print(f"{title}: background {img_background:.0f}, 99.9th percentile {upper:.0f}, "
          f"max {values.max():.0f} counts above background")
    normalised[title] = values / upper
bin_edges = np.linspace(min(v.min() for v in normalised.values()), 1, 101)
for title, values in normalised.items():
    axs[2].hist(values, bins=bin_edges, density=True, histtype="step", linewidth=1.5, label=title)
axs[2].set_yscale("log")
axs[2].set_title("Normalised to the 99.9th percentile")
axs[2].set_xlabel("Intensity above background / 99.9th percentile")
axs[2].set_ylabel("Density")
axs[2].legend()
fig.tight_layout()
fig.savefig(output_prefix + '_intensity_histograms.png', dpi=300, bbox_inches='tight')
plt.close(fig)
