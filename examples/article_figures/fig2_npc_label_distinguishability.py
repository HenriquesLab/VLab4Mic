"""Fig. 2: how much does the linkage error of each probe enlarge the nuclear pore, and in which imaging method can it be seen?

Nup96 (7R5K, C-terminal motif ELAVGSL) is labelled with four probes of
increasing linkage error: a direct label (a Linker at 0 nm from the
epitope: no linkage error, the reference), a SNAP-tag, an anti-GFP
nanobody (GFP + nanobody) and an anti-GFP antibody (GFP fused to Nup96 as
primary, at labelling efficiency 1, and an Antibody bound to it), as in
figS2_npc_three_labels.py. The reference dataset (Thevathasan et al. 2019)
gives ring radii of 53.7 +/- 2.1, 55.0 +/- 1.9 and 64.3 +/- 2.6 nm for the
last three; the Nup96 C-terminus of the model is at EPITOPE_RADIUS_NM.

For each probe, N independent realisations (one pore each, at the centre
of a 1000 x 1000 nm field, seen from above and rotated in plane) are
imaged by the five imaging methods of the dataset: Widefield, Confocal,
AiryScan, STED and SMLM. Every realisation is imaged by all of them, so
within a probe the methods see the same pores. The modality parameters
are fixed here (nothing is downloaded): resolution and noise of
fig2_npc_panel.py for Widefield, Confocal and AiryScan (with the default
Widefield resolution), STED as in figS2_npc_three_labels.py and SMLM with
its defaults. The exposure of each probe is set so that a pore gives
SIGNAL_PHOTONS in each camera modality, so all probes are equally bright
despite their different numbers of fluorophores.

Every image is measured from its pixels (as in
fig2_npc_integrity_distinguishability.py):
- radius_nm, width_nm: the pixels within REGION_NM of the image centre
  above an Otsu threshold are taken as points of the ring; radius of a
  geometric circle fit to them and s.d. of their distances to its centre.
  In Widefield, Confocal and AiryScan the ring is not resolved and this is
  the apparent radius of the blob.
- blob_fwhm_nm (Widefield, Confocal and AiryScan only): FWHM of the blob,
  the mean of the FWHMs of the projections of the image onto x and y, with
  the background taken BACKGROUND_NM at each end of the projection.

Within each imaging method, the increase of the mean radius (and blob
FWHM) from the direct label is the enlargement due to the linkage error
of each probe, with a 95% interval; the published increase is the
published radius minus EPITOPE_RADIUS_NM. Every pair of probes is scored
per imaging method with distinguishability_from_features: out-of-fold ROC
AUC with a 95% bootstrap interval, a null control (two random halves of
the first probe, expected 0.5) and a permutation p-value, for each feature
set (radius alone, radius and width, blob FWHM alone).

Outputs, in the folder <date>fig2_npc_labels of the VLab4Mic output
directory: npc_label_figure.png, the main figure: (a) one image per
imaging method (rows) and probe (columns), with its fitted circle, (b) the
increase of the radius from the direct label per probe, one line per
imaging method, with the published increase, (c) the AUC of each probe
against the direct label from the radius alone, (d) the increase of the
blob FWHM in Widefield, Confocal and AiryScan; features.csv (one row per
image), increase_from_direct.csv, distinguishability.csv (one row per pair
of probes, feature set and imaging method), npc_label_figure_values.csv
(the values of b, c and d, one row per measure, imaging method and probe),
settings.json and examples.npz (the example images of a and their fitted
circles, so --features plots them again).

Run:   python fig2_npc_label_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_npc_label_distinguishability.py --features <folder>/features.csv
(the outputs are written next to the feature table)

Environment variables: VLAB4MIC_N_REALISATIONS (100 per probe; e.g. 5 for
a quick run), VLAB4MIC_N_PERMUTATIONS (200), VLAB4MIC_N_WORKERS (1;
probes and scoring in parallel processes), VLAB4MIC_SAVE_IMAGES (0; 1
saves the images of every probe and imaging method as a TIFF stack, one
realisation per frame, in <folder>/images).
"""

import itertools
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import to_hex
from tqdm.auto import tqdm

from vlab4mic import experiments
from vlab4mic.analysis import particle_measures
from vlab4mic.analysis.distinguishability import distinguishability_from_features
from vlab4mic.experiments import ExperimentParametrisation, output_path

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 44
OUTPUT_NAME = "fig2_npc_labels"  # base name of the output folder (after the date)

STRUCTURE = "7R5K"
LABELLING_EFFICIENCY = 0.4  # Supp. Note 4
# 1000 nm (the detector field of view of every modality): the Widefield blob
# (FWHM ~320 nm) leaves background at the edges to measure its size
SAMPLE_DIMENSIONS = [1000, 1000, 100]

TARGET = dict(probe_target_type="Sequence", probe_target_value="ELAVGSL")  # Nup96 C-terminus
# probes in order of linkage error; the first one is the reference
LABELS = {
    "direct": [dict(probe_template="Linker", probe_name="direct", probe_distance_to_epitope=0,
                    **TARGET)],
    "SNAP-tag": [dict(probe_template="SNAP-tag", probe_name="SNAP-tag", **TARGET)],
    "anti-GFP nanobody": [dict(probe_template="GFP_w_nanobody", probe_name="GFP_w_nanobody",
                               probe_DoL=2, **TARGET)],
    "anti-GFP antibody": [
        dict(probe_template="GFP", probe_name="Nup96-GFP", as_primary=True,
             labelling_efficiency=1.0, **TARGET),
        dict(probe_template="Antibody", probe_name="anti-GFP",
             probe_target_type="Primary", probe_target_value="Nup96-GFP"),
    ],
}
REFERENCE = "direct"
PROBE_NAMES = {"direct": "direct", "SNAP-tag": "SNAP-tag", "anti-GFP nanobody": "nanobody",
               "anti-GFP antibody": "antibody"}
PUBLISHED_RADIUS_NM = {
    "SNAP-tag": (53.7, 2.1),
    "anti-GFP nanobody": (55.0, 1.9),
    "anti-GFP antibody": (64.3, 2.6),
}
EPITOPE_RADIUS_NM = 53.5  # Nup96 C-terminus of 7R5K (pore diameter 107 nm at the epitopes)

MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]  # increasing resolution
DIFFRACTION_LIMITED = ["Widefield", "Confocal", "AiryScan"]  # the ring is not resolved
CAMERA_MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED"]  # exposure set per probe
# resolution: fig2_npc_panel.py (Widefield: default) and figS2_npc_three_labels.py (STED)
MODALITY_PARAMS = {
    "Confocal": dict(lateral_resolution_nm=100),
    "AiryScan": dict(lateral_resolution_nm=60),
    "STED": dict(lateral_resolution_nm=15, axial_resolution_nm=20, psf_voxel_nm=5,
                 depth_of_field_nm=100, pixelsize_nm=15),
}
# detector noise, same sources: binomial p, base level and s.d. of the Gaussian (read) noise
NOISE = {
    "Widefield": dict(binomial=0.8, baselevel=4000, gaussian=4000 ** 0.5),
    "AiryScan": dict(binomial=0.8, baselevel=600, gaussian=25),
    "STED": dict(binomial=0.7, gaussian=0.2),
}
NOISE_KEYS = dict(binomial="p", baselevel="bl", gaussian="sigma")
# exposure of each probe: a pore gives SIGNAL_PHOTONS (noiseless image), estimated from
# CALIBRATION_REALISATIONS realisations at EXP_TIME. EXP_TIME is the exposure of
# fig2_npc_panel.py (STED: figS2_npc_three_labels.py), at which a nanobody pore gives about
# SIGNAL_PHOTONS (STED: STED_SIGNAL_PHOTONS of figS2_npc_three_labels.py)
EXP_TIME = {"Widefield": 0.01, "Confocal": 0.001, "AiryScan": 0.01, "STED": 0.0002}
SIGNAL_PHOTONS = {"Widefield": 17000, "Confocal": 1900, "AiryScan": 20000, "STED": 2000}
CALIBRATION_REALISATIONS = 10

# measurement of the pore: pixels further than REGION_NM from the image centre are background
# (the ring blurred by the PSF: about the ring radius + 2 PSF s.d.)
REGION_NM = {"Widefield": 300, "Confocal": 250, "AiryScan": 200, "STED": 100, "SMLM": 100}
BACKGROUND_NM = 200  # at each end of a projection, for the blob FWHM

FEATURES = ["radius_nm", "width_nm"]
FEATURE_SETS = {"radius": ["radius_nm"], "all": FEATURES, "blob size": ["blob_fwhm_nm"]}
CONTROLS = {feature: feature for feature in FEATURES}  # each alone
# images without a measurable ring (e.g. too few labels): nothing detected
FILL_VALUES = dict(radius_nm=0, width_nm=0)
PAIRS = list(itertools.combinations(LABELS, 2))

# figure
CROP_NM = 400  # example images: around the image centre
SCALEBAR_NM = 100
MODALITY_COLOURS = dict(zip(MODALITIES, [to_hex(colour) for colour in
                                         plt.cm.viridis(np.linspace(0, 0.9, len(MODALITIES)))]))


def pixelsize_nm(experiment, modality):
    detector = experiment.imaging_modalities[modality]["detector"]
    return detector["pixelsize"] * detector["scale"] / 1e-9


def set_acquisition(experiment, exp_times):
    """Exposure time of the camera modalities and detector noise."""
    for modality, exp_time in exp_times.items():
        experiment.set_modality_acq(modality_name=modality, exp_time=exp_time)
    for modality, noise in NOISE.items():
        noise_model = experiment.imager.modalities[modality]["detector"]["noise_model"]
        for name, value in noise.items():
            noise_model[name][NOISE_KEYS[name]] = value


def measure_npc(image, pixelsize, region_nm):
    """Ring radius and width from a circle fit to the bright pixels around the image centre,
    and the centre of the circle (nm, from the image centre)."""
    missing = dict(radius_nm=np.nan, width_nm=np.nan, n_pixels=0, centre_x_nm=np.nan,
                   centre_y_nm=np.nan)
    image = np.asarray(image, dtype=float)
    rows_px, cols_px = np.indices(image.shape)
    xy = np.c_[cols_px.ravel() + 0.5, rows_px.ravel() + 0.5] * pixelsize  # pixel centres, nm
    values = image.ravel()
    image_centre = np.array(image.shape[::-1]) * pixelsize / 2
    region = np.linalg.norm(xy - image_centre, axis=1) <= region_nm
    xy, values = xy[region], values[region]
    if np.ptp(values) == 0:
        return missing
    points = xy[values > particle_measures.compute_threshold(values, method="otsu")]
    if len(points) < 3:
        return dict(missing, n_pixels=len(points))
    centre, radius = particle_measures.fit_circle_geometric(points)
    width = np.linalg.norm(points - centre, axis=1).std()
    offset = np.asarray(centre) - image_centre
    return dict(radius_nm=radius, width_nm=width, n_pixels=len(points),
                centre_x_nm=offset[0], centre_y_nm=offset[1])


def blob_fwhm(image, pixelsize):
    """FWHM of the blob (nm): mean of the FWHMs of the projections of the image onto x and y,
    each with its background taken BACKGROUND_NM at both ends. NaN if one is not measurable."""
    image = np.asarray(image, dtype=float)
    background_px = max(1, int(round(BACKGROUND_NM / pixelsize)))
    widths = []
    for profile in (image.sum(axis=0), image.sum(axis=1)):
        x = (np.arange(profile.size) + 0.5) * pixelsize
        found = particle_measures._fwhm_1d(profile, x, background_px, min_snr=0, mode="single")
        if found is None or found[0] is None:
            return np.nan
        widths.append(found[0])
    return float(np.mean(widths))


def crop_centre(image, pixelsize):
    """CROP_NM around the image centre, and its extent (left, right, bottom, top) in nm from
    the image centre, for imshow with the origin at the top."""
    image = np.asarray(image, dtype=float)
    half = int(round(CROP_NM / 2 / pixelsize))
    rows, cols = np.array(image.shape) // 2
    crop = image[rows - half:rows + half, cols - half:cols + half]
    centre_x, centre_y = np.array(image.shape[::-1]) * pixelsize / 2
    extent = np.array([(cols - half) * pixelsize - centre_x, (cols + half) * pixelsize - centre_x,
                       (rows + half) * pixelsize - centre_y, (rows - half) * pixelsize - centre_y])
    return crop.astype(np.float32), extent


def example_key(label, modality, what):
    return f"{label}_{modality}_{what}".replace(" ", "_")


def simulate(job):
    """One probe: N_REALISATIONS pores, each imaged by every imaging method and measured (one
    row per image), and the first image of every method with a fitted ring (cropped)."""
    label, seed, images_dir = job
    __, __, experiment = experiments.image_vsample(
        structure=STRUCTURE,
        probe_list=[dict({"labelling_efficiency": LABELLING_EFFICIENCY}, **probe)
                    for probe in LABELS[label]],
        multimodal=MODALITIES,
        number_of_particles=1,
        sample_dimensions=SAMPLE_DIMENSIONS,
        sample_inital_orientation=[0, 0, 1],
        run_simulation=False,
        clear_experiment=True,
        random_seed=seed,
    )
    experiment.set_virtualsample_params(random_orientations=False, random_rotations=True)
    for modality, params in MODALITY_PARAMS.items():
        experiment.update_modality(modality_name=modality, **params)
    # equal brightness for all probes: scale the exposure to the probe's signal
    set_acquisition(experiment, EXP_TIME)
    calibration = experiment.run_replicates(CALIBRATION_REALISATIONS, modality="All")
    exp_times = {}
    for modality in CAMERA_MODALITIES:
        signal = np.mean([np.asarray(result["images_noiseless"][modality]["ch0"])[0].sum()
                          for result in calibration])
        exp_times[modality] = EXP_TIME[modality] * SIGNAL_PHOTONS[modality] / signal
    set_acquisition(experiment, exp_times)
    pixelsizes = {modality: pixelsize_nm(experiment, modality) for modality in MODALITIES}

    rows, examples = [], {}
    stacks = {modality: [] for modality in MODALITIES}
    for result in experiment.run_replicates(N_REALISATIONS, modality="All"):
        for modality in MODALITIES:
            image = np.asarray(result["images"][modality]["ch0"])[0]
            pixelsize = pixelsizes[modality]
            measures = measure_npc(image, pixelsize, REGION_NM[modality])
            blob = blob_fwhm(image, pixelsize) if modality in DIFFRACTION_LIMITED else np.nan
            rows.append(dict(label=label, modality_name=modality, replicate=result["replicate"],
                             seed=seed, exp_time=exp_times.get(modality, np.nan),
                             pixelsize_nm=pixelsize, blob_fwhm_nm=blob, **measures))
            if example_key(label, modality, "image") not in examples \
                    and np.isfinite(measures["radius_nm"]):
                crop, extent = crop_centre(image, pixelsize)
                examples[example_key(label, modality, "image")] = crop
                examples[example_key(label, modality, "extent")] = extent
                examples[example_key(label, modality, "circle")] = np.array(
                    [measures["centre_x_nm"], measures["centre_y_nm"], measures["radius_nm"]])
            if images_dir is not None:
                stacks[modality].append(image)
    if images_dir is not None:
        os.makedirs(images_dir, exist_ok=True)
        for modality, stack in stacks.items():
            tifffile.imwrite(os.path.join(images_dir, f"{label}_{modality}.tif".replace(" ", "_")),
                             np.array(stack, dtype=np.float32))
    return pd.DataFrame(rows), examples


def increase_from_direct(features, column):
    """Per imaging method and probe: mean, s.d. and number of the measured images, and the
    increase of the mean from the direct label with a 95% interval; with the published radius
    and its increase from EPITOPE_RADIUS_NM (0 for the direct label, at the epitope)."""
    data = features.dropna(subset=[column])
    rows = []
    for modality in MODALITIES:
        summary = data[data.modality_name == modality].groupby("label")[column] \
            .agg(["mean", "std", "size"])
        if REFERENCE not in summary.index:
            continue
        reference = summary.loc[REFERENCE]
        for label in LABELS:
            if label not in summary.index:
                continue
            row = summary.loc[label]
            increase = row["mean"] - reference["mean"]
            half = np.nan
            if label != REFERENCE:
                half = 1.96 * np.sqrt(row["std"] ** 2 / row["size"]
                                      + reference["std"] ** 2 / reference["size"])
            published, published_sd = PUBLISHED_RADIUS_NM.get(label, (np.nan, np.nan))
            rows.append(dict(
                metric=column, modality_name=modality, label=label, n=int(row["size"]),
                mean=row["mean"], sd=row["std"], reference_mean=reference["mean"],
                increase=increase, increase_low=increase - half, increase_high=increase + half,
                published_radius=published, published_sd=published_sd,
                published_increase=0.0 if label == REFERENCE else published - EPITOPE_RADIUS_NM))
    return pd.DataFrame(rows)


def score(features):
    """Distinguishability of every pair of probes per imaging method and feature set."""
    features = features.fillna(FILL_VALUES)
    tables = []
    for (a, b), (name, columns) in itertools.product(PAIRS, FEATURE_SETS.items()):
        print(f"{a} vs {b}, features: {name}")
        data = features[features.label.isin([a, b])]
        if name == "blob size":
            data = data[data.modality_name.isin(DIFFRACTION_LIMITED)].dropna(subset=columns)
        table = distinguishability_from_features(
            data, "label", (a, b), columns, by=["modality_name"],
            controls=CONTROLS if name == "all" else None,
            replicate_column="replicate", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        table.insert(0, "probe_b", b)
        table.insert(0, "probe_a", a)
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def figure_values(increase, scores):
    """The values of b, c and d, one row per measure (radius_nm: b and c; blob_fwhm_nm: d),
    imaging method and probe: the increase from the direct label and the AUC of the probe
    against the direct label from that measure alone (95% bootstrap interval, permutation
    p-value). For the direct label, null_auc is the mean of the null controls (two random
    halves of the direct label) of its comparisons."""
    feature_set = {"radius_nm": "radius", "blob_fwhm_nm": "blob size"}
    against = scores[scores.probe_a == REFERENCE]
    rows = []
    for __, row in increase.iterrows():
        cell = against[(against.feature_set == feature_set[row.metric])
                       & (against.modality_name == row.modality_name)]
        values = dict(auc=np.nan, auc_low=np.nan, auc_high=np.nan, p_value=np.nan,
                      null_auc=cell.null_auc.mean())
        found = cell[cell.probe_b == row.label]
        if row.label != REFERENCE and len(found):
            values = found.iloc[0][["auc", "auc_low", "auc_high", "null_auc", "p_value"]].to_dict()
        rows.append(dict(row.to_dict(), **values))
    return pd.DataFrame(rows)


def error_bars(value, low, high):
    """Lengths of the error bars below and above value, never negative (rounding can put
    value a hair outside [low, high], which matplotlib refuses); none where low or high is
    missing."""
    value, low, high = (np.asarray(x, dtype=float) for x in (value, low, high))
    low, high = np.where(np.isnan(low), value, low), np.where(np.isnan(high), value, high)
    return [np.clip(value - low, 0, None), np.clip(high - value, 0, None)]


def probe_axis(ax):
    """Probes in order of linkage error."""
    ax.set_xticks(range(len(LABELS)), [PROBE_NAMES[label] for label in LABELS], fontsize=8)
    ax.set_xlim(-0.4, len(LABELS) - 0.6)
    ax.spines[["top", "right"]].set_visible(False)


def plot_published(ax, values):
    """The published increase of the radius from the Nup96 C-terminus, black and dashed."""
    published = values.drop_duplicates("label").set_index("label").reindex(LABELS)
    ax.errorbar(range(len(LABELS)), published.published_increase,
                yerr=published.published_sd.fillna(0), fmt="s--", color="k", capsize=3,
                markersize=4, linewidth=1,
                label=f"published radius − {EPITOPE_RADIUS_NM:g} nm (epitopes)")


def plot_figure(examples, values, prefix):
    """The main figure: (a) one image per imaging method (rows) and probe (columns) with its
    fitted circle, (b) the increase of the radius from the direct label, (c) the AUC of each
    probe against the direct label from the radius alone and (d) the increase of the blob
    FWHM in the diffraction-limited methods, against the probe, one line per method."""
    labels = list(LABELS)
    x = np.arange(len(labels))
    dodge = dict(zip(MODALITIES, np.linspace(-0.12, 0.12, len(MODALITIES))))
    radius = values[values.metric == "radius_nm"]
    blob = values[values.metric == "blob_fwhm_nm"]
    fig = plt.figure(figsize=[11.5, 8.2])
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 1.1], wspace=0.28)

    left = grid[0].subgridspec(len(MODALITIES), len(labels), wspace=0.05, hspace=0.08)
    for i, modality in enumerate(MODALITIES):
        for k, label in enumerate(labels):
            ax = fig.add_subplot(left[i, k])
            ax.set_xticks([])
            ax.set_yticks([])
            image = examples.get(example_key(label, modality, "image"))
            extent = np.array([-CROP_NM / 2, CROP_NM / 2, CROP_NM / 2, -CROP_NM / 2])
            if image is not None:
                extent = examples[example_key(label, modality, "extent")]
                ax.imshow(image, cmap="magma", extent=extent)
                centre_x, centre_y, r = examples[example_key(label, modality, "circle")]
                ax.add_patch(plt.Circle((centre_x, centre_y), r, fill=False, linestyle=":",
                                        color="w", linewidth=1))
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            if i == 0:
                published = PUBLISHED_RADIUS_NM.get(label)
                detail = f"{published[0]} ± {published[1]} nm" if published else "0 nm linker"
                ax.set_title(f"{PROBE_NAMES[label]}\n{detail}", fontsize=8)
            if k == 0:
                ax.set_ylabel(modality, fontsize=8)
                ax.spines["left"].set_color(MODALITY_COLOURS[modality])
                ax.spines["left"].set_linewidth(4)
            if i == 0 and k == 0:
                ax.set_title("a", loc="left", fontweight="bold", x=-0.3)
            if i == len(MODALITIES) - 1 and k == 0:
                x0, y0 = extent[0] + 25, extent[2] - 30
                ax.plot([x0, x0 + SCALEBAR_NM], [y0, y0], "-", color="w", lw=2)

    right = grid[1].subgridspec(3, 1, hspace=0.5)
    ax_b = fig.add_subplot(right[0])
    for modality in MODALITIES:
        data = radius[radius.modality_name == modality].set_index("label").reindex(labels)
        ax_b.errorbar(x + dodge[modality], data.increase,
                      yerr=error_bars(data.increase, data.increase_low, data.increase_high),
                      fmt="o-", color=MODALITY_COLOURS[modality], capsize=2, markersize=4,
                      label=modality)
    plot_published(ax_b, radius)
    ax_b.axhline(0, color="0.5", linestyle=":", linewidth=1)
    ax_b.set_ylabel("Radius increase (nm)")
    ax_b.set_title("b", loc="left", fontweight="bold", x=-0.14)
    ax_b.set_title("Apparent ring radius, from the direct label", fontsize=9)
    ax_b.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))

    ax_c = fig.add_subplot(right[1])
    for j, modality in enumerate(MODALITIES):
        data = radius[radius.modality_name == modality].set_index("label").reindex(labels)
        colour = MODALITY_COLOURS[modality]
        ax_c.errorbar(x[1:] + dodge[modality], data.auc[1:],
                      yerr=error_bars(data.auc[1:], data.auc_low[1:], data.auc_high[1:]),
                      fmt="o-", color=colour, capsize=2, markersize=4)
        ax_c.plot(x[0] + dodge[modality], data.null_auc.iloc[0], "x", color=colour,
                  label="null (direct vs direct)" if j == 0 else None)
    ax_c.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_c.set_ylim(min(0.3, np.nanmin(radius.null_auc) - 0.05), 1.05)
    ax_c.set_ylabel("ROC AUC")
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.14)
    ax_c.set_title("Probe vs direct label, ring radius", fontsize=9)
    ax_c.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))

    ax_d = fig.add_subplot(right[2])
    for modality in DIFFRACTION_LIMITED:
        data = blob[blob.modality_name == modality].set_index("label").reindex(labels)
        if data.increase.isna().all():
            continue
        ax_d.errorbar(x + dodge[modality], data.increase,
                      yerr=error_bars(data.increase, data.increase_low, data.increase_high),
                      fmt="o-", color=MODALITY_COLOURS[modality], capsize=2, markersize=4,
                      label=f"{modality}: FWHM {data.reference_mean.iloc[0]:.0f} nm, "
                            f"antibody AUC {data.auc.iloc[-1]:.2f}")
    plot_published(ax_d, radius)
    ax_d.axhline(0, color="0.5", linestyle=":", linewidth=1)
    ax_d.set_ylabel("FWHM increase (nm)")
    ax_d.set_title("d", loc="left", fontweight="bold", x=-0.14)
    ax_d.set_title("Apparent blob size, from the direct label", fontsize=9)
    ax_d.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    for ax in (ax_b, ax_c, ax_d):
        probe_axis(ax)
    ax_d.set_xlabel("Probe (increasing linkage error)")
    fig.savefig(prefix + "npc_label_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
        examples = {}
        if os.path.exists(prefix + "examples.npz"):
            examples = dict(np.load(prefix + "examples.npz"))
    else:
        output_dir = os.path.join(str(output_path),
                                  ExperimentParametrisation().date_as_string + OUTPUT_NAME)
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, "")  # the outputs are named <output_dir>/<name>
        images_dir = prefix + "images" if SAVE_IMAGES else None
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(LABELS))
        jobs = [(label, int(seed), images_dir) for label, seed in zip(LABELS, seeds)]
        print(f"{len(jobs)} probes, {N_REALISATIONS} realisations each, imaged by "
              f"{len(MODALITIES)} methods, {N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(min(N_WORKERS, len(jobs))) as pool:
                outputs = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="probes"))
        else:
            outputs = [simulate(job) for job in tqdm(jobs, desc="probes")]
        features = pd.concat([table for table, __ in outputs], ignore_index=True)
        features.to_csv(prefix + "features.csv", index=False)
        examples = {key: value for __, found in outputs for key, value in found.items()}
        np.savez_compressed(prefix + "examples.npz", **examples)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    increase = pd.concat([
        increase_from_direct(features, "radius_nm"),
        increase_from_direct(features[features.modality_name.isin(DIFFRACTION_LIMITED)],
                             "blob_fwhm_nm"),
    ], ignore_index=True)
    increase.to_csv(prefix + "increase_from_direct.csv", index=False)
    scores = score(features)
    scores.to_csv(prefix + "distinguishability.csv", index=False)
    values = figure_values(increase, scores)
    values.to_csv(prefix + "npc_label_figure_values.csv", index=False)
    if examples:
        plot_figure(examples, values, prefix)
    print(values[["metric", "modality_name", "label", "n", "mean", "sd", "increase",
                  "published_increase", "auc", "null_auc"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
