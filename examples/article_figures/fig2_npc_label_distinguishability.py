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
directory: npc_label_figure.png, the main figure, in four sections:
1, the model: (a) the structure seen from above with the Nup96 C-termini
(epitopes), (b) the fluorophores of one pore per probe, with the
fluorophore ring radius of the model (emitter_radius_nm: a description of
the model, never used as a measurement); 2, SMLM, as in
figS2_npc_three_labels.py: (c) one image per probe with its fitted circle,
(d) histograms of the fitted radii, (e) simulated against published radius;
3, SMLM distinguishability, with the probe as a proxy for the linkage
error: (f) the AUC of every pair of probes from the radius alone, (g) the
radius increase from the direct label (the estimated linkage error) against
the published one; 4, the other imaging methods, as an extension: (h) one
image per method and probe, (i) the AUC of each probe against the direct
label, (j) the increase of the blob FWHM in Widefield, Confocal and
AiryScan. Also features.csv (one row per image), increase_from_direct.csv,
distinguishability.csv (one row per pair of probes, feature set and
imaging method), smlm_auc_matrix.csv (the values of f),
npc_label_figure_values.csv (the values of the other panels, one row per
measure, imaging method and probe), settings.json and examples.npz (the
structure, fluorophores and example images with their fitted circles, so
--features plots them again).

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
MODEL_VIEW_NM = 170  # model panels: width of the view around the pore centre
MODEL_SCALEBAR_NM = 50
N_STRUCTURE_ATOMS = 20000  # atoms of the structure shown (random subset)
OTHER_MODALITIES = [modality for modality in MODALITIES if modality != "SMLM"]
MODALITY_COLOURS = dict(zip(MODALITIES, [to_hex(colour) for colour in
                                         plt.cm.viridis(np.linspace(0, 0.9, len(MODALITIES)))]))
# probes in order of linkage error (checked for colour vision deficiency: every pair apart)
PROBE_COLOURS = dict(zip(LABELS, ["#3b3b3b", "#0072b2", "#e69f00", "#d55e00"]))


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


def structure_view(structure, target, seed):
    """N_STRUCTURE_ATOMS atoms of the structure and its epitopes (target), in nm from the
    centre of the assembly."""
    to_nm = structure.scale / 1e-9
    atoms = np.asarray(structure.assembly_atoms, dtype=float)
    subset = np.random.default_rng(seed).choice(len(atoms), min(N_STRUCTURE_ATOMS, len(atoms)),
                                                replace=False)
    centre = np.asarray(structure.assembly_refpt, dtype=float)
    epitopes = np.asarray(structure.get_target_coords(target), dtype=float)
    return dict(structure_atoms=((atoms[subset] - centre) * to_nm)[:, :2],
                structure_epitopes=((epitopes - centre) * to_nm)[:, :2])


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
    if label == REFERENCE:
        examples.update(structure_view(experiment.structure, LABELS[label][0]["probe_name"], seed))
    stacks = {modality: [] for modality in MODALITIES}
    for result in experiment.run_replicates(N_REALISATIONS, modality="All"):
        # the fluorophores of the pore, as the model places them (not a measurement)
        emitters = list(result["positions"]["SMLM"]["ch0"].values())[0]["emitters"]
        fluorophore_ring = particle_measures.ring_measures(emitters)
        if example_key(label, "SMLM", "emitters") not in examples:
            examples[example_key(label, "SMLM", "emitters")] = \
                np.asarray(emitters, dtype=float)[:, :2] - fluorophore_ring["centre"]
        for modality in MODALITIES:
            image = np.asarray(result["images"][modality]["ch0"])[0]
            pixelsize = pixelsizes[modality]
            measures = measure_npc(image, pixelsize, REGION_NM[modality])
            blob = blob_fwhm(image, pixelsize) if modality in DIFFRACTION_LIMITED else np.nan
            rows.append(dict(label=label, modality_name=modality, replicate=result["replicate"],
                             seed=seed, exp_time=exp_times.get(modality, np.nan),
                             pixelsize_nm=pixelsize, blob_fwhm_nm=blob,
                             emitter_radius_nm=fluorophore_ring["radius"]
                             if modality == "SMLM" else np.nan, **measures))
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


def figure_values(increase, scores, features):
    """The values of the figure, one row per measure (radius_nm; blob_fwhm_nm: Widefield,
    Confocal and AiryScan only), imaging method and probe: mean, s.d. and number of the
    measured images (the SMLM histograms), the increase from the direct label and the AUC of
    the probe against the direct label from that measure alone (95% bootstrap interval,
    permutation p-value), and, for SMLM, the fluorophore ring radius of the model (mean and
    s.d. over pores). For the direct label, null_auc is the mean of the null controls (two
    random halves of the direct label) of its comparisons."""
    feature_set = {"radius_nm": "radius", "blob_fwhm_nm": "blob size"}
    against = scores[scores.probe_a == REFERENCE]
    fluorophores = features.dropna(subset=["emitter_radius_nm"]) \
        .groupby("label").emitter_radius_nm.agg(["mean", "std"])
    rows = []
    for __, row in increase.iterrows():
        cell = against[(against.feature_set == feature_set[row.metric])
                       & (against.modality_name == row.modality_name)]
        values = dict(auc=np.nan, auc_low=np.nan, auc_high=np.nan, p_value=np.nan,
                      null_auc=cell.null_auc.mean())
        found = cell[cell.probe_b == row.label]
        if row.label != REFERENCE and len(found):
            values = found.iloc[0][["auc", "auc_low", "auc_high", "null_auc", "p_value"]].to_dict()
        model = dict(fluorophore_radius_mean=np.nan, fluorophore_radius_sd=np.nan)
        if row.metric == "radius_nm" and row.modality_name == "SMLM" \
                and row.label in fluorophores.index:
            model = dict(fluorophore_radius_mean=fluorophores.loc[row.label, "mean"],
                         fluorophore_radius_sd=fluorophores.loc[row.label, "std"])
        rows.append(dict(row.to_dict(), **values, **model))
    return pd.DataFrame(rows)


def smlm_auc_matrix(scores):
    """AUC of every pair of probes in SMLM from the radius alone (panel f)."""
    smlm = scores[(scores.modality_name == "SMLM") & (scores.feature_set == "radius")]
    return smlm[["probe_a", "probe_b", "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc",
                 "p_value"]].reset_index(drop=True)


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


def show_example(ax, examples, label, modality, scale_bar=False):
    """The example image of a probe and imaging method with its fitted circle, without ticks."""
    ax.set_xticks([])
    ax.set_yticks([])
    image = examples.get(example_key(label, modality, "image"))
    extent = np.array([-CROP_NM / 2, CROP_NM / 2, CROP_NM / 2, -CROP_NM / 2])
    if image is not None:
        extent = examples[example_key(label, modality, "extent")]
        ax.imshow(image, cmap="magma", extent=extent)
        centre_x, centre_y, r = examples[example_key(label, modality, "circle")]
        ax.add_patch(plt.Circle((centre_x, centre_y), r, fill=False, linestyle=":", color="w",
                                linewidth=1))
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    if scale_bar:
        x0, y0 = extent[0] + 25, extent[2] - 30
        ax.plot([x0, x0 + SCALEBAR_NM], [y0, y0], "-", color="w", lw=2)


def model_axis(ax):
    """A model panel: the pore seen from above, in nm, with the epitope circle."""
    ax.add_patch(plt.Circle((0, 0), EPITOPE_RADIUS_NM, fill=False, linestyle="--", color="0.5",
                            linewidth=0.8))
    half = MODEL_VIEW_NM / 2
    ax.set(xlim=(-half, half), ylim=(-half, half), aspect="equal", xticks=[], yticks=[])


def panel_letter(ax, letter, x=-0.12):
    ax.set_title(letter, loc="left", fontweight="bold", x=x)


def section_title(fig, spec, text):
    """A section heading above a band of the figure."""
    box = spec.get_position(fig)
    fig.text(box.x0 - 0.03, box.y1 + 0.022, text, fontsize=11, fontweight="bold", ha="left")


def published_title(label):
    published = PUBLISHED_RADIUS_NM.get(label)
    return f"{published[0]} ± {published[1]} nm" if published else "0 nm linker"


def plot_figure(examples, features, values, matrix, prefix):
    """The main figure, in four sections. 1, the model: (a) the structure from above with the
    Nup96 C-termini (epitopes), (b) the fluorophores of one pore per probe. 2, SMLM: (c) one
    image per probe with its fitted circle, (d) histograms of the fitted radii, (e) simulated
    against published radius. 3, SMLM distinguishability: (f) AUC of every pair of probes from
    the radius, (g) the radius increase from the direct label (estimated linkage error)
    against the published one. 4, the other imaging methods: (h) one image per method and
    probe, (i) the AUC of each probe against the direct label, (j) the increase of the blob
    FWHM in Widefield, Confocal and AiryScan."""
    labels = list(LABELS)
    x = np.arange(len(labels))
    radius = values[values.metric == "radius_nm"]
    blob = values[values.metric == "blob_fwhm_nm"]
    smlm = radius[radius.modality_name == "SMLM"].set_index("label").reindex(labels)
    fig = plt.figure(figsize=[12, 17])
    bands = fig.add_gridspec(4, 1, height_ratios=[0.8, 1, 1, 1.9], hspace=0.42)

    # 1: the model
    top = bands[0].subgridspec(1, 1 + len(labels), wspace=0.1)
    ax = fig.add_subplot(top[0])
    if "structure_atoms" in examples:
        ax.scatter(*examples["structure_atoms"].T, s=0.2, color="0.75", alpha=0.3, linewidths=0,
                   rasterized=True)
        ax.scatter(*examples["structure_epitopes"].T, s=6, color=PROBE_COLOURS[REFERENCE])
    model_axis(ax)
    half = MODEL_VIEW_NM / 2
    ax.plot([-half + 10, -half + 10 + MODEL_SCALEBAR_NM], [-half + 10] * 2, "-", color="k", lw=2)
    ax.set_title(f"{STRUCTURE}, Nup96 C-terminus\n(epitopes, {EPITOPE_RADIUS_NM:g} nm)",
                 fontsize=8)
    panel_letter(ax, "a")
    for k, label in enumerate(labels):
        ax = fig.add_subplot(top[1 + k])
        emitters = examples.get(example_key(label, "SMLM", "emitters"))
        if emitters is not None:
            ax.scatter(*emitters.T, s=8, color=PROBE_COLOURS[label])
        model_axis(ax)
        ax.set_title(f"{PROBE_NAMES[label]}: fluorophores\n"
                     f"{smlm.fluorophore_radius_mean[label]:.1f} ± "
                     f"{smlm.fluorophore_radius_sd[label]:.1f} nm", fontsize=8)
        if k == 0:
            panel_letter(ax, "b")

    # 2: SMLM
    middle = bands[1].subgridspec(1, 3, width_ratios=[2.3, 1.15, 1], wspace=0.3)
    images = middle[0].subgridspec(1, len(labels), wspace=0.05)
    for k, label in enumerate(labels):
        ax = fig.add_subplot(images[k])
        show_example(ax, examples, label, "SMLM", scale_bar=k == 0)
        ax.set_title(f"{PROBE_NAMES[label]}\n{published_title(label)}", fontsize=8)
        ax.spines[:].set_color(PROBE_COLOURS[label])
        ax.spines[:].set_linewidth(2)
        if k == 0:
            panel_letter(ax, "c", x=-0.2)
    ax_d = fig.add_subplot(middle[1])
    measured = features[(features.modality_name == "SMLM")].dropna(subset=["radius_nm"])
    bins = np.arange(np.floor(measured.radius_nm.min()) - 1, np.ceil(measured.radius_nm.max()) + 2)
    for label in labels:
        data = measured[measured.label == label].radius_nm
        ax_d.hist(data, bins=bins, histtype="step", linewidth=1.5, color=PROBE_COLOURS[label],
                  label=PROBE_NAMES[label])
        if label in PUBLISHED_RADIUS_NM:
            ax_d.axvline(PUBLISHED_RADIUS_NM[label][0], color=PROBE_COLOURS[label],
                         linestyle="--", linewidth=1)
    ax_d.axvline(EPITOPE_RADIUS_NM, color="0.5", linestyle=":", linewidth=1)
    ax_d.set_xlabel("Fitted radius, SMLM (nm)")
    ax_d.set_ylabel("Pores")
    ax_d.set_title("Fitted radii (dashed: published)", fontsize=9)
    ax_d.legend(frameon=False, fontsize=7)
    ax_d.spines[["top", "right"]].set_visible(False)
    panel_letter(ax_d, "d", x=-0.2)
    ax_e = fig.add_subplot(middle[2])
    for label in labels:
        row = smlm.loc[label]
        colour = PROBE_COLOURS[label]
        if label in PUBLISHED_RADIUS_NM:
            ax_e.errorbar(row.published_radius, row["mean"], xerr=row.published_sd,
                          yerr=row.sd, fmt="o", color=colour, capsize=2)
        else:  # the direct label: at the epitopes of the model
            ax_e.errorbar(EPITOPE_RADIUS_NM, row["mean"], yerr=row.sd, fmt="o", color=colour,
                          markerfacecolor="white", capsize=2)
            ax_e.annotate("direct\n(epitopes)", (EPITOPE_RADIUS_NM, row["mean"]),
                          xytext=(6, -14), textcoords="offset points", fontsize=7, color="0.2")
    lims = [EPITOPE_RADIUS_NM - 6, 70]
    ax_e.plot(lims, lims, "--", color="0.5", linewidth=1)
    ax_e.set(xlim=lims, ylim=lims, aspect="equal")
    ax_e.set_xlabel("Published radius (nm)")
    ax_e.set_ylabel("Simulated radius, SMLM (nm)")
    ax_e.set_title("Simulated vs experimental", fontsize=9)
    ax_e.spines[["top", "right"]].set_visible(False)
    panel_letter(ax_e, "e", x=-0.45)

    # 3: SMLM distinguishability
    third = bands[2].subgridspec(1, 2, width_ratios=[1, 1.4], wspace=0.45)
    ax_f = fig.add_subplot(third[0])
    # lower triangle only: rows are the second probe of a pair, columns the first
    rows_f, columns_f = labels[1:], labels[:-1]
    auc = np.full((len(rows_f), len(columns_f)), np.nan)
    for __, row in matrix.iterrows():
        i, j = rows_f.index(row.probe_b), columns_f.index(row.probe_a)
        auc[i, j] = row.auc
        colour = "w" if row.auc > 0.8 else "0.1"
        ax_f.text(j, i, f"{row.auc:.2f}\n[{row.auc_low:.2f}, {row.auc_high:.2f}]",
                  ha="center", va="center", fontsize=7, color=colour)
    shown = ax_f.imshow(np.ma.masked_invalid(auc), cmap="Blues", vmin=0.5, vmax=1)
    ax_f.set_xticks(range(len(columns_f)), [PROBE_NAMES[label] for label in columns_f],
                    fontsize=8)
    ax_f.set_yticks(range(len(rows_f)), [PROBE_NAMES[label] for label in rows_f], fontsize=8)
    ax_f.spines[:].set_visible(False)
    null = matrix[matrix.probe_a == REFERENCE].null_auc.mean()
    ax_f.set_title(f"ROC AUC from the fitted radius\n(null, direct vs direct: {null:.2f})",
                   fontsize=9)
    fig.colorbar(shown, ax=ax_f, fraction=0.046, pad=0.04, label="ROC AUC")
    panel_letter(ax_f, "f", x=-0.3)
    ax_g = fig.add_subplot(third[1])
    ax_g.plot(x, smlm.increase, "-", color="0.6", linewidth=1, zorder=1)
    for k, label in enumerate(labels):
        row = smlm.loc[label]
        ax_g.errorbar(k, row.increase, yerr=error_bars([row.increase], [row.increase_low],
                                                       [row.increase_high]),
                      fmt="o", color=PROBE_COLOURS[label], capsize=3, markersize=6, zorder=3)
        if label != REFERENCE:
            top_y = row.increase_high if np.isfinite(row.increase_high) else row.increase
            ax_g.annotate(f"AUC {row.auc:.2f}", (k, top_y), xytext=(0, 5),
                          textcoords="offset points", ha="center", fontsize=7, color="0.2")
    plot_published(ax_g, radius)
    ax_g.axhline(0, color="0.5", linestyle=":", linewidth=1)
    ax_g.set_ylabel("Radius increase from direct (nm)")
    ax_g.set_title("Estimated linkage error, SMLM", fontsize=9)
    ax_g.legend(frameon=False, fontsize=7, loc="upper left")
    probe_axis(ax_g)
    ax_g.set_xlabel("Probe (increasing linkage error)")
    panel_letter(ax_g, "g")

    # 4: the other imaging methods
    bottom = bands[3].subgridspec(1, 2, width_ratios=[1, 1.15], wspace=0.25)
    grid = bottom[0].subgridspec(len(OTHER_MODALITIES), len(labels), wspace=0.05, hspace=0.08)
    for i, modality in enumerate(OTHER_MODALITIES):
        for k, label in enumerate(labels):
            ax = fig.add_subplot(grid[i, k])
            show_example(ax, examples, label, modality,
                         scale_bar=i == len(OTHER_MODALITIES) - 1 and k == 0)
            if i == 0:
                ax.set_title(PROBE_NAMES[label], fontsize=8)
            if k == 0:
                ax.set_ylabel(modality, fontsize=8)
                ax.spines["left"].set_color(MODALITY_COLOURS[modality])
                ax.spines["left"].set_linewidth(4)
            if i == 0 and k == 0:
                panel_letter(ax, "h", x=-0.3)
    plots = bottom[1].subgridspec(2, 1, hspace=0.45)
    dodge = dict(zip(MODALITIES, np.linspace(-0.12, 0.12, len(MODALITIES))))
    ax_i = fig.add_subplot(plots[0])
    for j, modality in enumerate(MODALITIES):
        data = radius[radius.modality_name == modality].set_index("label").reindex(labels)
        colour = MODALITY_COLOURS[modality]
        reference = modality == "SMLM"
        ax_i.errorbar(x[1:] + dodge[modality], data.auc[1:],
                      yerr=error_bars(data.auc[1:], data.auc_low[1:], data.auc_high[1:]),
                      fmt="o--" if reference else "o-", color=colour, capsize=2,
                      markersize=3 if reference else 4, linewidth=1 if reference else 1.5,
                      label=f"{modality} (as f)" if reference else modality)
        ax_i.plot(x[0] + dodge[modality], data.null_auc.iloc[0], "x", color=colour,
                  label="null (direct vs direct)" if j == 0 else None)
    ax_i.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_i.set_ylim(min(0.3, np.nanmin(radius.null_auc) - 0.05), 1.05)
    ax_i.set_ylabel("ROC AUC vs direct")
    ax_i.set_title("Probe vs direct label, fitted radius", fontsize=9)
    ax_i.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    panel_letter(ax_i, "i")
    ax_j = fig.add_subplot(plots[1])
    for modality in DIFFRACTION_LIMITED:
        data = blob[blob.modality_name == modality].set_index("label").reindex(labels)
        if data.increase.isna().all():
            continue
        ax_j.errorbar(x + dodge[modality], data.increase,
                      yerr=error_bars(data.increase, data.increase_low, data.increase_high),
                      fmt="o-", color=MODALITY_COLOURS[modality], capsize=2, markersize=4,
                      label=f"{modality}: FWHM {data.reference_mean.iloc[0]:.0f} nm,\n"
                            f"  antibody AUC {data.auc.iloc[-1]:.2f}")
    plot_published(ax_j, radius)
    ax_j.axhline(0, color="0.5", linestyle=":", linewidth=1)
    ax_j.set_ylabel("FWHM increase (nm)")
    ax_j.set_title("Apparent blob size, from the direct label", fontsize=9)
    ax_j.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax_j.set_xlabel("Probe (increasing linkage error)")
    panel_letter(ax_j, "j")
    for axis in (ax_i, ax_j):
        probe_axis(axis)

    section_title(fig, bands[0], "1  Model: the nuclear pore, its epitopes and the four probes")
    section_title(fig, bands[1], "2  SMLM: images, fitted radii and the reference dataset")
    section_title(fig, bands[2], "3  SMLM: telling the probes apart (probe as a proxy for "
                                 "linkage error)")
    section_title(fig, bands[3], "4  The other imaging methods of the dataset")
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
    values = figure_values(increase, scores, features)
    values.to_csv(prefix + "npc_label_figure_values.csv", index=False)
    matrix = smlm_auc_matrix(scores)
    matrix.to_csv(prefix + "smlm_auc_matrix.csv", index=False)
    if examples:
        plot_figure(examples, features, values, matrix, prefix)
    print(values[["metric", "modality_name", "label", "n", "mean", "sd", "increase",
                  "published_increase", "auc", "null_auc"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
