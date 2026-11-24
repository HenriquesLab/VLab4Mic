"""Fig. S: can a break be detected in an experimental actin filament, and at which labelling efficiency?

The counterpart of figS_filaments_MT3D_smlm.py (microtubules) for F-actin labelled with
phalloidin, from experimental 3D SMLM localisations (ShareLoc.xyz, Zenodo 5534337).

The experimental segments are the CANDIDATE_SEGMENTS: 500 nm segments of straight, isolated,
non-overlapping actin filaments, selected once from the dataset (their centre and direction
are listed in CANDIDATE_SEGMENTS). The filament simulated is the segment REFERENCE_CANDIDATE
(numbered as on segment_overview.png), cropped from the localisations and aligned to the x
axis. Every localisation of the segment is a target, labelled directly (Linker probe, no
normals) and imaged by SMLM.

The analysis is that of fig2_filament_break_distinguishability.py with this filament: a
second model has an artificial break in the middle (every localisation within gap/2 of the
centre of the segment is removed). Each model is simulated with run_replicates at every
labelling efficiency and every image is measured with measure_filaments (widths, mask
continuity, number of segments, eccentricities, strip profiles), giving a feature table with
one row per image. A break is detected in an image when its mask splits into
MIN_SEGMENTS_BROKEN or more segments; the share of broken filaments where it is detected is
given with a 95% Wilson interval, and the same share of intact filaments is its false alarm
rate. The table is scored per break size and labelling efficiency with
distinguishability_from_features (out-of-fold ROC AUC with a 95% bootstrap interval, null
control and permutation p-value, for each feature set). A gap of 0 nm is a second,
independently simulated intact filament: a negative control.

Experimental segments, without simulation: the candidate segments (CANDIDATE_SEGMENTS) are
labelled intact or broken from half of their localisations and scored on the other half
(classify_by_frame_split). The localisations of a segment are split by frame: those of odd
frames are rendered as the simulated images (2 nm pixels, Gaussian of the 5 nm lateral
precision) and measured with measure_filaments, and the segment is broken when its mask splits
or its strip profiles show a break, intact when neither. The localisations of even frames are
rendered and measured the same way, and these features are the ones scored: a real break is in
both halves, while the noise that set the label is not carried over, so the scoring is not
circular. The share of segments where a break is detected is given per class, and the
classes are scored with distinguishability_from_features.

Experimental labelling efficiency, without simulation: the filament figure again, with the
labelling efficiency taken from the data. From the same candidates, a reference set of
N_SEGMENTS_PER_SET segments with no detected break, and sets of up to N_SEGMENTS_PER_SET
sparser segments (within DENSITY_TOLERANCE of DENSITY_TARGETS times its number of
localisations); the labelling efficiency of a set is its mean number of localisations over
that of the reference set. Every segment is measured intact and with a gap of DENSITY_GAPS_NM
cut in its middle (its localisations within gap/2 of the centre removed), and the break
detection and distinguishability are those of the simulated filaments. Its outputs are in the
subfolder experimental_labelling_efficiency, with a README.txt that describes every file.

Outputs, in the folder <date>figS_filaments_actin3D of the VLab4Mic output directory (created
if it does not exist): segment_overview.png (the dataset with the numbered candidate segments
and the reference segment), filament.png (the targets of the reference segment),
filament_figure.png, the main figure:
(a) one filament per break of FIGURE_EXAMPLE_GAPS_NM imaged at every labelling efficiency,
(b) the percentage of broken filaments whose break is detected and the false alarm rate of the
intact filaments, and (c) the distinguishability (AUC) of intact and broken filaments from the
number of segments alone; features.csv, distinguishability.csv, break_detection.csv,
filament_figure_values.csv, settings.json, strip_profiles.png, distinguishability_heatmap.png,
feature_distributions.png, auc_vs_labelling.png, examples_gap<break>.png, models.png,
models_grid.png and examples.npz, as in fig2_filament_break_distinguishability.py; for the
experimental segments candidates.csv, experimental_segments_figure.png ((a) the segments on
the dataset, (b) their images (even frames), (c) the percentage of segments with a detected
break per class and (d) the distinguishability of the classes), experimental_features.csv,
experimental_break_detection.csv, experimental_distinguishability.csv and
experimental_segments.npz; and the subfolder experimental_labelling_efficiency
(filament_figure.png, sets.png, distinguishability_heatmap.png, the tables and README.txt).

Run:   python figS_filaments_actin3D_smlm.py
Score an existing feature table again, without simulating:
       python figS_filaments_actin3D_smlm.py --features <folder>/features.csv
(the outputs are written next to the feature table)

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and condition; e.g. 5 for a
quick run), VLAB4MIC_N_PERMUTATIONS (200), VLAB4MIC_N_WORKERS (1; simulations and scoring in
parallel processes), VLAB4MIC_SAVE_IMAGES (0; 1 saves the images of every model and condition
as a TIFF stack in <folder>/images), VLAB4MIC_SMLM_PATH (the .smlm file of the localisations,
local or URL; default: Zenodo; CANDIDATE_SEGMENTS are those of the default dataset),
VLAB4MIC_REFERENCE_CANDIDATE (the reference segment, instead of REFERENCE_CANDIDATE),
VLAB4MIC_N_SEGMENTS_PER_SET (20; experimental segments per set of the experimental labelling
efficiency, as many as the data allows),
VLAB4MIC_N_EXPERIMENTAL_PER_CLASS (12; experimental segments per class, intact and broken).
"""

import copy
import io
import json
import os
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import tifffile
from matplotlib.patches import Polygon
from scipy import ndimage
from scipy.spatial import cKDTree
from tqdm.auto import tqdm

from vlab4mic import experiments
from vlab4mic.analysis import particle_measures
from vlab4mic.analysis.distinguishability import (
    distinguishability_from_features,
    plot_distinguishability_heatmaps,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24
OUTPUT_NAME = "figS_filaments_actin3D"  # base name of the output folder (after the date)

# Experimental 3D SMLM localisations of actin labelled with phalloidin (ShareLoc.xyz, Zenodo 5534337)
SMLM_EXPERIMENTAL_PATH = os.environ.get(
    "VLAB4MIC_SMLM_PATH",
    "https://zenodo.org/records/5534337/files/160623_C1C10_act3D_3342K/data.smlm")
# straight, isolated, non-overlapping actin filament segments of the dataset, the most selective
CANDIDATE_SEGMENTS = [
    (36815.24, 21558.60, 51), (37139.95, 27802.58, 90), (36093.80, 24014.69, 153),
    (36856.64, 30455.83, 3), (38871.12, 19332.44, 3), (35842.89, 21547.91, 96),
    (37066.70, 23037.86, 111), (15609.19, 20264.76, 60), (34588.68, 28558.73, 45),
    (12789.75, 31264.25, 99), (36295.59, 27323.27, 57), (35850.83, 20801.16, 69),
    (38894.56, 31262.66, 45), (38116.17, 30292.01, 36), (35148.38, 28010.38, 15),
    (31893.90, 32796.58, 171), (35614.84, 22562.64, 9), (39066.80, 22791.64, 144),
    (38303.68, 25792.93, 51), (33627.55, 24513.60, 33), (38860.33, 30244.91, 9),
    (33865.78, 25251.90, 177), (27859.52, 38331.69, 6), (33890.74, 32271.90, 63),
    (34870.05, 24000.59, 39), (35865.95, 31091.26, 3), (33361.95, 14842.98, 27),
    (39344.66, 24286.43, 117), (37877.52, 28258.29, 54), (37094.93, 31802.53, 150),
    (33117.60, 24996.12, 153), (33876.44, 13314.79, 162), (31612.24, 36340.64, 168),
    (36601.56, 28771.79, 18), (27104.23, 36530.49, 15), (24361.26, 23046.88, 168),
    (38578.62, 20260.90, 66), (31850.44, 26273.70, 33), (33892.31, 38301.32, 171),
    (37656.25, 24330.46, 135), (33334.60, 34296.92, 3), (35582.75, 32573.17, 63),
    (31107.07, 4305.80, 147), (20125.70, 23574.51, 3), (38331.83, 31050.47, 66),
    (39100.41, 32573.29, 15), (38133.74, 32045.93, 33), (34330.40, 26793.27, 150),
    (33359.26, 5305.12, 111), (37679.04, 20825.06, 84), (31840.46, 12812.17, 66),
    (36106.09, 29733.95, 159), (34167.81, 25807.29, 114), (36893.98, 32492.31, 42),
    (16108.28, 3288.03, 114), (37337.83, 26740.72, 114), (39038.35, 26015.40, 102),
    (39049.49, 28834.84, 54), (12407.63, 35570.27, 111), (17378.34, 21823.39, 156),
    (35659.84, 34123.13, 168), (36329.63, 25285.85, 36), (33103.63, 31298.09, 84),
    (36896.05, 24024.09, 42), (36089.21, 28356.80, 18), (36385.57, 33361.40, 156),
    (33339.58, 33600.48, 48), (33328.95, 13835.33, 15), (38642.32, 27526.86, 99),
    (30367.04, 34733.13, 177), (33865.99, 31067.47, 9), (15871.55, 27784.94, 57),
    (35664.97, 29262.76, 45), (30851.35, 27032.15, 117), (35336.30, 27036.67, 75),
    (39615.99, 30591.13, 42), (37326.62, 25838.68, 30), (31859.06, 15289.62, 0),
    (39659.93, 27061.14, 108), (13102.47, 36278.76, 69), (27863.03, 16542.29, 144),
    (37612.07, 23546.73, 174), (34096.04, 34527.29, 33), (31103.75, 15816.25, 156),
    (34862.60, 25283.09, 15), (37343.20, 33038.48, 57), (33368.29, 30246.40, 168),
    (14864.77, 27526.19, 99), (34618.66, 32291.15, 24), (30115.27, 10779.71, 156),
    (36853.24, 26326.45, 165), (19097.85, 21807.30, 66), (38326.66, 29046.12, 105),
    (38644.63, 34089.14, 159),
]
# the segment simulated (index in CANDIDATE_SEGMENTS, as numbered on segment_overview.png)
REFERENCE_CANDIDATE = int(os.environ.get("VLAB4MIC_REFERENCE_CANDIDATE", 42))
SEGMENT_LENGTH_NM = 500  # as the simulated filament in fig2_filaments_actin.py
SEGMENT_WIDTH_NM = 100  # narrower than for microtubules: the actin network is dense

STRUCTURE = "6T25"  # phalloidin-stabilised F-actin, as fig2_filaments_actin.py
TARGET = "Phalloidin_EEP"
TARGET_VALUE = dict(atoms=["O2"], residues=["EEP"])  # 5-hydroxyl of 4,5-dihydroxyleucine
SAMPLE_DIMENSIONS = [1000, 1000, 10]  # nm

LABELLING_EFFICIENCIES = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1]
GAPS_NM = [0, 5, 10, 25, 50, 100, 200]  # length of the break; 0: intact (negative control)
# breaks of the main figure, one line each in b and c; 0: the intact filament
FIGURE_EXAMPLE_GAPS_NM = [0, 50, 100]

# a break is detected when the mask splits into this many segments or more
MIN_SEGMENTS_BROKEN = 2

# example images: broken filaments whose break is detected first, then others
N_EXAMPLES = 8  # per break size and labelling efficiency
# crop of the example images around the image centre, shown with the filament horizontal
EXAMPLE_LENGTH_NM = 620  # along the filament
EXAMPLE_WIDTH_NM = 120  # across it

RENDERING_SIGMA_NM = 5  # as the simulated images: the lateral precision

# experimental segments: the candidates labelled from their odd frames and scored on
# their even frames (classify_by_frame_split), up to this many per class
N_EXPERIMENTAL_PER_CLASS = int(os.environ.get("VLAB4MIC_N_EXPERIMENTAL_PER_CLASS", 12))
MIN_HALF_LOCALISATIONS = 30  # in each half; segments with fewer are left out

# experimental labelling efficiency: sets of candidate segments. Set 0, the reference:
# the N_SEGMENTS_PER_SET densest segments where neither the mask nor the strip profiles detect
# a break. The other sets: up to N_SEGMENTS_PER_SET segments closest to DENSITY_TARGETS times
# the mean number of localisations of the reference, within DENSITY_TOLERANCE of it (chosen by
# density only); a set with fewer than DENSITY_MIN_SEGMENTS is left out. Each set's labelling
# efficiency is its mean number of localisations over that of the reference.
N_SEGMENTS_PER_SET = int(os.environ.get("VLAB4MIC_N_SEGMENTS_PER_SET", 20))
DENSITY_TARGETS = [1, 0.8, 0.6, 0.45, 0.33, 0.25]
DENSITY_TOLERANCE = 0.15  # relative to the target number of localisations
DENSITY_MIN_SEGMENTS = 5  # per set, for 5-fold scoring
# breaks cut into every segment; no 0 nm: without a second realisation it would be the intact
# segment itself (the null control of distinguishability_from_features is the negative control)
DENSITY_GAPS_NM = [5, 10, 25, 50, 100, 200]
DENSITY_FOLDER = "experimental_labelling_efficiency"  # subfolder of the outputs

# measure_filaments columns used as features (the others are settings)
MASK_FEATURES = ["mask_continuity", "mask_major_lenght", "mask_number_of_segments",
                 "min_eccentricity", "max_eccentricity", "mean_eccentricity"]
STRIP_FEATURES = ["strip_fraction_flat", "strip_n_breaks", "strip_longest_break",
                  "strip_min_amplitude"]
WIDTHS = ["mean_filament_width", "std_filament_width", "median_filament_width", "s_width"]
FEATURES = WIDTHS + MASK_FEATURES + STRIP_FEATURES
FEATURE_SETS = {"segments": ["mask_number_of_segments"], "mask": MASK_FEATURES,
                "strip profiles": STRIP_FEATURES, "all": FEATURES}
CONTROLS = {feature: feature for feature in FEATURES}  # each alone
# images without a measurable filament (e.g. too few labels)
FILL_VALUES = dict({width: 0 for width in WIDTHS}, strip_fraction_flat=1, strip_min_amplitude=0)
CONDITIONS = ["gap_nm", "labelling_efficiency"]

COMMON = dict(
    structure=STRUCTURE,
    sample_dimensions=SAMPLE_DIMENSIONS,
    number_of_particles=1,
    modality="SMLM",
    SMLM=dict(
        lateral_precision=5,  # nm
        axial_precision=5,  # nm
        n_localisations=30,
        precision_model="fixed",  # not photon-limited
    ),
    clear_experiment=True,
    random_seed=RANDOM_SEED,
)
FILAMENTS_PARAMS = dict(
    pixel_size=2,
    max_length_nm=SEGMENT_LENGTH_NM,
    modality="SMLM",
    channel="ch0",
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=1,  # px (2 nm), Gaussian before the triangle threshold
    threshold_method="triangle",
    strip_kwargs=dict(
        strip_length=10,  # nm, along the filament axis
        half_width=50,  # nm, profiles across the axis
        background_width=10,  # nm at each end of a profile
        flat_fraction=0.1,  # flat: amplitude below this fraction of the reference
    ),
)


def load_smlm_table(path):
    """Read a .smlm file (zip with manifest.json and a binary table) into a DataFrame."""
    if path.startswith("http"):
        response = requests.get(path)
        response.raise_for_status()
        archive = zipfile.ZipFile(io.BytesIO(response.content))
    else:
        archive = zipfile.ZipFile(path)
    manifest = json.loads(archive.read("manifest.json"))
    table_file = manifest["files"][0]
    table_format = manifest["formats"][table_file["format"]]
    dtype = np.dtype([(h, "<" + np.dtype(t).str[1:]) for h, t in zip(table_format["headers"], table_format["dtype"])])
    table = np.frombuffer(archive.read(table_file["name"]), dtype=dtype, count=table_file["rows"])
    return pd.DataFrame(table)


_localisation_cache = {}  # xy, z and the cKDTree of xy of a table of localisations


def localisation_arrays(localisations):
    """xy and z of the localisations (nm) and the cKDTree of xy, built once per table."""
    cached = _localisation_cache.get(id(localisations))
    if cached is None or cached[0] is not localisations:
        xy = localisations[["x", "y"]].to_numpy(dtype=float)
        cached = (localisations, xy, localisations["z"].to_numpy(dtype=float), cKDTree(xy))
        _localisation_cache[id(localisations)] = cached
    return cached[1:]


def crop_filament_segment(localisations, center, length, width, *, direction):
    """
    Crop a straight filament segment and align it to the x axis.

    center (xy, nm) and direction (unit xy) are those of the segment. Localisations within
    length/2 along and width/2 across its axis are kept, rotated so the filament is horizontal
    and centred at the origin (z centred on its median), in the order of the table. Returns the
    segment, the mask of the localisations kept, its centre and its direction.
    """
    center, direction = np.asarray(center, dtype=float), np.asarray(direction, dtype=float)
    xy, z, tree = localisation_arrays(localisations)
    index = np.sort(np.array(tree.query_ball_point(center, np.hypot(length / 2, width / 2)), dtype=int))
    rotation = np.array([[direction[0], direction[1]], [-direction[1], direction[0]]])
    along_across = (xy[index] - center) @ rotation.T
    inside = (np.abs(along_across[:, 0]) <= length / 2) & (np.abs(along_across[:, 1]) <= width / 2)
    keep = np.zeros(len(xy), dtype=bool)
    keep[index[inside]] = True
    segment_z = z[index[inside]]
    segment = np.column_stack([along_across[inside],
                               segment_z - (np.median(segment_z) if len(segment_z) else 0)])
    return segment, keep, center, direction


def experimental_template(fluorophore_positions):
    """Targets of the experimental filament, as particle source data: without normals the
    targets are labelled directly, each localisation becomes a fluorophore at its own
    position."""
    return dict(
        targets={TARGET: dict(coordinates=fluorophore_positions.copy(), normals=None)},
        scale=1e-9,
        reference_point=np.array([0.0, 0.0, 0.0]),
        axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0])),
    )


def linker_probe(labelling_efficiency=1):
    """Every target labelled directly: a linker of length 0 places the fluorophore on it."""
    return dict(
        probe_template="Linker",
        probe_name=TARGET,
        probe_distance_to_epitope=0,
        probe_target_type="Atom_residue",
        probe_target_value=copy.deepcopy(TARGET_VALUE),
        probe_steric_hindrance=0,  # localisations are not excluded by probe size
        labelling_efficiency=labelling_efficiency,
    )


### Break distinguishability (as fig2_filament_break_distinguishability.py)

def position_along_filament(coordinates):
    """Position of each point along the principal axis, from the median (nm)."""
    centred = coordinates - coordinates.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov(centred.T))
    along = centred @ eigenvectors[:, -1]
    return along - np.median(along)


def break_filament(template, gap_nm):
    """Copy of the template without the targets within gap_nm / 2 of its centre."""
    broken = copy.deepcopy(template)
    if gap_nm <= 0:
        return broken
    target = broken["targets"][TARGET]
    along = position_along_filament(target["coordinates"])
    span = np.ptp(along)
    assert 0.9 * SEGMENT_LENGTH_NM < span < 1.2 * SEGMENT_LENGTH_NM, \
        f"filament spans {span:.0f} along its axis: coordinates not in nm?"
    keep = np.abs(along) >= gap_nm / 2
    target["coordinates"] = target["coordinates"][keep]
    if target["normals"] is not None:
        target["normals"] = target["normals"][keep]
    return broken


def measure_filaments(replicates,
                      max_length_nm=500,
                      pixel_size=1,
                      modality="SMLM",
                      channel="ch0",
                      fill_holes=True,
                      smooth_sigma=1,
                      row_binning=1,
                      min_snr=2,
                      mode="single",
                      threshold_method="otsu",
                      threshold_method_kwargs=None,
                      strip_kwargs=None,):
    """Features of the first frame of every realisation, one row each: filament widths, mask
    continuity, number of mask segments and their eccentricities, and the strip profile
    measures. Also returns the masks and strip profiles of every image."""
    masks = []
    rows = []
    strips_list = []
    for result in replicates:
        image = result["images"][modality][channel][0]
        mask, table, __ = particle_measures.threshold_mask_and_measure(image, method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = particle_measures.measure_filament_width(image, pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = particle_measures.simple_width(image)
        strips = particle_measures.filament_strip_profiles(image, pixel_size=pixel_size, **(strip_kwargs or {}))
        strips_list.append(strips)
        masks.append(mask)
        n_segments = table.shape[0]
        if n_segments == 0:
            min_eccentricity = max_eccentricity = mean_eccentricity = 0
            major_lenght = 0
            apparent_continuity = 0
        else:
            min_eccentricity = np.min(table["eccentricity"])
            max_eccentricity = np.max(table["eccentricity"])
            mean_eccentricity = np.mean(table["eccentricity"])
            major_lenght = table.iloc[table['major_axis_length'].idxmax()]["major_axis_length"]
            apparent_continuity = (major_lenght / max_length_nm) / n_segments

        rows.append(dict(
            replicate=result["replicate"],
            mean_filament_width=filament_metrics["mean"],
            std_filament_width=filament_metrics["std"],
            median_filament_width=filament_metrics["median"],
            mask_continuity=apparent_continuity,
            mask_major_lenght=major_lenght,
            mask_number_of_segments=n_segments,
            mask_prefilter_smooth_sigma=smooth_sigma,
            pixelsize=pixel_size,
            row_binning=row_binning,
            s_width=s_width,
            min_eccentricity=min_eccentricity,
            max_eccentricity=max_eccentricity,
            mean_eccentricity=mean_eccentricity,
            strip_fraction_flat=strips["fraction_flat"],
            strip_n_breaks=strips["n_breaks"],
            strip_longest_break=strips["longest_break"],
            strip_min_amplitude=strips["min_amplitude"],
        ))
    df = pd.DataFrame(rows)
    return df, masks, strips_list


def build_experiment(labelling_efficiency, seed, particle_data, run_simulation=True,
                     modules=("coordinate_field", "imager")):
    """An experiment whose particle is particle_data, its targets labelled directly (Linker)
    with labelling_efficiency. image_vsample first images a sample of STRUCTURE, which is
    discarded; run_simulation=True keeps that simulation so the realisations that follow draw
    the same random numbers as before."""
    __, __, experiment = experiments.image_vsample(
        probe_list=[linker_probe(labelling_efficiency)],
        run_simulation=run_simulation,
        **dict(COMMON, random_seed=seed),
    )
    experiment.particle.source = dict()
    experiment.particle.load_source(**copy.deepcopy(particle_data))
    experiment.set_virtualsample_params(
        sample_dimensions=SAMPLE_DIMENSIONS,
        update_mode=False,
    )
    experiment.particle.generate_instance()
    experiment.build(modules=list(modules))
    return experiment


def simulate_replicates(labelling_efficiency, seed, particle_data, n_replicates=N_REALISATIONS):
    """Independent realisations of a filament (run_replicates output): its targets labelled
    directly (Linker) with labelling_efficiency."""
    experiment = build_experiment(labelling_efficiency, seed, particle_data)
    return experiment.run_replicates(n_replicates=n_replicates, modality="SMLM")


def break_detected(table):
    """Rows where the mask splits into MIN_SEGMENTS_BROKEN or more segments."""
    return table.mask_number_of_segments >= MIN_SEGMENTS_BROKEN


def crop_filament(image):
    """The image centre, EXAMPLE_LENGTH_NM along the filament by EXAMPLE_WIDTH_NM across,
    with the filament horizontal: it runs along the rows of the image, so the crop is
    transposed."""
    image = np.asarray(image, dtype=float)
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    rows, cols = np.array(image.shape) // 2
    half_length = int(round(EXAMPLE_LENGTH_NM / pixel_size / 2))
    half_width = int(round(EXAMPLE_WIDTH_NM / pixel_size / 2))
    crop = image[rows - half_length:rows + half_length, cols - half_width:cols + half_width]
    return crop.T.astype(np.float32)


EXTENT = [-EXAMPLE_LENGTH_NM / 2, EXAMPLE_LENGTH_NM / 2, -EXAMPLE_WIDTH_NM / 2,
          EXAMPLE_WIDTH_NM / 2]


def show(ax, image):
    """A cropped filament image, without ticks."""
    ax.imshow(image, cmap="magma", extent=EXTENT, origin="lower")
    ax.set_xticks([])
    ax.set_yticks([])


def scale_bar(ax, colour="w"):
    """100 nm, at the lower left of a cropped filament image."""
    x0, y0 = -EXAMPLE_LENGTH_NM / 2 + 15, -EXAMPLE_WIDTH_NM / 2 + 15
    ax.plot([x0, x0 + 100], [y0, y0], "-", color=colour, lw=2)


def model_key(model, gap, le):
    return f"model_{model}_gap_{gap:g}_labelling_efficiency_{le:g}"


def epitope_key(model, gap):
    return f"model_{model}_gap_{gap:g}_epitopes"


def example_key(gap, le):
    return f"gap_{gap:g}_labelling_efficiency_{le:g}"


def model_label(model, gap):
    return "intact" if model == "intact" else f"break {gap:g} nm"


def wilson_interval(k, n, z=1.96):
    """95% interval of a proportion k / n."""
    p = k / n
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return max(centre - half, 0.0), min(centre + half, 1.0)


def error_bars(value, low, high):
    """Lengths of the error bars below and above value, never negative (rounding can put
    value a hair outside [low, high], which matplotlib refuses)."""
    value, low, high = (np.asarray(x, dtype=float) for x in (value, low, high))
    return [np.clip(value - low, 0, None), np.clip(high - value, 0, None)]


def labelling_axis(ax, les=LABELLING_EFFICIENCIES):
    """Labelling efficiency on a log scale, ticks at the values of les."""
    les = sorted(les)
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks(les, [f"{le:.2g}" for le in les], fontsize=8)
    ax.set_xlabel("Labelling efficiency")
    ax.spines[["top", "right"]].set_visible(False)


def simulate(job):
    """One model at one labelling efficiency: its feature table (one row per image), the
    first image (cropped) and, for the broken filament, N_EXAMPLES images (cropped), those
    whose break is detected first."""
    model, gap_nm, labelling_efficiency, seed, particle_data, images_dir = job
    replicates = simulate_replicates(labelling_efficiency, seed, particle_data)
    table = measure_filaments(replicates=replicates, **FILAMENTS_PARAMS)[0]
    table = table.assign(model=model, gap_nm=gap_nm, labelling_efficiency=labelling_efficiency,
                         seed=seed)
    # the frame measured by measure_filaments, one per realisation
    images = [np.asarray(result["images"]["SMLM"]["ch0"][0]) for result in replicates]
    if images_dir is not None:
        os.makedirs(images_dir, exist_ok=True)
        tifffile.imwrite(
            os.path.join(images_dir, f"{model}_gap_{gap_nm}_labelling_efficiency_"
                                     f"{labelling_efficiency}.tif"),
            np.array(images, dtype=np.float32))
    examples = {model_key(model, gap_nm, labelling_efficiency): crop_filament(images[0])}
    if model == "broken":
        replicate = table["replicate"].astype(int)
        detected = break_detected(table)
        replicas = pd.concat([replicate[detected], replicate[~detected]]).to_numpy()[:N_EXAMPLES]
        key = example_key(gap_nm, labelling_efficiency)
        examples[key] = np.array([crop_filament(images[r]) for r in replicas])
        examples[key + "_replicas"] = replicas
    return table, examples


def model_epitopes(models):
    """Targets of every model in the plane, from the centre of the intact filament (nm)."""
    centre = models[0][2]["targets"][TARGET]["coordinates"].mean(axis=0)
    return {epitope_key(model, gap):
            data["targets"][TARGET]["coordinates"][:, :2] - centre[:2]
            for model, gap, data in models}


def models_shown():
    """(model, break) of the intact filament and of every broken one."""
    return [("intact", 0)] + [("broken", gap) for gap in GAPS_NM]


def plot_models(examples, prefix):
    """The models, one per column: their targets from above (removed ones open) and an
    image at every labelling efficiency (rows)."""
    models = models_shown()
    les = sorted(LABELLING_EFFICIENCIES)
    all_epitopes = examples[epitope_key("intact", 0)]
    fig, axes = plt.subplots(1 + len(les), len(models), squeeze=False,
                             figsize=[2.4 * len(models), 0.75 * (1 + len(les)) + 0.6])
    for column, (model, gap) in zip(axes.T, models):
        epitopes = examples[epitope_key(model, gap)]
        column[0].scatter(*all_epitopes.T, s=0.5, facecolors="none", edgecolors="0.7",
                          linewidths=0.3)
        column[0].scatter(*epitopes.T, s=0.5, color="k")
        column[0].set(xlim=EXTENT[:2], ylim=EXTENT[2:], aspect="equal", xticks=[], yticks=[])
        column[0].set_title(f"{model_label(model, gap)}\n{len(epitopes)} localisations", fontsize=8)
        for ax, le in zip(column[1:], les):
            show(ax, examples[model_key(model, gap, le)])
    scale_bar(axes[0, 0], colour="k")
    axes[0, 0].set_ylabel("Targets", fontsize=8)
    for ax, le in zip(axes[1:, 0], les):
        scale_bar(ax)
        ax.set_ylabel(f"LE {le:g}", fontsize=8)
    fig.suptitle("Experimental actin filament models: targets (experimental localisations; black, "
                 "kept; open, removed) and one image at every labelling efficiency; scale bar "
                 "100 nm", x=0.02, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(prefix + "models.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strip_examples(template, prefix, labelling_efficiency=1, gaps=FIGURE_EXAMPLE_GAPS_NM):
    """Mask and strip profiles of one image of the intact and broken filaments."""
    flat_fraction = FILAMENTS_PARAMS["strip_kwargs"]["flat_fraction"]
    fig, axes = plt.subplots(len(gaps), 4, figsize=(16, 3.2 * len(gaps)),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.4, 2]))
    for row, gap in zip(axes, gaps):
        replicates = simulate_replicates(labelling_efficiency, RANDOM_SEED + gap,
                                         break_filament(template, gap), n_replicates=1)
        table, masks, strips = measure_filaments(replicates, **FILAMENTS_PARAMS)
        image, mask, strip = replicates[0]["images"]["SMLM"]["ch0"][0], masks[0], strips[0]
        plot_strip_row(row, image, mask, strip, table, f"gap {gap} nm, LE {labelling_efficiency}",
                       flat_fraction)
    fig.tight_layout()
    fig.savefig(prefix + "strip_profiles.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strip_row(row, image, mask, strip, table, title, flat_fraction):
    """Image and mask around the filament, with the axis and the flat strips, the line
    profile of every strip and their amplitude along the filament, on four axes."""
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    rows, cols = np.nonzero(image > 0.05 * image.max())
    crop = (slice(max(rows.min() - 30, 0), rows.max() + 30),
            slice(max(cols.min() - 30, 0), cols.max() + 30))
    origin = np.array([crop[0].start, crop[1].start])
    row[0].imshow(image[crop], cmap="grey", vmax=np.percentile(image[image > 0], 99.5))
    axis_points = (strip["centroid"] - origin
                   + np.outer(strip["positions"] / pixel_size, strip["direction"]))
    row[0].plot(axis_points[:, 1], axis_points[:, 0], c="C1", lw=0.8)
    row[0].scatter(axis_points[strip["flat"], 1], axis_points[strip["flat"], 0], c="red", s=6)
    row[0].set_title(title, loc="left", fontsize=9)
    row[1].imshow(mask[crop], cmap="grey")
    row[1].set_title(f"mask: {int(table.mask_number_of_segments[0])} segment(s)",
                     loc="left", fontsize=9)
    for ax in row[:2]:
        ax.axis("off")
    # one line profile per strip
    extent = [strip["offsets"][0], strip["offsets"][-1],
              strip["positions"][-1], strip["positions"][0]]
    row[2].imshow(strip["profiles"], cmap="grey", aspect="auto", extent=extent)
    row[2].set_xlabel("across the axis (nm)")
    row[2].set_ylabel("along the axis (nm)")
    row[2].set_title("line profile of each strip", loc="left", fontsize=9)
    # amplitude of the profiles along the filament
    reference = strip["amplitudes"].max() if strip["amplitudes"].max() > 0 else 1
    row[3].plot(strip["positions"], strip["amplitudes"] / reference, "o-", c="k", ms=2, lw=0.8)
    row[3].scatter(strip["positions"][strip["flat"]],
                   strip["amplitudes"][strip["flat"]] / reference, c="red", s=12, zorder=3,
                   label="flat")
    row[3].axhline(flat_fraction * np.percentile(strip["amplitudes"], 90) / reference,
                   ls="--", c="grey", lw=0.8, label="flat threshold")
    row[3].set_xlabel("along the axis (nm)")
    row[3].set_ylabel("profile amplitude (norm.)")
    row[3].set_title(f"{strip['n_breaks']} break(s), longest {strip['longest_break']:.0f} nm",
                     loc="left", fontsize=9)
    row[3].legend(frameon=False, fontsize=8)


def feature_table(tables, gaps=GAPS_NM):
    """Intact and broken rows for every break size of gaps."""
    features = pd.concat(tables, ignore_index=True)
    intact = features[features.model == "intact"]
    broken = features[features.model == "broken"]
    # the intact images are shared by every break size
    return pd.concat([broken] + [intact.assign(gap_nm=gap) for gap in gaps],
                     ignore_index=True)


def break_detection(features):
    """Per model, break size and labelling efficiency: the share of filaments whose mask splits
    into MIN_SEGMENTS_BROKEN or more segments (broken: detected breaks; intact: false alarms),
    and the median, mean and s.d. of their number of segments."""
    intact = features[features.model == "intact"]
    # the intact filaments appear once per break size: keep one copy
    intact = intact[intact.gap_nm == intact.gap_nm.min()].assign(gap_nm=0)
    rows = []
    for (model, gap, le), data in pd.concat([intact, features[features.model == "broken"]]) \
            .groupby(["model", "gap_nm", "labelling_efficiency"]):
        split = int(break_detected(data).sum())
        low, high = wilson_interval(split, len(data))
        rows.append(dict(model=model, gap_nm=gap, labelling_efficiency=le,
                         n_filaments=len(data), n_split=split, fraction_split=split / len(data),
                         fraction_low=low, fraction_high=high,
                         median_segments=data.mask_number_of_segments.median(),
                         mean_segments=data.mask_number_of_segments.mean(),
                         std_segments=data.mask_number_of_segments.std()))
    return pd.DataFrame(rows).sort_values(["model", "gap_nm", "labelling_efficiency"],
                                          ascending=[False, True, True], ignore_index=True)


def score_feature_sets(features, by=None, n_workers=1, progress=True, label="features"):
    """Distinguishability of intact and broken filaments (distinguishability_from_features)
    for every feature set of FEATURE_SETS, each feature alone as controls of "all", in one
    table."""
    scores = []
    for name, columns in FEATURE_SETS.items():
        if progress:
            print(f"{label}: {name}")
        table = distinguishability_from_features(
            features, "model", ("intact", "broken"), columns, by=by,
            controls=CONTROLS if name == "all" else None,
            replicate_column="replicate", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=n_workers, progress=progress,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    return pd.concat(scores, ignore_index=True)


def figure_values(detection, scores):
    """The values of b and c of the main figure, one row per model (the intact filament and every
    break of FIGURE_EXAMPLE_GAPS_NM) and labelling efficiency: b, the percentage of filaments whose
    mask splits into MIN_SEGMENTS_BROKEN or more segments (95% Wilson interval), with the mean and
    s.d. of their number of segments; c, the AUC of intact against broken filaments from the
    number of segments alone (95% bootstrap interval), its null control (intact against intact)
    and permutation p-value."""
    breaks = [gap for gap in FIGURE_EXAMPLE_GAPS_NM if gap > 0]
    shown = detection[(detection.model == "intact")
                      | ((detection.model == "broken") & detection.gap_nm.isin(breaks))]
    values = pd.DataFrame(dict(
        model=shown.model, gap_nm=shown.gap_nm, labelling_efficiency=shown.labelling_efficiency,
        n_filaments=shown.n_filaments, b_percent_detected=100 * shown.fraction_split,
        b_percent_low=100 * shown.fraction_low, b_percent_high=100 * shown.fraction_high,
        b_mean_segments=shown.mean_segments, b_std_segments=shown.std_segments))
    segments = scores.loc[scores.feature_set == "segments",
                          CONDITIONS + ["auc", "auc_low", "auc_high", "null_auc", "p_value"]]
    segments = segments.rename(columns=lambda column: column if column in CONDITIONS
                               else "c_" + column)
    # c compares the intact filament with each break: no AUC on the intact rows
    return values.merge(segments.assign(model="broken"), how="left", on=["model"] + CONDITIONS)


def plot_feature_distributions(features, prefix):
    """Mean +- s.d. of the continuity features against labelling efficiency, one line per break
    size (a break of 0 nm: the second intact filament)."""
    broken = features[features.model == "broken"]
    gaps = sorted(broken.gap_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(gaps)))
    columns = ["mask_number_of_segments", "mask_continuity", "strip_n_breaks",
               "strip_longest_break", "strip_fraction_flat", "mean_filament_width"]
    fig, axes = plt.subplots(1, len(columns), figsize=[3.5 * len(columns), 3])
    for ax, feature in zip(axes, columns):
        for gap, colour in zip(gaps, colours):
            summary = broken[np.isclose(broken.gap_nm, gap)] \
                .groupby("labelling_efficiency")[feature].agg(["mean", "std"])
            ax.errorbar(summary.index, summary["mean"], yerr=summary["std"], fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{gap:g}")
        labelling_axis(ax)
        ax.set_ylabel(feature.replace("_", " "))
    axes[-1].legend(title="Break (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "feature_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_auc_vs_labelling(scores, prefix):
    """AUC (95% interval) against labelling efficiency, one line per break size."""
    gaps = sorted(scores.gap_nm.unique())
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(gaps)))
    fig, axes = plt.subplots(1, len(FEATURE_SETS), figsize=[3.5 * len(FEATURE_SETS), 3],
                             sharey=True)
    for ax, name in zip(axes, FEATURE_SETS):
        data = scores[scores.feature_set == name]
        for gap, colour in zip(gaps, colours):
            sub = data[np.isclose(data.gap_nm, gap)].sort_values("labelling_efficiency")
            ax.errorbar(sub.labelling_efficiency, sub.auc,
                        yerr=error_bars(sub.auc, sub.auc_low, sub.auc_high), fmt="o-",
                        color=colour, capsize=2, markersize=3, label=f"{gap:g}")
        ax.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
        ax.set_ylim(0.3, 1.05)
        labelling_axis(ax)
        ax.set_title(name, loc="left")
    axes[0].set_ylabel("ROC AUC (intact vs broken)")
    axes[-1].legend(title="Break (nm)", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "auc_vs_labelling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_examples(features, examples, prefix):
    """Per break size: example broken filaments, one column per labelling efficiency, those
    whose break is detected first (green frame), each with its segments and longest break."""
    broken = features[features.model == "broken"]
    les = sorted(LABELLING_EFFICIENCIES)
    for gap in GAPS_NM:
        # fewer rows when there are fewer realisations than N_EXAMPLES
        n_rows = max([len(examples.get(example_key(gap, le), [])) for le in les] + [1])
        fig, axes = plt.subplots(n_rows, len(les), squeeze=False,
                                 figsize=[2.4 * len(les), 0.75 * n_rows + 0.8])
        for column, le in zip(axes.T, les):
            key = example_key(gap, le)
            images = examples.get(key, [])
            replicas = examples.get(key + "_replicas", [])
            data = broken[np.isclose(broken.gap_nm, gap)
                          & np.isclose(broken.labelling_efficiency, le)]
            measured = data.set_index(data.replicate.astype(int))
            for k, ax in enumerate(column):
                if k >= len(images):
                    ax.axis("off")
                    continue
                show(ax, images[k])
                if k < len(replicas) and replicas[k] in measured.index:
                    annotate_break(ax, measured.loc[replicas[k]])
            if len(images):
                scale_bar(column[0])
            column[0].set_title(f"LE {le:g}\n{int(break_detected(data).sum())} / {len(data)} "
                                f"split", fontsize=9)
        fig.suptitle(f"Experimental actin filament, break of {gap:g} nm: labelling efficiency per "
                     f"column; green, break detected ({MIN_SEGMENTS_BROKEN} or more mask "
                     f"segments); scale bar 100 nm", x=0.02, ha="left", fontsize=10)
        fig.tight_layout()
        fig.savefig(prefix + f"examples_gap{gap:g}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def annotate_break(ax, found):
    """Segments and longest break of a filament image, framed green when its break is
    detected."""
    detected = break_detected(found)
    ax.text(0.02, 0.95, f"{found.mask_number_of_segments:.0f} segment(s), "
                        f"longest break {found.strip_longest_break:.0f} nm",
            color="C2" if detected else "w", fontsize=5, ha="left", va="top",
            transform=ax.transAxes)
    if detected:
        ax.spines[:].set_color("C2")
        ax.spines[:].set_linewidth(2.5)


def plot_models_grid(examples, scores, prefix):
    """The pair of images compared in every cell of the distinguishability map (intact above
    broken), laid out as the map, with the AUC of all features."""
    gaps = sorted(GAPS_NM, reverse=True)  # as the map: smallest break at the bottom
    les = sorted(LABELLING_EFFICIENCIES)
    ratios = ([1, 1, 0.35] * len(gaps))[:-1]  # a pair of panels per break, then a gap
    fig, axes = plt.subplots(len(ratios), len(les), squeeze=False,
                             figsize=[2.0 * len(les) + 0.8, 0.52 * sum(ratios) + 0.8],
                             gridspec_kw=dict(height_ratios=ratios, wspace=0.05, hspace=0.08))
    auc = scores[scores.feature_set == "all"]
    for i, gap in enumerate(gaps):
        intact_row, broken_row = axes[3 * i], axes[3 * i + 1]
        for k, le in enumerate(les):
            show(intact_row[k], examples[model_key("intact", 0, le)])
            show(broken_row[k], examples[model_key("broken", gap, le)])
            cell = auc[np.isclose(auc.gap_nm, gap) & np.isclose(auc.labelling_efficiency, le)]
            if len(cell):
                broken_row[k].text(0.98, 0.9, f"AUC {cell.auc.iloc[0]:.2f}", color="w",
                                   fontsize=6, ha="right", va="top",
                                   transform=broken_row[k].transAxes)
            if i == 0:
                intact_row[k].set_title(f"LE {le:g}", fontsize=8)
        intact_row[0].set_ylabel("intact", fontsize=6)
        broken_row[0].set_ylabel(f"{gap:g} nm", fontsize=6)
    for row in axes[2::3]:
        for ax in row:
            ax.axis("off")
    scale_bar(axes[-1, 0])
    fig.supylabel("Break length (nm)", fontsize=9)
    fig.supxlabel("Labelling efficiency", fontsize=9)
    fig.suptitle("Experimental actin filament: filaments compared in each cell of the "
                 "distinguishability map (intact above broken; AUC, all features); scale bar "
                 "100 nm", x=0.02, ha="left", fontsize=10)
    fig.savefig(prefix + "models_grid.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_figure(examples, detection, scores, path, panel_d=None, les=LABELLING_EFFICIENCIES):
    """The main figure: (a) one filament per break of FIGURE_EXAMPLE_GAPS_NM (rows; 0 nm: the intact
    filament) at every labelling efficiency of les, (b) the share of broken filaments whose
    break is detected, and the false alarm rate of intact filaments, (c) the AUC of intact
    against broken filaments from the number of mask segments alone, both against labelling
    efficiency, and (d), if given, drawn by panel_d(fig, grid cell)."""
    les = sorted(les)
    breaks = [gap for gap in FIGURE_EXAMPLE_GAPS_NM if gap > 0]  # one line each in b and c
    colours = {gap: f"C{k}" for k, gap in enumerate(breaks)}
    colours[0] = "0.5"  # the intact filament, as in b
    n_rows = 2 if panel_d is None else 3
    fig = plt.figure(figsize=[11, (3.9 if panel_d is None else 6.8)
                              + 0.55 * len(FIGURE_EXAMPLE_GAPS_NM)])
    grid = fig.add_gridspec(n_rows, 1, height_ratios=[0.25 * len(FIGURE_EXAMPLE_GAPS_NM) + 0.15, 1, 1][:n_rows],
                            hspace=0.55)
    top = grid[0].subgridspec(len(FIGURE_EXAMPLE_GAPS_NM), len(les), wspace=0.05, hspace=0.08)
    for i, gap in enumerate(FIGURE_EXAMPLE_GAPS_NM):
        model = "broken" if gap > 0 else "intact"
        for k, le in enumerate(les):
            ax = fig.add_subplot(top[i, k])
            show(ax, examples[model_key(model, gap, le)])
            if i == 0:
                ax.set_title(f"LE {le:.2g}", fontsize=9)
            if k == 0:
                ax.set_ylabel(f"{gap:g} nm" if gap > 0 else "no break", fontsize=8,
                              color=colours[gap])
            if i == 0 and k == 0:
                ax.set_title("a", loc="left", fontweight="bold", x=-0.12)
            if i == len(FIGURE_EXAMPLE_GAPS_NM) - 1 and k == 0:
                scale_bar(ax)
    bottom = grid[1].subgridspec(1, 2, wspace=0.3)

    ax = fig.add_subplot(bottom[0])
    for gap in breaks:
        data = detection[(detection.model == "broken") & np.isclose(detection.gap_nm, gap)] \
            .sort_values("labelling_efficiency")
        percent = 100 * data.fraction_split
        ax.errorbar(data.labelling_efficiency, percent,
                    yerr=error_bars(percent, 100 * data.fraction_low, 100 * data.fraction_high),
                    fmt="o-", color=colours[gap], capsize=3, label=f"{gap:g} nm break")
    intact = detection[detection.model == "intact"].sort_values("labelling_efficiency")
    percent = 100 * intact.fraction_split
    ax.errorbar(intact.labelling_efficiency, percent,
                yerr=error_bars(percent, 100 * intact.fraction_low, 100 * intact.fraction_high),
                fmt="s:", color="0.5", capsize=3, markersize=4, label="intact (false positive)")
    ax.set_ylim(-3, 103)
    ax.set_ylabel("Break detected (% of filaments)")
    ax.set_title("b", loc="left", fontweight="bold", x=-0.22)
    n_min, n_max = detection.n_filaments.min(), detection.n_filaments.max()
    ax.set_title(f"Mask split into {MIN_SEGMENTS_BROKEN}+ segments (n = "
                 f"{n_max if n_min == n_max else f'{n_min}-{n_max}'} per point)", fontsize=9)
    ax.legend(frameon=False, fontsize=7)

    ax_c = fig.add_subplot(bottom[1])
    data = scores[scores.feature_set == "segments"]
    for k, gap in enumerate(breaks):
        sub = data[np.isclose(data.gap_nm, gap)].sort_values("labelling_efficiency")
        ax_c.errorbar(sub.labelling_efficiency, sub.auc,
                      yerr=error_bars(sub.auc, sub.auc_low, sub.auc_high), fmt="o-",
                      color=colours[gap], capsize=3, label=f"{gap:g} nm break")
        ax_c.plot(sub.labelling_efficiency, sub.null_auc, "x", color=colours[gap],
                  label="null (intact vs intact)" if k == 0 else None)
    ax_c.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_c.set_ylim(min(0.3, np.nanmin(data.null_auc) - 0.05), 1.05)
    ax_c.set_ylabel("ROC AUC")
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.22)
    ax_c.set_title("Intact vs broken, number of mask segments", fontsize=9)
    ax_c.legend(frameon=False, fontsize=7, loc="lower right")
    for axis in (ax, ax_c):
        labelling_axis(axis, les)

    if panel_d is not None:
        panel_d(fig, grid[2])
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


### Experimental segments: intact vs broken from a frame split, without simulation

def render_segment(segment):
    """Image of the localisations of an aligned segment, as the simulated images: a field of
    SAMPLE_DIMENSIONS with FILAMENTS_PARAMS["pixel_size"] pixels, the filament centred along
    the rows, each localisation a Gaussian of RENDERING_SIGMA_NM."""
    pixel_size = FILAMENTS_PARAMS["pixel_size"]
    edges = [np.arange(-size / 2, size / 2 + pixel_size / 2, pixel_size)
             for size in SAMPLE_DIMENSIONS[:2]]
    image, __, __ = np.histogram2d(segment[:, 0], segment[:, 1], bins=edges)
    return ndimage.gaussian_filter(image, RENDERING_SIGMA_NM / pixel_size)


def segment_outline(center, direction, length=SEGMENT_LENGTH_NM, width=SEGMENT_WIDTH_NM):
    """Corners (xy, nm) of a cropped segment in the dataset."""
    across = np.array([-direction[1], direction[0]])
    return np.array([center + a * length / 2 * direction + b * width / 2 * across
                     for a, b in [(-1, -1), (1, -1), (1, 1), (-1, 1)]])


def measure_image(image, replicate=0):
    """measure_filaments of one rendered image (one row)."""
    return measure_filaments([dict(replicate=replicate, images={"SMLM": {"ch0": [image]}})],
                             **FILAMENTS_PARAMS)[0]


def classify_by_frame_split(localisations_all, crops):
    """
    The experimental segments: the candidates (their crops from candidate_table) labelled
    intact or broken from the localisations of their odd frames, and measured on those of
    their even frames.

    The localisations of each candidate are split by the parity of their frame (they are merged
    over consecutive frames, so a blink is in one half only). The odd half is rendered and
    measured: the segment is broken when its mask splits into MIN_SEGMENTS_BROKEN or more
    segments or its strip profiles show a break, intact when neither. The even half is rendered
    and measured, and these are the features scored: a real break is in both halves, the noise
    that set the label is not. Candidates with fewer than MIN_HALF_LOCALISATIONS localisations
    in either half are left out, and up to N_EXPERIMENTAL_PER_CLASS are kept per class, in the
    order of the candidates (most selective first), whatever their scored features.

    Returns the features (measure_filaments of the even half, one row per segment, with the
    features of the odd half prefixed label_) and the arrays of the segments: images (even
    half, cropped), centres, directions and models.
    """
    frames = localisations_all["frame"].to_numpy()
    labelled = []
    for k, (segment, index, centre, direction) in enumerate(tqdm(crops, desc="frame split")):
        odd = frames[index] % 2 == 1
        if min(odd.sum(), (~odd).sum()) < MIN_HALF_LOCALISATIONS:
            continue
        label = measure_image(render_segment(segment[odd])).iloc[0]
        broken = (label.mask_number_of_segments >= MIN_SEGMENTS_BROKEN
                  or label.strip_n_breaks >= 1)
        labelled.append(dict(model="broken" if broken else "intact", candidate=k, label=label,
                             segment=segment, odd=odd, centre=centre, direction=direction))
    rows, segments = [], dict(images=[], centres=[], directions=[], models=[])
    for model in ("intact", "broken"):
        found = [s for s in labelled if s["model"] == model]
        print(f"frame split: {len(found)} {model} segments, {min(len(found), N_EXPERIMENTAL_PER_CLASS)} kept")
        for k, s in enumerate(found[:N_EXPERIMENTAL_PER_CLASS]):
            image = render_segment(s["segment"][~s["odd"]])
            table = measure_image(image, replicate=k)
            rows.append(table.assign(
                model=model, candidate=s["candidate"], segment_x=s["centre"][0],
                segment_y=s["centre"][1], n_localisations=len(s["segment"]),
                n_label=int(s["odd"].sum()), n_scored=int((~s["odd"]).sum()),
                **{"label_" + feature: s["label"][feature] for feature in FEATURES}))
            segments["images"].append(crop_filament(image))
            segments["centres"].append(s["centre"])
            segments["directions"].append(s["direction"])
            segments["models"].append(model)
    features = pd.concat(rows, ignore_index=True)
    return features, {key: np.array(value) for key, value in segments.items()}


def experimental_break_detection(features):
    """Per class of experimental segments: the share whose mask splits into MIN_SEGMENTS_BROKEN
    or more segments, and the share with a break in the strip profiles (95% Wilson
    intervals), with the mean and s.d. of their segments, strip breaks and longest break."""
    rows = []
    for model in ("intact", "broken"):
        data = features[features.model == model]
        row = dict(model=model, n_segments=len(data))
        for name, found in [("mask_split", break_detected(data)),
                            ("strip_break", data.strip_n_breaks >= 1)]:
            k = int(found.sum())
            low, high = wilson_interval(k, len(data))
            row.update({f"n_{name}": k, f"fraction_{name}": k / len(data),
                        f"fraction_{name}_low": low, f"fraction_{name}_high": high})
        for column in ["mask_number_of_segments", "strip_n_breaks", "strip_longest_break",
                       "n_localisations"]:
            row[f"mean_{column}"] = data[column].mean()
            row[f"std_{column}"] = data[column].std()
        rows.append(row)
    return pd.DataFrame(rows)


def plot_experimental_figure(features, segments, detection, scores, prefix):
    """(a) The intact and broken experimental segments (and the reference segment, if chosen) on
    the dataset, (b) their images (even frames), with their segments and longest break, framed
    green when the break is detected, (c) the share of segments with a detected break per class
    and (d) the distinguishability of the classes per feature set. The classes are those of
    classify_by_frame_split: labelled from the odd frames, measured on the even frames."""
    colours = {"intact": "#2b8cbe", "broken": "#d7301f"}
    if "reference_centre" in segments:
        colours["reference"] = "#F58B19"
    models = segments["models"]
    fig = plt.figure(figsize=[15, 11])
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.5], height_ratios=[1.6, 1],
                            hspace=0.25, wspace=0.15)
    fig.suptitle("Experimental actin segments: labelled intact or broken from the localisations "
                 "of odd frames, imaged and scored on those of even frames", x=0.02, ha="left",
                 fontsize=10)

    ax = fig.add_subplot(grid[0, 0])
    overview, extent = segments["overview"], segments["overview_extent"]
    ax.imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
              extent=extent)
    outlines = []
    if "reference_centre" in segments:
        outlines.append(("reference", 0, segments["reference_centre"],
                         segments["reference_direction"]))
    counters = {"intact": 0, "broken": 0}
    for model, centre, direction in zip(models, segments["centres"], segments["directions"]):
        outlines.append((model, counters[model], centre, direction))
        counters[model] += 1
    for model, k, centre, direction in outlines:
        ax.add_patch(Polygon(segment_outline(centre, direction), closed=True, fill=False,
                             edgecolor=colours[model], linewidth=1.2))
        ax.text(centre[0] + 400, centre[1] - 400, "ref" if model == "reference" else f"{k}",
                color=colours[model], fontsize=6)
    for model in colours:
        ax.plot([], [], "s", mfc="none", color=colours[model], label=model)
    ax.legend(frameon=False, fontsize=7, labelcolor="w", loc="lower right")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.tick_params(labelsize=7)
    ax.set_title("a", loc="left", fontweight="bold")
    ax.set_title(f"Experimental segments ({SEGMENT_LENGTH_NM} x {SEGMENT_WIDTH_NM} nm)", fontsize=9)

    # intact segments above, broken below
    n_columns = 3
    n_rows = int(np.ceil(max(counters.values()) / n_columns))
    gallery = grid[0, 1].subgridspec(2 * n_rows + 1, n_columns, wspace=0.05, hspace=0.25,
                                     height_ratios=[1] * n_rows + [0.4] + [1] * n_rows)
    for model, offset in [("intact", 0), ("broken", n_rows + 1)]:
        index = np.flatnonzero(models == model)
        measured = features[features.model == model].set_index("replicate")
        for k, i in enumerate(index):
            ax_k = fig.add_subplot(gallery[offset + k // n_columns, k % n_columns])
            show(ax_k, segments["images"][i])
            annotate_break(ax_k, measured.loc[k])
            ax_k.text(0.98, 0.05, f"{k}", color=colours[model], fontsize=6, ha="right",
                      transform=ax_k.transAxes)
            if k == 0:
                ax_k.set_title(f"{model} ({len(index)})", fontsize=9, loc="left",
                               color=colours[model])
                scale_bar(ax_k)
                if model == "intact":
                    ax_k.text(-0.08, 1.15, "b", fontweight="bold", fontsize=12, ha="right",
                              transform=ax_k.transAxes)

    bottom = grid[1, :].subgridspec(1, 2, wspace=0.3)
    ax_c = fig.add_subplot(bottom[0])
    width = 0.38
    for k, (name, label) in enumerate([("mask_split", f"mask split ({MIN_SEGMENTS_BROKEN}+ segments)"),
                                       ("strip_break", "strip profiles (1+ breaks)")]):
        positions = np.arange(len(detection)) + (k - 0.5) * width
        percent = 100 * detection[f"fraction_{name}"]
        ax_c.bar(positions, percent, width, color=[colours[m] for m in detection.model],
                 alpha=1 if k == 0 else 0.5, hatch=None if k == 0 else "//",
                 edgecolor="k", linewidth=0.5, label=label,
                 yerr=error_bars(percent, 100 * detection[f"fraction_{name}_low"],
                                 100 * detection[f"fraction_{name}_high"]), capsize=3)
    ax_c.set_xticks(np.arange(len(detection)),
                    [f"{m}\n(n = {n})" for m, n in zip(detection.model, detection.n_segments)])
    ax_c.set_ylim(0, 105)
    ax_c.set_ylabel("Break detected (% of segments)")
    ax_c.legend(frameon=False, fontsize=7, loc="upper left")
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.set_title("c", loc="left", fontweight="bold", x=-0.15)

    ax_d = fig.add_subplot(bottom[1])
    positions = np.arange(len(scores))
    ax_d.errorbar(positions, scores.auc, yerr=error_bars(scores.auc, scores.auc_low, scores.auc_high),
                  fmt="o", color="k", capsize=3, label="intact vs broken")
    ax_d.plot(positions, scores.null_auc, "x", color="0.5", label="null (intact vs intact)")
    for x, (auc, p) in enumerate(zip(scores.auc, scores.p_value)):
        ax_d.text(x + 0.08, auc, f"p = {p:.3g}", fontsize=6, va="center")
    ax_d.axhline(0.5, color="0.5", linestyle=":", linewidth=1)
    ax_d.set_xticks(positions, scores.feature_set)
    ax_d.set_xlim(-0.5, len(scores) - 0.3)
    ax_d.set_ylim(min(0.3, np.nanmin(scores.null_auc) - 0.05), 1.05)
    ax_d.set_ylabel("ROC AUC")
    ax_d.legend(frameon=False, fontsize=7, loc="lower right")
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.set_title("d", loc="left", fontweight="bold", x=-0.15)
    ax_d.set_title("Intact vs broken experimental segments", fontsize=9)
    fig.savefig(prefix + "experimental_segments_figure.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def experimental_analysis(features, segments, prefix):
    """Break detection and distinguishability of the experimental segments, with their
    figure."""
    counts = features.model.value_counts()
    if min(counts.get("intact", 0), counts.get("broken", 0)) < 2:
        print(f"Experimental segments not scored: {counts.to_dict()} (2 or more per class needed)")
        return
    features = features.fillna(FILL_VALUES)
    detection = experimental_break_detection(features)
    detection.to_csv(prefix + "experimental_break_detection.csv", index=False)
    scores = score_feature_sets(features, progress=False)
    scores.to_csv(prefix + "experimental_distinguishability.csv", index=False)
    plot_experimental_figure(features, segments, detection, scores, prefix)
    print(detection.round(3).to_string(index=False))
    print(scores[["feature_set", "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc",
                  "p_value"]].round(3).to_string(index=False))


### Experimental labelling efficiency: sets of segments of decreasing density, with cut gaps

def cut_gap(segment, gap_nm):
    """Localisations of an aligned segment (filament along x, centred) without those within
    gap_nm / 2 of its centre."""
    return segment[np.abs(segment[:, 0]) >= gap_nm / 2]


def candidate_table(localisations_all):
    """
    The candidate segments (CANDIDATE_SEGMENTS), each cropped, rendered and measured
    (measure_filaments) as the images of the sets.

    Returns one row per candidate, in the order of CANDIDATE_SEGMENTS (the most selective
    first): the centre (segment_x, segment_y) and direction (direction_x, direction_y) of the
    crop, n_localisations, and the mask segments and strip breaks of its image; and the crop of
    every candidate (crop_filament_segment, with the index of the localisations kept instead of
    their mask).
    """
    centres = np.array([(x, y) for x, y, __ in CANDIDATE_SEGMENTS], dtype=float)
    angles = np.deg2rad([angle for __, __, angle in CANDIDATE_SEGMENTS])
    candidates = pd.DataFrame(dict(segment_x=centres[:, 0], segment_y=centres[:, 1],
                                   direction_x=np.cos(angles), direction_y=np.sin(angles)))
    crops, n_localisations, segments, breaks = [], [], [], []
    for k, row in enumerate(tqdm(candidates.itertuples(), total=len(candidates),
                                 desc="candidate segments")):
        segment, keep, centre, direction = crop_filament_segment(
            localisations_all, (row.segment_x, row.segment_y), SEGMENT_LENGTH_NM,
            SEGMENT_WIDTH_NM, direction=(row.direction_x, row.direction_y))
        crops.append((segment, np.flatnonzero(keep), centre, direction))
        measured = measure_image(render_segment(segment), k).fillna(FILL_VALUES).iloc[0]
        n_localisations.append(len(segment))
        segments.append(measured.mask_number_of_segments)
        breaks.append(measured.strip_n_breaks)
    return candidates.assign(n_localisations=n_localisations, mask_segments=segments,
                             strip_breaks=breaks), crops


def select_density_sets(candidates):
    """
    Sets of experimental segments of decreasing density, from the candidates
    (candidate_table).

    Set 0 is the N_SEGMENTS_PER_SET densest candidates without a detected break (mask and strip
    profiles); set k the up to N_SEGMENTS_PER_SET others closest to DENSITY_TARGETS[k] times its
    mean number of localisations, within DENSITY_TOLERANCE.

    Returns the selection (one row per segment: candidate, set, target, the centre and direction
    of the crop) and the candidates with their set (-1 if none).
    """
    candidates = candidates.assign(set=-1)
    no_break = candidates[(candidates.mask_segments < MIN_SEGMENTS_BROKEN)
                          & (candidates.strip_breaks == 0)]
    reference = no_break.sort_values("n_localisations", ascending=False).head(N_SEGMENTS_PER_SET)
    candidates.loc[reference.index, "set"] = 0
    reference_mean = reference.n_localisations.mean()
    selection = [reference.assign(set=0, target=1.0)]
    for target in DENSITY_TARGETS[1:]:
        pool = candidates[candidates.set == -1]
        distance = np.abs(pool.n_localisations - target * reference_mean)
        within = distance[distance <= DENSITY_TOLERANCE * target * reference_mean]
        picked = pool.loc[within.sort_values().index[:N_SEGMENTS_PER_SET]]
        if len(picked) < N_SEGMENTS_PER_SET:
            print(f"Warning: {len(picked)} of {N_SEGMENTS_PER_SET} segments within "
                  f"{100 * DENSITY_TOLERANCE:.0f}% of the target {target:g}")
        if len(picked) < DENSITY_MIN_SEGMENTS:
            print(f"Warning: target {target:g} left out (fewer than {DENSITY_MIN_SEGMENTS} segments)")
            continue
        s = len(selection)
        candidates.loc[picked.index, "set"] = s
        selection.append(picked.assign(set=s, target=target))
    selection = pd.concat(selection).rename_axis("candidate").reset_index()
    for (s, target), data in selection.groupby(["set", "target"]):
        print(f"set {s} (target {target:g}): {len(data)} segments, "
              f"{data.n_localisations.mean():.0f} ± {data.n_localisations.std():.0f} localisations")
    return selection, candidates


def measure_density_sets(crops, selection):
    """The segments of the selection (select_density_sets; their crops from candidate_table),
    intact and with every gap of DENSITY_GAPS_NM cut: their features (measure_filaments, one
    row per image), with the labelling efficiency of their set (its mean number of
    localisations over that of the reference set), and the arrays of the sets (images of the
    intact segments, examples of the figure, centres, directions, numbers of localisations,
    targets)."""
    segments = []
    for s, data in selection.groupby("set"):
        for k, row in enumerate(data.itertuples()):
            segment, __, segment_centre, direction = crops[row.candidate]
            segments.append((s, k, segment, segment_centre, direction))
    set_index = np.array([s for s, *__ in segments])
    n_localisations = np.array([len(segment) for __, __, segment, __, __ in segments])
    n_sets = selection.set.nunique()
    means = np.array([n_localisations[set_index == s].mean() for s in range(n_sets)])
    les = np.round(means / means[0], 3)
    rows, images, examples = [], [], {}
    for s, k, segment, centre, __ in tqdm(segments, desc="experimental segments"):
        for model, gap in [("intact", 0)] + [("broken", gap) for gap in DENSITY_GAPS_NM]:
            image = render_segment(cut_gap(segment, gap))
            table = measure_image(image, replicate=k)
            rows.append(table.assign(model=model, gap_nm=gap, labelling_efficiency=les[s], set=s,
                                     n_localisations=len(segment), segment_x=centre[0],
                                     segment_y=centre[1]))
            if model == "intact" or k == 0:
                cropped = crop_filament(image)
            if model == "intact":
                images.append(cropped)
            if k == 0:
                examples["example_" + model_key(model, gap, les[s])] = cropped
    sets = dict(examples, images=np.array(images), set=set_index, segment=np.array([k for __, k, *__ in segments]),
                n_localisations=n_localisations, labelling_efficiency=les,
                centres=np.array([centre for *__, centre, __ in segments]),
                directions=np.array([direction for *__, direction in segments]),
                targets=selection.groupby("set").target.first().to_numpy())
    return pd.concat(rows, ignore_index=True), sets


def density_set_table(features, sets):
    """One row per segment of the sets: its set, labelling efficiency, centre, direction,
    number of localisations, and the mask segments and strip breaks of the intact segment."""
    intact = features[features.model == "intact"].set_index(["set", "replicate"])
    rows = []
    for s, k, n, centre, direction in zip(sets["set"], sets["segment"], sets["n_localisations"],
                                          sets["centres"], sets["directions"]):
        measured = intact.loc[(s, k)]
        rows.append(dict(set=int(s), segment=int(k), target=sets["targets"][s],
                         labelling_efficiency=sets["labelling_efficiency"][s],
                         centre_x_nm=centre[0], centre_y_nm=centre[1],
                         direction_x=direction[0], direction_y=direction[1], n_localisations=int(n),
                         intact_mask_segments=int(measured.mask_number_of_segments),
                         intact_strip_breaks=int(measured.strip_n_breaks)))
    return pd.DataFrame(rows)


def density_panel(table):
    """Panel d of the experimental figure: the localisations of every segment of each set,
    with their mean +- s.d. and the labelling efficiency of the set."""
    def panel(fig, cell):
        ax = fig.add_subplot(cell)
        rng = np.random.default_rng(0)
        colours = density_colours(table.set.nunique())
        for s, data in table.groupby("set"):
            ax.scatter(s + rng.uniform(-0.15, 0.15, len(data)), data.n_localisations, s=12,
                       color=colours[s])
            mean, std = data.n_localisations.mean(), data.n_localisations.std()
            ax.errorbar(s + 0.3, mean, yerr=std, fmt="o", c="k", capsize=3, markersize=4)
            ax.text(s, data.n_localisations.max() + 25,
                    f"LE {data.labelling_efficiency.iloc[0]:.2g}\n{mean:.0f} ± {std:.0f}\n"
                    f"n = {len(data)}",
                    ha="center", fontsize=7)
        ax.set_xticks(range(table.set.nunique()),
                      ["reference" if s == 0 else f"set {s}" for s in range(table.set.nunique())])
        ax.set_ylim(0, table.n_localisations.max() * 1.35)
        ax.set_ylabel("Localisations per segment")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("d", loc="left", fontweight="bold", x=-0.1)
        ax.set_title("Experimental labelling efficiency: localisations of the segments of each "
                     "set (mean ± s.d.) over those of the reference set", fontsize=9)
    return panel


def density_colours(n_sets):
    """A colour per set, visible on the dark overview."""
    return plt.cm.plasma(np.linspace(0.15, 0.95, n_sets))


def plot_density_sets(table, sets, path):
    """The segments of every set on the dataset, and the image of every intact segment (a row
    per set)."""
    n_sets = table.set.nunique()
    # at most 10 images per row, as many rows per set as needed
    n_columns = min(table.groupby("set").size().max(), 10)
    rows_per_set = table.groupby("set").size().apply(lambda n: int(np.ceil(n / n_columns)))
    first_row = np.concatenate([[0], np.cumsum(rows_per_set)[:-1]])
    colours = density_colours(n_sets)
    fig = plt.figure(figsize=[7 + 1.6 * n_columns, max(7, 0.75 * rows_per_set.sum())])
    grid = fig.add_gridspec(1, 2, width_ratios=[7, 1.6 * n_columns], wspace=0.08)
    ax = fig.add_subplot(grid[0])
    overview = sets["overview"]
    ax.imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
              extent=sets["overview_extent"])
    for row in table.itertuples():
        centre = np.array([row.centre_x_nm, row.centre_y_nm])
        ax.add_patch(Polygon(segment_outline(centre, np.array([row.direction_x, row.direction_y])),
                             closed=True, fill=False, edgecolor=colours[row.set], linewidth=1.2))
        ax.text(centre[0] + 300, centre[1] - 300, f"{row.set}.{row.segment}", color=colours[row.set],
                fontsize=5)
    for s in range(n_sets):
        ax.plot([], [], "s", mfc="none", color=colours[s],
                label=f"{'reference' if s == 0 else f'set {s}'} (LE {sets['labelling_efficiency'][s]:.2g})")
    ax.legend(frameon=False, fontsize=7, labelcolor="w", loc="lower right")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.tick_params(labelsize=7)
    gallery = grid[1].subgridspec(rows_per_set.sum(), n_columns, wspace=0.05, hspace=0.35)
    for i, row in enumerate(table.itertuples()):
        ax_k = fig.add_subplot(gallery[first_row[row.set] + row.segment // n_columns,
                                       row.segment % n_columns])
        show(ax_k, sets["images"][i])
        ax_k.set_title(f"{row.set}.{row.segment}: {row.n_localisations} loc., "
                       f"{row.intact_mask_segments} seg., {row.intact_strip_breaks} br.",
                       fontsize=5, color=colours[row.set] if row.set else "k")
        if row.segment == 0:
            ax_k.set_ylabel("reference" if row.set == 0 else f"set {row.set}", fontsize=7)
    fig.suptitle("Experimental labelling efficiency: segments of each set (set.segment; "
                 "localisations, mask segments and strip breaks of the intact segment)",
                 x=0.02, ha="left", fontsize=10)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_heatmap(scores, path, title):
    """Distinguishability heatmaps (ROC AUC) per feature set, break length against labelling
    efficiency."""
    fig = plot_distinguishability_heatmaps(scores, x="labelling_efficiency", y="gap_nm",
                                           facet="feature_set", facet_order=list(FEATURE_SETS))
    n_le = scores.labelling_efficiency.nunique()
    gaps = sorted(scores.gap_nm.unique())
    panels = fig.axes[:len(FEATURE_SETS)]
    for ax in panels:
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Break length (nm)")
    fig.suptitle(title, x=0.02, ha="left")
    fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * len(gaps))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def density_analysis(features, sets, folder):
    """Break detection and distinguishability of intact and cut experimental segments at every
    experimental labelling efficiency, with their figures, tables and README.txt, in folder."""
    features = features.fillna(FILL_VALUES)
    table = density_set_table(features, sets)
    table.to_csv(folder + "sets.csv", index=False)
    plot_density_sets(table, sets, folder + "sets.png")
    paired = feature_table([features], DENSITY_GAPS_NM)
    detection = break_detection(paired)
    detection.to_csv(folder + "break_detection.csv", index=False)
    scores = score_feature_sets(paired, by=CONDITIONS, n_workers=N_WORKERS,
                                label="experimental labelling efficiency, features")
    scores.to_csv(folder + "distinguishability.csv", index=False)
    figure_values(detection, scores).to_csv(folder + "figure_values.csv", index=False)
    save_heatmap(scores, folder + "distinguishability_heatmap.png",
                 "Experimental segments: intact vs cut (ROC AUC), experimental labelling efficiency")
    examples = {key[len("example_"):]: value for key, value in sets.items()
                if key.startswith("example_")}
    plot_figure(examples, detection, scores, folder + "filament_figure.png", density_panel(table),
                les=list(sets["labelling_efficiency"]))
    write_density_readme(folder, table, detection, scores,
                         int(sets.get("n_segments_per_set", N_SEGMENTS_PER_SET)))
    print(detection.round(3).to_string(index=False))


def write_density_readme(folder, table, detection, scores, requested):
    """README.txt of the experimental labelling efficiency folder: what every file holds, with
    the sets of this run (requested: segments asked for per set)."""
    summary = table.groupby("set").agg(
        n=("segment", "size"), target=("target", "first"),
        labelling_efficiency=("labelling_efficiency", "first"),
        mean=("n_localisations", "mean"), std=("n_localisations", "std"),
        minimum=("n_localisations", "min"), maximum=("n_localisations", "max"),
        mask_split=("intact_mask_segments", lambda x: int((x >= MIN_SEGMENTS_BROKEN).sum())),
        strip_break=("intact_strip_breaks", lambda x: int((x > 0).sum())))
    set_lines = "\n".join(
        f"  {'reference' if s == 0 else f'set {s}':9s}  LE {row['labelling_efficiency']:.3f}  "
        f"(target {row['target']:.2f})  {row['n']:.0f} segments"
        f"{'' if row['n'] >= requested else f' (of {requested} asked for)'}  "
        f"{row['mean']:6.1f} ± {row['std']:5.1f} localisations "
        f"(range {row['minimum']:.0f}-{row['maximum']:.0f}); intact segments with a mask split: "
        f"{row['mask_split']:.0f}/{row['n']:.0f}, with a strip break: {row['strip_break']:.0f}/{row['n']:.0f}"
        for s, row in summary.iterrows())
    missing = [f"{t:g}" for t in DENSITY_TARGETS if not np.isclose(summary.target, t).any()]
    segments_auc = scores[scores.feature_set == "segments"].pivot_table(
        index="gap_nm", columns="labelling_efficiency", values="auc").round(2).to_string()
    detected = detection.pivot_table(index=["model", "gap_nm"], columns="labelling_efficiency",
                                     values="fraction_split").mul(100).round(0).to_string()
    text = f"""EXPERIMENTAL LABELLING EFFICIENCY
=================================

Written by figS_filaments_actin3D_smlm.py. Nothing in this folder is simulated: every
image is rendered from experimental 3D SMLM localisations of actin labelled with phalloidin
(ShareLoc.xyz, Zenodo 5534337) and measured as the simulated images of the parent folder.

What was done
-------------
1. Segments. The {len(CANDIDATE_SEGMENTS)} candidate segments of CANDIDATE_SEGMENTS (centre and direction in the
   dataset), selected once from the dataset: straight (the localisations
   concentrate in a narrow strip along the filament and their position across it varies little
   along it), isolated (few localisations away from the axis, i.e. no crossing filament),
   without a cluster of localisations along it, and not overlapping. Each is cropped as the
   reference segment of the parent folder ({SEGMENT_LENGTH_NM} x {SEGMENT_WIDTH_NM} nm, along x across, aligned to the
   x axis; crop_filament_segment). The candidates are listed in candidates.csv.
2. Sets. {requested} segments were asked for per set (N_SEGMENTS_PER_SET; set the environment
   variable VLAB4MIC_N_SEGMENTS_PER_SET to change it). The reference set is the densest
   candidates where neither the mask (fewer than {MIN_SEGMENTS_BROKEN} segments) nor the strip profiles
   (no flat strip between filament) detect a break. The other sets are the candidates
   closest to {', '.join(f'{t:g}' for t in DENSITY_TARGETS[1:])} times the mean number of localisations of the
   reference set, within {100 * DENSITY_TOLERANCE:.0f}% of that number, chosen by density only (they may show
   apparent breaks); a set has fewer segments than asked for when there are not enough
   candidates within the tolerance, and is left out below {DENSITY_MIN_SEGMENTS}.{f' Left out: targets {", ".join(missing)}.' if missing else ''}
   The labelling efficiency (LE) of a set is its mean number of localisations over that of
   the reference set: a density proxy, not a labelling efficiency simulated in VLab4Mic.

{set_lines}

3. Breaks. Every segment is measured intact and with a gap cut in its middle: the
   localisations within gap/2 of its centre are removed (cut_gap), for gaps of
   {', '.join(f'{g:g}' for g in DENSITY_GAPS_NM)} nm. There is no 0 nm gap: without a second realisation it would be
   the intact segment itself.
4. Images. The localisations are rendered as the simulated images: {FILAMENTS_PARAMS['pixel_size']} nm pixels,
   each localisation a Gaussian of {RENDERING_SIGMA_NM} nm, the filament centred along the rows of a
   {SAMPLE_DIMENSIONS[0]} x {SAMPLE_DIMENSIONS[1]} nm field (render_segment).
5. Measures. measure_filaments, as in the parent folder: filament widths, mask (Gaussian
   smoothing of {FILAMENTS_PARAMS['smooth_sigma']} px, {FILAMENTS_PARAMS['threshold_method']} threshold) continuity, number of segments and eccentricities,
   and strip profiles ({FILAMENTS_PARAMS['strip_kwargs']['strip_length']} nm strips along the axis; a strip is flat below
   {FILAMENTS_PARAMS['strip_kwargs']['flat_fraction']:g} of the reference amplitude, and flat strips between filament are a break).
6. Break detection. A break is detected when the mask splits into {MIN_SEGMENTS_BROKEN} or more
   segments; the share of segments where it is, with a 95% Wilson interval, is given for
   the cut segments (detected breaks) and for the intact ones (false alarms).
7. Distinguishability. For every gap and LE, the intact segments of the set are
   compared with the same segments cut (distinguishability_from_features): out-of-fold
   ROC AUC of a classifier with a 95% bootstrap interval, a null control (two random halves
   of the intact segments, expected 0.5) and a permutation p-value ({N_PERMUTATIONS} permutations),
   for each feature set: segments (mask_number_of_segments alone), mask, strip profiles,
   all.

Break detected (% of segments, mask split), rows model and gap (nm), columns LE:
{detected}

AUC, intact vs cut, number of mask segments alone; rows gap (nm), columns LE:
{segments_auc}

Caveats
-------
- The two classes are the same segments, intact and cut (paired); the out-of-fold
  scoring does not use the pairing.
- LE is the density of localisations relative to the reference set. Sparser segments
  may also differ in labelling, blinking, axial position or crowding, not only in the
  number of labels.
- Phalloidin labels every filamentous actin: a segment may be a bundle of filaments or part
  of a stress fibre rather than a single filament, and denser segments are more likely to
  be bundles.

Files
-----
filament_figure.png
    The filament figure of the parent folder with the experimental LE:
    (a) one segment per LE (columns: the first segment of each set), intact (top row,
        "no break") and with a gap of {', '.join(f'{g:g}' for g in FIGURE_EXAMPLE_GAPS_NM if g > 0)} nm cut (rows); {EXAMPLE_LENGTH_NM} x {EXAMPLE_WIDTH_NM} nm, scale
        bar 100 nm.
    (b) Break detected (% of segments whose mask splits into {MIN_SEGMENTS_BROKEN} or more segments,
        95% Wilson interval, n = the segments of the set) against LE, one line per gap, and for the
        intact segments (grey, false positives).
    (c) ROC AUC of intact vs cut segments from the number of mask segments alone (95%
        bootstrap interval) against LE, one line per gap; crosses: null control.
    (d) Localisations per segment of each set (dots), their mean +- s.d. (black), and the
        LE of the set.
sets.png
    Left: the dataset with the segments of every set (boxes, labelled set.segment, colour
    per set). Right: the image of every intact segment, one row per set, with its number
    of localisations, mask segments and strip breaks.
distinguishability_heatmap.png
    ROC AUC of intact vs cut, break length against LE, one panel per feature set.
sets.csv
    One row per segment: set (0: reference), segment (index in the set), target (fraction
    of the localisations of the reference set aimed at), labelling_efficiency, centre_x_nm and centre_y_nm (centre of the crop in the dataset),
    direction_x and direction_y (filament direction), n_localisations (in the {SEGMENT_LENGTH_NM} x {SEGMENT_WIDTH_NM} nm
    crop), intact_mask_segments and intact_strip_breaks (of the intact segment).
features.csv
    One row per image (each segment intact, gap_nm 0, and with every gap cut): the
    measure_filaments columns (mean/std/median_filament_width and s_width, nm; mask_continuity,
    mask_major_lenght, mask_number_of_segments, min/max/mean_eccentricity; strip_fraction_flat,
    strip_n_breaks, strip_longest_break (nm), strip_min_amplitude; and settings:
    mask_prefilter_smooth_sigma, pixelsize, row_binning), with replicate (segment index in
    its set), model (intact or broken), gap_nm, labelling_efficiency, set, n_localisations
    (of the intact segment), segment_x and segment_y (nm).
break_detection.csv
    One row per model, gap and LE: n_filaments (segments), n_split and fraction_split
    (mask split into {MIN_SEGMENTS_BROKEN}+ segments), fraction_low and fraction_high (95% Wilson
    interval), median/mean/std_segments (number of mask segments). The intact segments are
    listed once, with gap_nm 0.
distinguishability.csv
    One row per feature set, gap and LE: features, n_a and n_b (intact and cut segments),
    accuracy and auc with their 95% intervals (_low, _high), null_auc, p_value, and for the
    feature set "all" the AUC of every feature alone (<feature>_auc, _low, _high).
figure_values.csv
    The values of b and c: one row per model (intact, and cut with every gap of the figure)
    and LE; b_percent_detected (with b_percent_low, b_percent_high), b_mean_segments,
    b_std_segments; c_auc (with c_auc_low, c_auc_high), c_null_auc, c_p_value.
candidates.csv
    One row per candidate segment (CANDIDATE_SEGMENTS: straight, isolated, non-overlapping),
    the most selective first: segment_x and segment_y (centre of the crop, nm), direction_x and
    direction_y, n_localisations, mask_segments and strip_breaks (of its image), set (-1: not
    in a set).
sets.npz
    Arrays to plot the figures again without the localisations (--features): images (the
    intact segments, in the order of sets.csv), set, segment, n_localisations, centres,
    directions, labelling_efficiency and targets (per set), n_segments_per_set (asked for),
    example_<model>_gap_<gap>_labelling_efficiency_<LE>
    (the images of panel a), overview and overview_extent (the dataset, 40 nm pixels).
README.txt
    This file.
"""
    with open(folder + "README.txt", "w") as handle:
        handle.write(text)


### Reference segment: overview and targets

def dataset_overview(localisations_all):
    """The localisations in 40 nm pixels, and the extent of the image (imshow, origin upper)."""
    overview, xedges, yedges = np.histogram2d(
        localisations_all.x, localisations_all.y, bins=int(localisations_all.x.max() / 40))
    return overview.astype(np.float32), np.array([xedges[0], xedges[-1], yedges[-1], yedges[0]])


def plot_overview(overview, extent, candidates, prefix, localisations_all, in_segment,
                  fluorophore_positions):
    """The dataset with the candidate segments (numbered as in candidates.csv, to choose
    REFERENCE_CANDIDATE) and the reference segment, cropped and aligned to x."""
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    axs[0].imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
                  extent=extent)
    for k, row in enumerate(candidates.itertuples()):
        centre = np.array([row.segment_x, row.segment_y])
        axs[0].add_patch(Polygon(segment_outline(centre, np.array([row.direction_x, row.direction_y])),
                                 closed=True, fill=False, edgecolor="#2b8cbe", linewidth=0.6))
        axs[0].text(centre[0] + 250, centre[1] - 250, f"{k}", color="#2b8cbe", fontsize=4)
    axs[0].set_title("Experimental localisations: candidate segments (blue, numbered as in "
                     "candidates.csv), reference segment in orange")
    axs[0].scatter(localisations_all.x.iloc[in_segment], localisations_all.y.iloc[in_segment],
                   s=0.5, c="#F58B19")
    zoom, xedges, yedges = np.histogram2d(
        fluorophore_positions[:, 0], fluorophore_positions[:, 1],
        bins=[np.arange(-SEGMENT_LENGTH_NM / 2, SEGMENT_LENGTH_NM / 2 + 5, 5),
              np.arange(-SEGMENT_WIDTH_NM / 2, SEGMENT_WIDTH_NM / 2 + 5, 5)]
    )
    axs[1].imshow(zoom.T, cmap="grey", origin="lower",
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    axs[1].set_title("Cropped segment, aligned to x (nm)")
    fig.savefig(prefix + "segment_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_targets(particle_data, prefix):
    """The targets of the reference segment, from above and from the side."""
    experiment = build_experiment(1, RANDOM_SEED, particle_data, run_simulation=False,
                                  modules=("coordinate_field",))
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    experiment.coordinate_field.show_field(view_init=[90, 0, 0], emitters_plotsize=1, axis_object=ax1)
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    experiment.coordinate_field.show_field(view_init=[0, 0, 0], emitters_plotsize=1, axis_object=ax2)
    fig.savefig(prefix + "filament.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
        examples = {}
        if os.path.exists(prefix + "examples.npz"):
            examples = dict(np.load(prefix + "examples.npz"))
        if os.path.exists(prefix + "experimental_features.csv") \
                and os.path.exists(prefix + "experimental_segments.npz"):
            experimental_analysis(pd.read_csv(prefix + "experimental_features.csv"),
                                  dict(np.load(prefix + "experimental_segments.npz")), prefix)
        density_dir = os.path.join(prefix + DENSITY_FOLDER, "")
        if os.path.exists(density_dir + "features.csv") and os.path.exists(density_dir + "sets.npz"):
            density_analysis(pd.read_csv(density_dir + "features.csv"),
                             dict(np.load(density_dir + "sets.npz")), density_dir)
    else:
        output_dir = os.path.join(str(output_path),
                                  ExperimentParametrisation().date_as_string + OUTPUT_NAME)
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, "")  # the outputs are named <output_dir>/<name>
        images_dir = prefix + "images" if SAVE_IMAGES else None

        localisations_all = load_smlm_table(SMLM_EXPERIMENTAL_PATH)
        overview, overview_extent = dataset_overview(localisations_all)
        candidates, crops = candidate_table(localisations_all)
        candidates.to_csv(prefix + "candidates.csv", index_label="candidate")
        fluorophore_positions, in_segment, reference_centre, reference_direction = \
            crops[REFERENCE_CANDIDATE]
        print(f"{fluorophore_positions.shape[0]} localisations used as fluorophore positions")
        plot_overview(overview, overview_extent, candidates, prefix, localisations_all,
                      in_segment, fluorophore_positions)

        # experimental segments, without simulation: labelled from their odd frames, scored on
        # their even frames
        experimental_features, segments = classify_by_frame_split(localisations_all, crops)
        segments.update(overview=overview, overview_extent=overview_extent,
                        reference_centre=reference_centre, reference_direction=reference_direction)
        experimental_features.to_csv(prefix + "experimental_features.csv", index=False)
        np.savez_compressed(prefix + "experimental_segments.npz", **segments)
        experimental_analysis(experimental_features, segments, prefix)

        # experimental labelling efficiency: sets of segments of decreasing density, in a subfolder
        density_dir = os.path.join(prefix + DENSITY_FOLDER, "")
        os.makedirs(density_dir, exist_ok=True)
        selection, density_candidates = select_density_sets(candidates)
        density_candidates.to_csv(density_dir + "candidates.csv", index=False)
        density_features, sets = measure_density_sets(crops, selection)
        sets.update(overview=overview, overview_extent=overview_extent,
                    n_segments_per_set=N_SEGMENTS_PER_SET)
        density_features.to_csv(density_dir + "features.csv", index=False)
        np.savez_compressed(density_dir + "sets.npz", **sets)
        density_analysis(density_features, sets, density_dir)

        template = experimental_template(fluorophore_positions)
        plot_targets(template, prefix)
        plot_strip_examples(template, prefix)
        models = [("intact", 0, template)]
        n_targets = len(fluorophore_positions)
        for gap in GAPS_NM:
            broken = break_filament(template, gap)
            n = len(broken["targets"][TARGET]["coordinates"])
            print(f"gap {gap} nm: {n} of {n_targets} localisations")
            models.append(("broken", gap, broken))
        blocks = [(model, gap, le) for model, gap, __ in models for le in LABELLING_EFFICIENCIES]
        data = {(model, gap): particle_data for model, gap, particle_data in models}
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [(model, gap, le, int(seed), data[(model, gap)], images_dir)
                for (model, gap, le), seed in zip(blocks, seeds)]
        print(f"{len(jobs)} simulations, {N_REALISATIONS} realisations each, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                outputs = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="simulations"))
        else:
            outputs = [simulate(job) for job in tqdm(jobs, desc="simulations")]
        features = feature_table([table for table, __ in outputs])
        features.to_csv(prefix + "features.csv", index=False)
        examples = {key: images for __, found in outputs for key, images in found.items()}
        examples.update(model_epitopes(models))
        np.savez_compressed(prefix + "examples.npz", **examples)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict, tuple))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    plot_feature_distributions(features, prefix)
    features = features.fillna(FILL_VALUES)
    if examples:
        plot_examples(features, examples, prefix)
    scores = score_feature_sets(features, by=CONDITIONS, n_workers=N_WORKERS)
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    save_heatmap(scores, prefix + "distinguishability_heatmap.png",
                 "Experimental actin filament: intact vs broken filament (ROC AUC)")
    plot_auc_vs_labelling(scores, prefix)
    detection = break_detection(features)
    detection.to_csv(prefix + "break_detection.csv", index=False)
    figure_values(detection, scores).to_csv(prefix + "filament_figure_values.csv", index=False)
    print(detection.round(3).to_string(index=False))
    if epitope_key("intact", 0) in examples:  # the targets of the models, from examples.npz
        plot_models(examples, prefix)
        plot_models_grid(examples, scores, prefix)
        plot_figure(examples, detection, scores, prefix + "filament_figure.png")
    print(scores.pivot_table(index=["feature_set", "gap_nm"], columns="labelling_efficiency",
                             values="auc", sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
