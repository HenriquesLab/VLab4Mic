"""Fig. 2, actin: can a break in an actin filament be detected, and at which labelling efficiency?

The counterpart of fig2_filament_break_distinguishability.py (microtubule) for F-actin labelled
with phalloidin.

The filament. 6T25 (phalloidin-Alexa Fluor 546 stabilised F-actin) is a segment of 5 actin
subunits (chains A-E), each with its phalloidin (chains F-J). The target is the terminal
5-hydroxyl (atom O2) of the 4,5-dihydroxyleucine of every phalloidin (residue EEP), where the
dye is conjugated; a linker of length 0 places the fluorophore on it (direct labelling). 
Superposing every actin subunit onto the next gives one screw about the
filament axis (166.86 degrees, 2.752 nm). Actin is a one-start helix that this screw extends
without gaps or clashes, but no pure translation repeats it. The filament is the EEP site of the
first subunit moved by 0, 1, ... N - 1 steps of the screw, N = FILAMENT_LENGTH_NM / rise + 1
(182 sites), laid flat (axis along x) and loaded as one particle of the virtual sample.

Break detection and distinguishability, as fig2_filament_break_distinguishability.py. A second
model has an artificial break in the middle: every EEP site within gap/2 of the centre of the
filament is removed. Each model is simulated with run_replicates at every labelling efficiency
(of the phalloidin) and every image is measured with measure_filaments (widths, mask continuity,
number of segments, eccentricities, strip profiles), giving a feature table with one row per
image. A break is detected in an image when its mask splits into MIN_SEGMENTS_BROKEN or more
segments; the share of broken filaments where it is detected is given with a 95% Wilson
interval, and the same share of intact filaments is its false alarm rate. The table is scored
per break size and labelling efficiency with distinguishability_from_features (out-of-fold ROC
AUC of intact vs broken, 95% bootstrap interval, null control and permutation p-value) for each
feature set (segments, mask, strip profiles, all). A gap of 0 nm is a second, independently
simulated intact filament: a negative control (expected 0.5).

Outputs, in the folder <date>fig2_filaments_actin of the VLab4Mic output directory (created if
it does not exist), as fig2_filament_break_distinguishability.py: filament_figure.png (the main
figure), features.csv, distinguishability.csv, break_detection.csv, filament_figure_values.csv,
settings.json, strip_profiles.png, distinguishability_heatmap.png, feature_distributions.png,
auc_vs_labelling.png, examples_gap<break>.png, models.png, models_grid.png and examples.npz.

Run:   python fig2_filaments_actin.py
Score an existing feature table again, without simulating:
       python fig2_filaments_actin.py --features <folder>/features.csv
(the outputs are written next to the feature table)

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and condition; e.g. 5 for a quick
run), VLAB4MIC_N_PERMUTATIONS (200), VLAB4MIC_N_WORKERS (1; simulations and scoring in parallel
processes), VLAB4MIC_SAVE_IMAGES (0; 1 saves the images of every model and condition as a TIFF
stack, one realisation per frame, in <folder>/images).
"""

import copy
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.spatial.transform import Rotation
from tqdm.auto import tqdm

from vlab4mic import experiments
from vlab4mic.analysis import particle_measures
from vlab4mic.analysis.distinguishability import (
    distinguishability_from_features,
    plot_distinguishability_heatmaps,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path
from vlab4mic.utils.transform.points_transforms import rotate_pts_by_vector

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
N_PERMUTATIONS = int(os.environ.get("VLAB4MIC_N_PERMUTATIONS", 200))
N_WORKERS = int(os.environ.get("VLAB4MIC_N_WORKERS", 1))
SAVE_IMAGES = os.environ.get("VLAB4MIC_SAVE_IMAGES", "0") not in ("0", "", "false", "False")
RANDOM_SEED = 24
OUTPUT_NAME = "fig2_filaments_actin"  # base name of the output folder (after the date)

STRUCTURE = "6T25"
TARGET = "Phalloidin_EEP"
TARGET_VALUE = dict(atoms=["O2"], residues=["EEP"])  # 5-hydroxyl of 4,5-dihydroxyleucine
MIN_SUBUNIT_RESIDUES = 300  # chains with more residues than this are actin subunits
SAMPLE_DIMENSIONS = [1000, 1000, 10]  # nm
PARTICLE_POSITION = np.array([[0.5, 0.5, 0.0]])  # the filament at the centre of the sample
FILAMENT_LENGTH_NM = 500

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
EXAMPLE_WIDTH_NM = 160  # across it

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
    structure_global_normal_orientation="local_plane",  # normals from a local plane fit
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
    max_length_nm=FILAMENT_LENGTH_NM,
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


def linker_probe(labelling_efficiency=1):
    """Phalloidin: a linker of length 0 on the O2 atom of every EEP, labelled directly."""
    return dict(
        probe_template="Linker",
        probe_name=TARGET,
        probe_distance_to_epitope=0,  # the fluorophore on the target atom
        probe_target_type="Atom_residue",
        probe_target_value=copy.deepcopy(TARGET_VALUE),
        labelling_efficiency=labelling_efficiency,
    )


def labelled_segment():
    """The labelled segment (6T25 with its EEP targets), not imaged."""
    __, __, experiment = experiments.image_vsample(
        probe_list=[linker_probe()],
        **COMMON,
        run_simulation=False,
    )
    return experiment


def segment_points(experiment):
    """The parsed structure of the segment (Biopython model), its targets (nm) and the factor
    from its units to nm."""
    segment = experiment.structure
    to_nm = segment.scale / 1e-9
    targets = np.asarray(segment.get_target_coords(TARGET), dtype=float) * to_nm
    return segment.struct[0], targets, to_nm


def kabsch(p, q):
    """Rotation r and translation t taking points p onto q (least squares): q = p r^T + t."""
    pc, qc = p.mean(axis=0), q.mean(axis=0)
    u, __, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, qc - r @ pc


def subunit_screw(model, to_nm):
    """The screw between consecutive actin subunits: superposition of the CA atoms of every
    subunit onto the next, the subunits ordered along the axis of their centres. Returns its
    twist (degrees, about the axis), rise (nm, along the axis, positive), axis and a point of
    it; the actin chains in order, and the axis, twist and rise of every pair of consecutive
    subunits."""
    subunits = [chain for chain in model
                if sum(1 for residue in chain if residue.id[0] == " ") > MIN_SUBUNIT_RESIDUES]
    alpha = {chain.id: {residue.id[1]: residue["CA"].coord * to_nm for residue in chain
                        if residue.id[0] == " " and "CA" in residue} for chain in subunits}
    centres = np.array([np.mean(list(alpha[chain.id].values()), axis=0) for chain in subunits])
    __, eigenvectors = np.linalg.eigh(np.cov((centres - centres.mean(axis=0)).T))
    names = [subunits[i].id for i in np.argsort(centres @ eigenvectors[:, -1])]

    rotvecs, translations = [], []
    for a, b in zip(names, names[1:]):
        common_residues = sorted(set(alpha[a]) & set(alpha[b]))
        r, t = kabsch(np.array([alpha[a][i] for i in common_residues]),
                      np.array([alpha[b][i] for i in common_residues]))
        rotvecs.append(Rotation.from_matrix(r).as_rotvec())
        translations.append(t)
    axis = rotvecs[0] / np.linalg.norm(rotvecs[0])
    if np.mean([t @ axis for t in translations]) < 0:  # the axis points along the rise
        axis = -axis
    twists = np.array([np.degrees(rotvec @ axis) for rotvec in rotvecs])
    rises = np.array([t @ axis for t in translations])
    twist, rise = float(twists.mean()), float(rises.mean())
    r = Rotation.from_rotvec(np.radians(twist) * axis).as_matrix()
    t = np.mean(translations, axis=0)
    point = np.linalg.lstsq(np.eye(3) - r, t - (t @ axis) * axis, rcond=None)[0]
    point = point - (point @ axis) * axis  # the point of the axis closest to the origin
    return dict(twist_deg=twist, rise_nm=rise, axis=axis.tolist(), axis_point_nm=point.tolist(),
                subunits=names,
                pair_axes=[(rotvec / np.linalg.norm(rotvec)).tolist() for rotvec in rotvecs],
                pair_twists_deg=twists.tolist(), pair_rises_nm=rises.tolist())


def screw(points, ops, n):
    """Points moved by n steps of the screw: n twists about the axis, n rises along it."""
    axis, point = np.asarray(ops["axis"]), np.asarray(ops["axis_point_nm"])
    turn = Rotation.from_rotvec(np.radians(n * ops["twist_deg"]) * axis)
    return turn.apply(points - point) + point + n * ops["rise_nm"] * axis


def along_axis(points, ops):
    return (points - np.asarray(ops["axis_point_nm"])) @ np.asarray(ops["axis"])


def helix(first, ops, start, stop):
    """Sites of the helix of first, steps start ... stop - 1 of the screw."""
    return np.vstack([screw(first, ops, j) for j in range(start, stop)])


def extend_by_screw(targets, ops, length):
    """The filament: the EEP site of the first subunit along the axis moved by 0 ... N - 1
    steps of the screw, N = floor(length / rise) + 1 subunits."""
    n_subunits = int(np.floor(length / ops["rise_nm"] + 1e-9)) + 1
    start = int(np.argmin(along_axis(targets, ops)))
    return helix(targets[[start]], ops, 0, n_subunits)


def lay_flat(points, ops):
    """The filament with its axis along x, its middle at the origin."""
    axis = np.asarray(ops["axis"])
    middle = np.asarray(ops["axis_point_nm"]) + axis * np.mean(
        [along_axis(points, ops).min(), along_axis(points, ops).max()])
    return rotate_pts_by_vector(points - middle, axis, np.array([1.0, 0, 0]), np.zeros(3))


def particle_source(coordinates):
    """The filament as particle source data: its EEP sites (nm), centred at the reference point;
    no normals, a linker of length 0 labels the targets directly."""
    return dict(
        targets={TARGET: dict(coordinates=np.array(coordinates, dtype=float), normals=None)},
        scale=1e-9,
        reference_point=np.zeros(3),
        axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([1.0, 0.0, 0.0])),
    )


def load_filament(experiment, particle_data):
    """The filament as the single particle of the virtual sample of experiment, at the centre
    of the sample."""
    experiment.particle.source = dict()
    experiment.particle.load_source(**copy.deepcopy(particle_data))
    experiment.set_virtualsample_params(
        sample_dimensions=SAMPLE_DIMENSIONS,
        particle_positions=PARTICLE_POSITION,
        update_mode=False,
    )
    experiment.particle.generate_instance()
    experiment.build(modules=["coordinate_field", "imager"])
    return experiment


def filament_template():
    """The filament: the segment extended by the screw between its subunits over
    FILAMENT_LENGTH_NM and laid flat, as particle source data."""
    model, targets, to_nm = segment_points(labelled_segment())
    ops = subunit_screw(model, to_nm)
    coordinates = lay_flat(extend_by_screw(targets, ops, FILAMENT_LENGTH_NM), ops)
    print(f"filament: {len(coordinates)} EEP sites, {ops['twist_deg']:.2f} degrees and "
          f"{ops['rise_nm']:.3f} nm per subunit, laid flat along x")
    return particle_source(coordinates)


def position_along_filament(coordinates):
    """Position of each point along the principal axis, from the median (nm)."""
    centred = coordinates - coordinates.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov(centred.T))
    along = centred @ eigenvectors[:, -1]
    return along - np.median(along)


def break_filament(template, gap_nm):
    """Copy of the template without the EEP sites within gap_nm / 2 of its centre."""
    broken = copy.deepcopy(template)
    if gap_nm <= 0:
        return broken
    target = broken["targets"][TARGET]
    along = position_along_filament(target["coordinates"])
    span = np.ptp(along)
    assert 0.9 * FILAMENT_LENGTH_NM < span < 1.2 * FILAMENT_LENGTH_NM, \
        f"filament spans {span:.0f} along its axis: coordinates not in nm?"
    keep = np.abs(along) >= gap_nm / 2
    target["coordinates"] = target["coordinates"][keep]
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
    measures. Also returns the width metrics, masks, smoothed images, segment tables and strip
    profiles of every image."""
    tables = []
    masks = []
    filament_metrics_list = []
    rows = []
    blurred_list = []
    strips_list = []
    for result in replicates:
        mask, table, blurred = particle_measures.threshold_mask_and_measure(result["images"][modality][channel][0], method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = particle_measures.measure_filament_width(result["images"][modality][channel][0], pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = particle_measures.simple_width(result["images"][modality][channel][0])
        strips = particle_measures.filament_strip_profiles(result["images"][modality][channel][0], pixel_size=pixel_size, **(strip_kwargs or {}))
        strips_list.append(strips)
        filament_metrics_list.append(filament_metrics)
        blurred_list.append(blurred)
        tables.append(table)
        masks.append(mask)
        n_segments = table.shape[0]
        if n_segments == 0:
            min_eccentricity = 0
            max_eccentricity = 0
            mean_eccentricity = 0
            major_lenght = 0
            apparent_continuity = 0
        else:
            min_eccentricity = np.min(table["eccentricity"])
            max_eccentricity = np.max(table["eccentricity"])
            if len(list(table["eccentricity"])) > 1:
                mean_eccentricity = np.mean(table["eccentricity"])
            else:
                mean_eccentricity = table["eccentricity"][0]
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
    return df, filament_metrics_list, masks, blurred_list, tables, strips_list



def simulate_replicates(labelling_efficiency, seed, particle_data, n_replicates=N_REALISATIONS):
    """Independent realisations of a filament (run_replicates output), labelled with phalloidin
    at labelling_efficiency."""
    __, __, experiment = experiments.image_vsample(
        probe_list=[linker_probe(labelling_efficiency)],
        **dict(COMMON, random_seed=seed),
        run_simulation=False,
    )
    load_filament(experiment, particle_data)
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


def labelling_axis(ax):
    """Labelling efficiency on a log scale, ticks at the simulated values."""
    les = sorted(LABELLING_EFFICIENCIES)
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks(les, [f"{le:g}" for le in les], fontsize=8)
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
    """Epitopes of every model in the plane, from the centre of the intact filament (nm)."""
    centre = models[0][2]["targets"][TARGET]["coordinates"].mean(axis=0)
    return {epitope_key(model, gap):
            data["targets"][TARGET]["coordinates"][:, :2] - centre[:2]
            for model, gap, data in models}


def models_shown():
    """(model, break) of the intact filament and of every broken one."""
    return [("intact", 0)] + [("broken", gap) for gap in GAPS_NM]


def plot_models(examples, prefix):
    """The models, one per column: their epitopes from above (removed ones open) and an
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
        column[0].set_title(f"{model_label(model, gap)}\n{len(epitopes)} epitopes", fontsize=8)
        for ax, le in zip(column[1:], les):
            show(ax, examples[model_key(model, gap, le)])
    scale_bar(axes[0, 0], colour="k")
    axes[0, 0].set_ylabel("Epitopes", fontsize=8)
    for ax, le in zip(axes[1:, 0], les):
        scale_bar(ax)
        ax.set_ylabel(f"LE {le:g}", fontsize=8)
    fig.suptitle("Actin filament models: EEP sites (black, kept; open, removed) and one image at "
                 "every labelling efficiency; scale bar 100 nm", x=0.02, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(prefix + "models.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strip_examples(template, prefix, labelling_efficiency=1, gaps=(0, 50, 100)):
    """Mask and strip profiles of one image of the intact and broken filaments."""
    flat_fraction = FILAMENTS_PARAMS["strip_kwargs"]["flat_fraction"]
    fig, axes = plt.subplots(len(gaps), 4, figsize=(16, 3.2 * len(gaps)),
                             gridspec_kw=dict(width_ratios=[1, 1, 1.4, 2]))
    for row, gap in zip(axes, gaps):
        replicates = simulate_replicates(labelling_efficiency, RANDOM_SEED + gap,
                                         break_filament(template, gap), n_replicates=1)
        table, __, masks, __, __, strips = measure_filaments(replicates, **FILAMENTS_PARAMS)
        image, mask, strip = replicates[0]["images"]["SMLM"]["ch0"][0], masks[0], strips[0]
        pixel_size = FILAMENTS_PARAMS["pixel_size"]
        # image and mask around the filament, with the axis and the flat strips
        rows, cols = np.nonzero(image > 0.05 * image.max())
        crop = (slice(max(rows.min() - 30, 0), rows.max() + 30),
                slice(max(cols.min() - 30, 0), cols.max() + 30))
        origin = np.array([crop[0].start, crop[1].start])
        row[0].imshow(image[crop], cmap="grey", vmax=np.percentile(image[image > 0], 99.5))
        axis_points = (strip["centroid"] - origin
                       + np.outer(strip["positions"] / pixel_size, strip["direction"]))
        row[0].plot(axis_points[:, 1], axis_points[:, 0], c="C1", lw=0.8)
        row[0].scatter(axis_points[strip["flat"], 1], axis_points[strip["flat"], 0], c="red", s=6)
        row[0].set_title(f"gap {gap} nm, LE {labelling_efficiency}", loc="left", fontsize=9)
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
    fig.tight_layout()
    fig.savefig(prefix + "strip_profiles.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def feature_table(tables):
    """Intact and broken rows for every break size."""
    features = pd.concat(tables, ignore_index=True)
    intact = features[features.model == "intact"]
    broken = features[features.model == "broken"]
    # the intact images are shared by every break size
    return pd.concat([broken] + [intact.assign(gap_nm=gap) for gap in GAPS_NM],
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
                    found = measured.loc[replicas[k]]
                    detected = found.mask_number_of_segments >= MIN_SEGMENTS_BROKEN
                    ax.text(0.02, 0.95, f"{found.mask_number_of_segments:.0f} segment(s), "
                                        f"longest break {found.strip_longest_break:.0f} nm",
                            color="C2" if detected else "w", fontsize=5, ha="left", va="top",
                            transform=ax.transAxes)
                    if detected:
                        ax.spines[:].set_color("C2")
                        ax.spines[:].set_linewidth(2.5)
            if len(images):
                scale_bar(column[0])
            column[0].set_title(f"LE {le:g}\n{int(break_detected(data).sum())} / {len(data)} "
                                f"split", fontsize=9)
        fig.suptitle(f"Actin filament, break of {gap:g} nm: labelling efficiency per column; "
                     f"green, break detected ({MIN_SEGMENTS_BROKEN} or more mask segments); "
                     f"scale bar 100 nm", x=0.02, ha="left", fontsize=10)
        fig.tight_layout()
        fig.savefig(prefix + f"examples_gap{gap:g}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


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
    fig.suptitle("Actin filament: filaments compared in each cell of the distinguishability map "
                 "(intact above broken; AUC, all features); scale bar 100 nm",
                 x=0.02, ha="left", fontsize=10)
    fig.savefig(prefix + "models_grid.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_figure(examples, detection, scores, prefix):
    """The main figure: (a) one filament per break of FIGURE_EXAMPLE_GAPS_NM (rows; 0 nm: the intact
    filament) at every labelling efficiency, (b) the share of broken filaments whose break is
    detected, and the false alarm rate of intact filaments, (c) the AUC of intact against
    broken filaments from the number of mask segments alone, both against labelling
    efficiency."""
    les = sorted(LABELLING_EFFICIENCIES)
    breaks = [gap for gap in FIGURE_EXAMPLE_GAPS_NM if gap > 0]  # one line each in b and c
    colours = {gap: f"C{k}" for k, gap in enumerate(breaks)}
    colours[0] = "0.5"  # the intact filament, as in b
    fig = plt.figure(figsize=[11, 3.4 + 0.55 * len(FIGURE_EXAMPLE_GAPS_NM)])
    grid = fig.add_gridspec(2, 1, height_ratios=[0.25 * len(FIGURE_EXAMPLE_GAPS_NM) + 0.15, 1],
                            hspace=0.3)
    top = grid[0].subgridspec(len(FIGURE_EXAMPLE_GAPS_NM), len(les), wspace=0.05, hspace=0.08)
    for i, gap in enumerate(FIGURE_EXAMPLE_GAPS_NM):
        model = "broken" if gap > 0 else "intact"
        for k, le in enumerate(les):
            ax = fig.add_subplot(top[i, k])
            show(ax, examples[model_key(model, gap, le)])
            if i == 0:
                ax.set_title(f"LE {le:g}", fontsize=9)
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
    ax.set_title(f"Mask split into {MIN_SEGMENTS_BROKEN}+ segments "
                 f"(n = {detection.n_filaments.max()} per point)", fontsize=9)
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
        labelling_axis(axis)
    fig.savefig(prefix + "filament_figure.png", dpi=300, bbox_inches="tight")
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
        template = filament_template()
        plot_strip_examples(template, prefix)
        models = [("intact", 0, template)]
        for gap in GAPS_NM:
            broken = break_filament(template, gap)
            n = len(broken["targets"][TARGET]["coordinates"])
            print(f"gap {gap} nm: {n} of "
                  f"{len(template['targets'][TARGET]['coordinates'])} EEP sites")
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
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    plot_feature_distributions(features, prefix)
    features = features.fillna(FILL_VALUES)
    if examples:
        plot_examples(features, examples, prefix)
    scores = []
    for name, columns in FEATURE_SETS.items():
        print(f"features: {name}")
        table = distinguishability_from_features(
            features, "model", ("intact", "broken"), columns, by=CONDITIONS,
            controls=CONTROLS if name == "all" else None,
            replicate_column="replicate", n_permutations=N_PERMUTATIONS,
            random_state=RANDOM_SEED, n_workers=N_WORKERS,
        )
        table.insert(0, "feature_set", name)
        scores.append(table)
    scores = pd.concat(scores, ignore_index=True)
    scores.to_csv(prefix + "distinguishability.csv", index=False)

    fig = plot_distinguishability_heatmaps(scores, x="labelling_efficiency", y="gap_nm",
                                           facet="feature_set", facet_order=list(FEATURE_SETS))
    n_le = scores.labelling_efficiency.nunique()
    gaps = sorted(scores.gap_nm.unique())
    panels = fig.axes[:len(FEATURE_SETS)]
    for ax in panels:
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Break length (nm)")
    fig.suptitle("Actin filament: intact vs broken (ROC AUC)", x=0.02, ha="left")
    fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * len(gaps))
    fig.savefig(prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    plot_auc_vs_labelling(scores, prefix)
    detection = break_detection(features)
    detection.to_csv(prefix + "break_detection.csv", index=False)
    figure_values(detection, scores).to_csv(prefix + "filament_figure_values.csv", index=False)
    print(detection.round(3).to_string(index=False))
    if epitope_key("intact", 0) in examples:  # the epitopes of the models, from examples.npz
        plot_models(examples, prefix)
        plot_models_grid(examples, scores, prefix)
        plot_figure(examples, detection, scores, prefix)
    print(scores.pivot_table(index=["feature_set", "gap_nm"], columns="labelling_efficiency",
                             values="auc", sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
