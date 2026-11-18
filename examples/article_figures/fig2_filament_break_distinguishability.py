"""Fig. 2: can a break in a microtubule be detected, and at which labelling efficiency?

A 500 nm microtubule (5SYF placed along a line, as in fig2_filaments.py) is
labelled with a primary and an NHS secondary antibody and imaged by SMLM.
A second model has an artificial break in the middle: every epitope within
gap/2 of the centre of the filament is removed. Each model is simulated
with run_replicates at every labelling efficiency (of the secondary) and
every image is measured with measure_filaments (widths, mask continuity,
number of segments, eccentricities...), giving a feature table with one row
per image, as the "merged" table of fig2_filaments.py.

Continuity is measured in two ways:
- mask: the image is smoothed (Gaussian, 4 px = 8 nm) and thresholded
  (triangle method), and the connected segments are measured. With these
  settings an intact filament is a single segment at labelling efficiency
  1, and a 100 nm break splits it in two.
- strip profiles (particle_measures.filament_strip_profiles): the image is
  divided into strips of 10 nm along the central axis of the filament and
  each strip is projected into a line profile across the axis. A flat
  profile (amplitude below 10% of the reference) is a strip that no
  filament crosses; flat strips between strips that are not flat are an
  apparent break.

The table is scored per break size and labelling efficiency with
distinguishability_from_features: out-of-fold ROC AUC of a classifier
(intact vs broken) with a 95% bootstrap interval, a null control (two
random halves of the intact filaments, expected 0.5) and a permutation
p-value, for each feature set (mask, strip profiles, all). Each metric is
also scored alone (<metric>_auc, with all features). A gap of 0 nm is
a second, independently simulated intact filament: a negative control
(expected 0.5).

Outputs, in the VLab4Mic output directory (<date>fig2_filament_break_...):
features.csv (one row per image), distinguishability.csv (one row per
break size, labelling efficiency and feature set), settings.json,
models.png, strip_profiles.png and distinguishability_heatmap.png.

Run:   python fig2_filament_break_distinguishability.py
Score an existing feature table again, without simulating:
       python fig2_filament_break_distinguishability.py --features <...features.csv>

Environment variables: VLAB4MIC_N_REALISATIONS (100 per model and
condition; e.g. 5 for a quick run), VLAB4MIC_N_PERMUTATIONS (200),
VLAB4MIC_N_WORKERS (1; simulations and scoring in parallel processes).
"""

import copy
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
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
RANDOM_SEED = 24

STRUCTURE = "5SYF"
PROBE_MODEL = "1IGT"  # used for every primary and secondary antibody
SAMPLE_DIMENSIONS = [1000, 1000, 10]  # nm
FILAMENT_LENGTH_NM = 500
SPACING_NM = 7  # distance between consecutive 5SYF along the filament
N_EPITOPES = 1872  # epitopes of the intact filament

LABELLING_EFFICIENCIES = [0.01, 0.05, 0.2, 0.5, 1]
GAPS_NM = [0, 25, 50, 100, 200]  # length of the break; 0: intact (negative control)
REFERENCE_GAP_NM = 50

# measure_filaments columns used as features (the others are settings)
MASK_FEATURES = ["mask_continuity", "mask_major_lenght", "mask_number_of_segments",
                 "min_eccentricity", "max_eccentricity", "mean_eccentricity"]
STRIP_FEATURES = ["strip_fraction_flat", "strip_n_breaks", "strip_longest_break",
                  "strip_min_amplitude"]
WIDTHS = ["mean_filament_width", "std_filament_width", "median_filament_width", "s_width"]
FEATURES = WIDTHS + MASK_FEATURES + STRIP_FEATURES
FEATURE_SETS = {"mask": MASK_FEATURES, "strip profiles": STRIP_FEATURES, "all": FEATURES}
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
    smooth_sigma=4,  # px; with triangle: one segment for an intact filament at LE 1
    threshold_method="triangle",
    strip_kwargs=dict(
        strip_length=10,  # nm, along the filament axis
        half_width=50,  # nm, profiles across the axis
        background_width=10,  # nm at each end of a profile
        flat_fraction=0.1,  # flat: amplitude below this fraction of the reference
    ),
)


def line_of_points2(
    x_dimension,
    y_dimension,
    spacing,
    npoints=None,
    line_position=(0.5, 0.5),
    orientation=0,
    length=None,
):
    if (npoints is None) == (length is None):
        raise ValueError("Specify exactly one of npoints or length.")
    if length is not None:
        npoints = int(np.floor(length / spacing + 1e-9)) + 1
    angle = np.deg2rad(orientation)
    # signed distances from the centre, symmetric on both sides
    distances = (np.arange(npoints) - (npoints - 1) / 2) * spacing
    points = np.zeros((npoints, 3))
    points[:, 0] = line_position[0] + distances * np.cos(angle) / x_dimension
    points[:, 1] = line_position[1] + distances * np.sin(angle) / y_dimension
    return points


def filament_template():
    """Epitopes and normals of the intact filament, as particle source data."""
    primary = dict(
        probe_template="Antibody",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope=0,
        probe_model=PROBE_MODEL,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        labelling_efficiency=1,
    )
    pts = line_of_points2(SAMPLE_DIMENSIONS[0], SAMPLE_DIMENSIONS[1], spacing=SPACING_NM,
                          line_position=(0.5, 0.5), orientation=0, length=FILAMENT_LENGTH_NM)
    __, __, experiment = experiments.image_vsample(
        probe_list=[primary],
        **COMMON,
        random_orientations=True,
        xz_orientations=[90, ],
    )
    experiment.clear_virtual_sample()
    experiment.set_virtualsample_params(particle_positions=pts)
    experiment.build(modules=["coordinate_field", "imager"])
    targets = [molecule.source["targets"]["ATubulin_epitopes"]
               for molecule in experiment.coordinate_field.molecules]
    coordinates = np.vstack([target["coordinates"] for target in targets])
    normals = np.vstack([target["normals"] for target in targets])
    assert coordinates.shape == (N_EPITOPES, 3)
    return dict(
        targets={"ATubulin_epitopes": dict(coordinates=coordinates, normals=normals)},
        scale=1e-9,
        reference_point=experiment.coordinate_field.get_field_param("absolute_reference_point").copy(),
        axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0])),
    )


def position_along_filament(coordinates):
    """Position of each point along the principal axis, from the median (nm)."""
    centred = coordinates - coordinates.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov(centred.T))
    along = centred @ eigenvectors[:, -1]
    return along - np.median(along)


def break_filament(template, gap_nm):
    """Copy of the template without the epitopes within gap_nm / 2 of its centre."""
    broken = copy.deepcopy(template)
    if gap_nm <= 0:
        return broken
    target = broken["targets"]["ATubulin_epitopes"]
    along = position_along_filament(target["coordinates"])
    span = np.ptp(along)
    # the epitopes extend a few nm beyond the line of 5SYF positions
    assert 0.9 * FILAMENT_LENGTH_NM < span < 1.2 * FILAMENT_LENGTH_NM, \
        f"filament spans {span:.0f} along its axis: coordinates not in nm?"
    keep = np.abs(along) >= gap_nm / 2
    target["coordinates"] = target["coordinates"][keep]
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
    tables = []
    masks = []
    eccentricities = []
    major_axis_lengths = []
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
            major_axis_lengths.append(major_lenght)
            apparent_continuity = 0
        else:
            min_eccentricity = np.min(table["eccentricity"])
            max_eccentricity = np.max(table["eccentricity"])
            if len(list(table["eccentricity"])) > 1:
                mean_eccentricity = np.mean(table["eccentricity"])
            else:
                mean_eccentricity = table["eccentricity"][0]
            major_lenght = table.iloc[table['major_axis_length'].idxmax()]["major_axis_length"]
            major_axis_lengths.append(major_lenght)
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
    """Independent realisations of a filament (run_replicates output)."""
    primary = dict(
        probe_template="Antibody",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope=0,
        probe_model=PROBE_MODEL,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        labelling_efficiency=1,
        probe_paratope="TITCH",
        probe_secondary_epitope="EKKNW",
    )
    secondary = dict(
        probe_template="Antibody",
        probe_name="Secondary_L1",
        probe_model=PROBE_MODEL,
        probe_target_type="Primary",
        probe_target_value="ATubulin_epitopes",
        probe_conjugation_target_info=dict(
            type="Atom_residue",
            value=dict(atoms=["CA"], residues=["LYS"]),
        ),
        probe_paratope="TITCH",
        labelling_efficiency=labelling_efficiency,
    )
    __, __, experiment = experiments.image_vsample(
        primary_probe=primary,
        secondary_probe=secondary,
        **dict(COMMON, random_seed=seed),
    )
    experiment.particle.source = dict()
    experiment.particle.load_source(**copy.deepcopy(particle_data))
    experiment.set_virtualsample_params(
        sample_dimensions=SAMPLE_DIMENSIONS,
        update_mode=False,
    )
    experiment.particle.generate_instance()
    experiment.build(modules=["coordinate_field", "imager"])
    return experiment.run_replicates(n_replicates=n_replicates, modality="SMLM")


def simulate(job):
    """One model at one labelling efficiency: its feature table (one row per image)."""
    model, gap_nm, labelling_efficiency, seed, particle_data = job
    replicates = simulate_replicates(labelling_efficiency, seed, particle_data)
    table = measure_filaments(replicates=replicates, **FILAMENTS_PARAMS)[0]
    return table.assign(model=model, gap_nm=gap_nm, labelling_efficiency=labelling_efficiency,
                        seed=seed)


def plot_models(template, prefix):
    """Epitopes of the intact filament and of the reference break, from above."""
    models = {"intact": template,
              f"broken ({REFERENCE_GAP_NM} nm)": break_filament(template, REFERENCE_GAP_NM)}
    fig, axes = plt.subplots(len(models), 1, figsize=(10, 3), sharex=True, sharey=True)
    for ax, (name, data) in zip(axes, models.items()):
        coordinates = data["targets"]["ATubulin_epitopes"]["coordinates"]
        ax.scatter(coordinates[:, 0], coordinates[:, 1], s=1, c="k")
        ax.set_aspect("equal")
        ax.set_title(f"{name}: {len(coordinates)} epitopes", loc="left", fontsize=9)
    axes[-1].set_xlabel("x (nm)")
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


def main():
    if "--features" in sys.argv:
        features_path = sys.argv[sys.argv.index("--features") + 1]
        features = pd.read_csv(features_path)
        prefix = features_path.replace("features.csv", "")
    else:
        prefix = os.path.join(str(output_path),
                              ExperimentParametrisation().date_as_string + "fig2_filament_break_")
        template = filament_template()
        plot_models(template, prefix)
        plot_strip_examples(template, prefix)
        models = [("intact", 0, template)]
        for gap in GAPS_NM:
            broken = break_filament(template, gap)
            n = len(broken["targets"]["ATubulin_epitopes"]["coordinates"])
            print(f"gap {gap} nm: {n} of {N_EPITOPES} epitopes")
            models.append(("broken", gap, broken))
        blocks = [(model, gap, le) for model, gap, __ in models for le in LABELLING_EFFICIENCIES]
        data = {(model, gap): particle_data for model, gap, particle_data in models}
        seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [(model, gap, le, int(seed), data[(model, gap)])
                for (model, gap, le), seed in zip(blocks, seeds)]
        print(f"{len(jobs)} simulations, {N_REALISATIONS} realisations each, "
              f"{N_WORKERS} worker(s)")
        if N_WORKERS > 1:
            with ProcessPoolExecutor(N_WORKERS) as pool:
                tables = list(tqdm(pool.map(simulate, jobs), total=len(jobs), desc="simulations"))
        else:
            tables = [simulate(job) for job in tqdm(jobs, desc="simulations")]
        features = feature_table(tables)
        features.to_csv(prefix + "features.csv", index=False)
        settings = {k: v for k, v in globals().items()
                    if k.isupper() and isinstance(v, (int, float, str, list, dict))}
        with open(prefix + "settings.json", "w") as handle:
            json.dump(settings, handle, indent=1, default=str)

    features = features.fillna(FILL_VALUES)
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
        if REFERENCE_GAP_NM in gaps:
            # outline the row of the reference break
            ax.add_patch(Rectangle((-0.5, gaps.index(REFERENCE_GAP_NM) - 0.5), n_le, 1,
                                   fill=False, linewidth=2, edgecolor="red"))
        ax.set_xlabel("Labelling efficiency")
    panels[0].set_ylabel("Break length (nm)")
    fig.suptitle("Intact vs broken filament (ROC AUC)", x=0.02, ha="left")
    fig.set_size_inches(1.0 + len(panels) * (0.8 + 0.45 * n_le), 1.5 + 0.5 * len(gaps))
    fig.savefig(prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(scores.pivot_table(index=["feature_set", "gap_nm"], columns="labelling_efficiency",
                             values="auc", sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
