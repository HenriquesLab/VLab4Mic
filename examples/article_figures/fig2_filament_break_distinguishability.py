"""Fig. 2: can a break in a microtubule be detected, and at which labelling efficiency?

A 500 nm microtubule (5SYF placed along a line, as in fig2_filaments.py) is
labelled with a primary and an NHS secondary antibody and imaged by SMLM.
A second model has an artificial break in the middle: every epitope within
gap/2 of the centre of the filament is removed. Each model is simulated
with run_replicates at every labelling efficiency (of the secondary) and
every image is measured with measure_filaments (widths, mask continuity,
number of segments, eccentricities...), giving a feature table with one row
per image, as the "merged" table of fig2_filaments.py.

The table is scored per break size and labelling efficiency with
distinguishability_from_features: out-of-fold ROC AUC of a classifier
(intact vs broken) with a 95% bootstrap interval, a null control (two
random halves of the intact filaments, expected 0.5) and a permutation
p-value. Each metric is also scored alone (<metric>_auc). A gap of 0 nm is
a second, independently simulated intact filament: a negative control
(expected 0.5).

Outputs, in the VLab4Mic output directory (<date>fig2_filament_break_...):
features.csv (one row per image), distinguishability.csv (one row per
break size and labelling efficiency), settings.json, models.png and
distinguishability_heatmap.png.

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

LABELLING_EFFICIENCIES = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
GAPS_NM = [0, 10, 25, 50, 100]  # length of the break; 0: intact (negative control)
REFERENCE_GAP_NM = 50

# measure_filaments columns used as features (the others are settings)
FEATURES = ["mean_filament_width", "std_filament_width", "median_filament_width",
            "mask_continuity", "mask_major_lenght", "mask_number_of_segments", "s_width",
            "min_eccentricity", "max_eccentricity", "mean_eccentricity"]
CONTROLS = {feature: feature for feature in FEATURES}  # each alone
WIDTHS = ["mean_filament_width", "std_filament_width", "median_filament_width", "s_width"]
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
    smooth_sigma=0.0,
    threshold_method="otsu",
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
                      threshold_method_kwargs=None,):
    tables = []
    masks = []
    eccentricities = []
    major_axis_lengths = []
    filament_metrics_list = []
    rows = []
    blurred_list = []
    for result in replicates:
        mask, table, blurred = particle_measures.threshold_mask_and_measure(result["images"][modality][channel][0], method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = particle_measures.measure_filament_width(result["images"][modality][channel][0], pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = particle_measures.simple_width(result["images"][modality][channel][0])
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
            mean_eccentricity=mean_eccentricity
        ))
    df = pd.DataFrame(rows)
    return df, filament_metrics_list, masks, blurred_list, tables


def simulate(job):
    """One model at one labelling efficiency: its feature table (one row per image)."""
    model, gap_nm, labelling_efficiency, seed, particle_data = job
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
    replicates = experiment.run_replicates(n_replicates=N_REALISATIONS, modality="SMLM")
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

    # no measurable width (e.g. too few labels): width 0
    features[WIDTHS] = features[WIDTHS].fillna(0)
    scores = distinguishability_from_features(
        features, "model", ("intact", "broken"), FEATURES, by=CONDITIONS, controls=CONTROLS,
        replicate_column="replicate", n_permutations=N_PERMUTATIONS,
        random_state=RANDOM_SEED, n_workers=N_WORKERS,
    )
    scores.to_csv(prefix + "distinguishability.csv", index=False)
    fig = plot_distinguishability_heatmaps(scores, x="labelling_efficiency", y="gap_nm")
    ax = fig.axes[0]
    gaps = sorted(scores.gap_nm.unique())
    if REFERENCE_GAP_NM in gaps:
        # outline the row of the reference break
        ax.add_patch(Rectangle((-0.5, gaps.index(REFERENCE_GAP_NM) - 0.5),
                               scores.labelling_efficiency.nunique(), 1,
                               fill=False, linewidth=2, edgecolor="red"))
    ax.set_xlabel("Labelling efficiency")
    ax.set_ylabel("Break length (nm)")
    ax.set_title("Intact vs broken filament", loc="left")
    fig.set_size_inches(1.0 + 0.6 * len(LABELLING_EFFICIENCIES), 1.0 + 0.6 * len(GAPS_NM))
    fig.savefig(prefix + "distinguishability_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(scores[CONDITIONS + ["auc", "auc_low", "auc_high", "null_auc", "p_value"]].to_string())


if __name__ == "__main__":
    main()
