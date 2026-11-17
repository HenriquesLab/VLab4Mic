"""Supp. Fig. S4: corners detected on PCNA across labelling efficiencies.

PCNA (4ZTD) labelled at SER-186 of each chain (site-specific, no linkage),
about 6 nm apart, i.e. three sites on a ring of
~3.3 nm radius. Same corner counting as fig2_npc_corner_detection.py
with 3 segments instead of 8.
"""

import numpy as np
import matplotlib.pyplot as plt
from vlab4mic import experiments
from vlab4mic.analysis.particle_measures import (
    plot_ring_corners,
    ring_corners_from_replicates,
)
import pandas as pd
import os
import seaborn as sns


###########

RANDOM_SEED = 61
PRECISION_NM = 1
N_SEGMENTS = 3  # PCNA trimer: three labelled sites
MEAN_RADIUS_PCNA_NM = 3.3  # SER-186 sites, ~6 nm apart
R_MIN_NM = 1.5
R_MAX_NM = 5
PIXEL_NM = 0.5
nlocalisations = 50
RENDERING_KERNEL_NM = 1.0

labelling_efficiencies = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
replicates_per_le = dict()
parameters_le = dict()
n_replicates_le = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
import copy

for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="Linker",
        probe_name="pcna_ser186",
        probe_target_type="Atom_residue",
        probe_target_value={"atoms": ["CA"], "residues": ["SER"], "position": 186},
        probe_distance_to_epitope=0,
        labelling_efficiency=labeff,
    )
    _, _, my_experiment_le = experiments.image_vsample(
        structure="4ZTD",
        **primary_le,
        random_rotations=True,
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[100, 100, 50],
        sample_inital_orientation=[0, 0, 1],  # trimer seen from above
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED,
    )
    my_experiment_le.update_modality(
        modality_name="SMLM",
        pixelsize_nm=PIXEL_NM,
        psf_voxel_nm=PIXEL_NM,
        simulate_localistations=True,
        nlocalisations=nlocalisations,
        lateral_precision=PRECISION_NM,
        axial_precision=PRECISION_NM,
        rendering_kernel_nm=RENDERING_KERNEL_NM,
    )

    replicates_le = my_experiment_le.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le[labeff] = copy.copy(replicates_le)
    parameters_le[labeff] = my_experiment_le.get_parameters()



#######
counting_measures_per_le = dict()
split_z = False
min_locs_per_segment = 4
sectorcounting_params_le = dict(
    modality="SMLM",
    channel="ch0",
    split_z = split_z,
    n_corners=N_SEGMENTS,
    radius=MEAN_RADIUS_PCNA_NM,
    r_min=R_MIN_NM,
    r_max=R_MAX_NM,
    min_per_section=min_locs_per_segment,
    min_locs_to_circlefit=10,
    max_radius=R_MAX_NM,
    resolve_spots=True,  # corners must be separate spots, not one cloud
)

for i in labelling_efficiencies:
    results, images_list, rows = ring_corners_from_replicates(
    replicates=replicates_per_le[i],
    **sectorcounting_params_le
    )
    counting_measures_per_le[i] = dict()
    counting_measures_per_le[i]["dataframe"] = results
    counting_measures_per_le[i]["images"] = images_list
    counting_measures_per_le[i]["rows"] = rows

dfs = dict()
for lef in labelling_efficiencies:
    dfs[lef] = counting_measures_per_le[lef]["dataframe"]

merged = pd.concat(
    [d.assign(source=name) for name, d in dfs.items()],
    ignore_index=True,
)
######## plot 1
fig, ax = plt.subplots(figsize=(5, 3.5))
sns.boxplot(data=merged, x="source", y="n_corners", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Number of Corners detected")
plt.title("Corners detected on PCNA simuations" + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_PCNA_simuations.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

#### plot 2 and 3
for lab_eff_reference in [0.3, 0.4, 0.7, 1]:
    counts = counting_measures_per_le[lab_eff_reference]["dataframe"].n_corners
    hist =  np.bincount(counts, minlength=N_SEGMENTS + 1)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(np.arange(N_SEGMENTS + 1), hist, label="detected")
    ax.set_xlabel("corners per PCNA")
    ax.set_ylabel("PCNAs")
    ax.set_xticks(np.arange(N_SEGMENTS + 1))
    ax.legend()
    ax.set_title("Corners detected on PCNA simuations. LE = " + str(lab_eff_reference) + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
    filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_PCNA_simuations_' + str(lab_eff_reference) + 'LE.png'
    filename2 = os.path.join(my_experiment_le.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()



### plot 4

n_le = len(labelling_efficiencies)
nrows = 5
fig, ax = plt.subplots(nrows, n_le, figsize=(30,10))
for row in range(nrows):
    for i, lef in enumerate(labelling_efficiencies):
        ax[row,i].imshow(counting_measures_per_le[lef]["images"][row][0], cmap = "grey")
        ax[row,i].set_title(f"LE = {lef}")
        ax[row,i].set_xticks([])
        ax[row,i].set_yticks([])
filename = my_experiment_le.date_as_string + 'vlab4mic_PCNAs_across_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.suptitle("Precision: " + str(PRECISION_NM) + " nm")
plt.tight_layout()
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()


## plot 5
# one example particle per number of detected corners
threshold_for_plot = [1, 2, 3]
lab_eff_reference = 0.3
_df = counting_measures_per_le[lab_eff_reference]["dataframe"]
for ncorners in threshold_for_plot:
    ids = list(_df.loc[_df["n_corners"] == ncorners].replicate)
    if not ids:
        continue
    fig = plot_ring_corners(counting_measures_per_le[lab_eff_reference]["rows"][ids[0]])
    filename = my_experiment_le.date_as_string + 'vlab4mic_pcna_corners_example_' + str(ncorners) + '.png'
    filename2 = os.path.join(my_experiment_le.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()

# plot 6
fig, ax = plt.subplots(figsize=(5, 3.5))
sns.boxplot(data=merged, x="source", y="radius_fit", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Radius fit")
filename = my_experiment_le.date_as_string + 'vlab4mic_pcna_radius_fit_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()


########### random orientations
# Same simulations and corner counting, but every replicate gets a random 3D
# orientation instead of the trimer seen from above. Tilted rings project
# as ellipses, so corners and radius are compared with the top-view case.

replicates_per_le_ro = dict()
parameters_le_ro = dict()
for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="Linker",
        probe_name="pcna_ser186",
        probe_target_type="Atom_residue",
        probe_target_value={"atoms": ["CA"], "residues": ["SER"], "position": 186},
        probe_distance_to_epitope=0,
        labelling_efficiency=labeff,
    )
    _, _, my_experiment_ro = experiments.image_vsample(
        structure="4ZTD",
        **primary_le,
        random_orientations=True,
        random_rotations=True,
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[100, 100, 50],
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED,
    )
    my_experiment_ro.update_modality(
        modality_name="SMLM",
        pixelsize_nm=PIXEL_NM,
        psf_voxel_nm=PIXEL_NM,
        simulate_localistations=True,
        nlocalisations=nlocalisations,
        lateral_precision=PRECISION_NM,
        axial_precision=PRECISION_NM,
        rendering_kernel_nm=RENDERING_KERNEL_NM,
    )

    replicates_ro = my_experiment_ro.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le_ro[labeff] = copy.copy(replicates_ro)
    parameters_le_ro[labeff] = my_experiment_ro.get_parameters()

counting_measures_per_le_ro = dict()
for i in labelling_efficiencies:
    results, images_list, rows = ring_corners_from_replicates(
    replicates=replicates_per_le_ro[i],
    **sectorcounting_params_le
    )
    counting_measures_per_le_ro[i] = dict()
    counting_measures_per_le_ro[i]["dataframe"] = results
    counting_measures_per_le_ro[i]["images"] = images_list
    counting_measures_per_le_ro[i]["rows"] = rows

merged_ro = pd.concat(
    [counting_measures_per_le_ro[lef]["dataframe"].assign(source=lef) for lef in labelling_efficiencies],
    ignore_index=True,
)
merged_orientations = pd.concat(
    [merged.assign(orientation="top view"), merged_ro.assign(orientation="random")],
    ignore_index=True,
)

print(
    merged_orientations.groupby(["orientation", "source"])[["n_corners", "radius_fit"]]
    .mean()
    .unstack("orientation")
)

# plot 7: corners vs LE, top view vs random orientations
fig, ax = plt.subplots(figsize=(7, 3.5))
sns.boxplot(data=merged_orientations, x="source", y="n_corners", hue="orientation", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Number of Corners detected")
plt.title("Corners detected on PCNA simuations, top view vs random orientations" + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
filename = my_experiment_ro.date_as_string + 'vlab4mic_corners_detected_on_PCNA_simuations_random_orientations.png'
filename2 = os.path.join(my_experiment_ro.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 8: corner histograms at reference LEs
width = 0.4
for lab_eff_reference in [0.3, 0.4, 0.7, 1]:
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for k, (label, measures) in enumerate(
        [("top view", counting_measures_per_le), ("random", counting_measures_per_le_ro)]
    ):
        counts = measures[lab_eff_reference]["dataframe"].n_corners
        hist = np.bincount(counts, minlength=N_SEGMENTS + 1)
        ax.bar(np.arange(N_SEGMENTS + 1) + (k - 0.5) * width, hist, width=width, label=label)
    ax.set_xlabel("corners per PCNA")
    ax.set_ylabel("PCNAs")
    ax.set_xticks(np.arange(N_SEGMENTS + 1))
    ax.legend()
    ax.set_title("Corners detected on PCNA simuations. LE = " + str(lab_eff_reference) + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
    filename = my_experiment_ro.date_as_string + 'vlab4mic_corners_detected_on_PCNA_simuations_random_orientations_' + str(lab_eff_reference) + 'LE.png'
    filename2 = os.path.join(my_experiment_ro.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()

# plot 9: radius fit vs LE, top view vs random orientations
fig, ax = plt.subplots(figsize=(7, 3.5))
sns.boxplot(data=merged_orientations, x="source", y="radius_fit", hue="orientation", ax=ax)
ax.axhline(MEAN_RADIUS_PCNA_NM, color="k", ls=":", lw=0.8)
plt.xlabel("Labelling efficiency")
plt.ylabel("Radius fit")
filename = my_experiment_ro.date_as_string + 'vlab4mic_pcna_radius_fit_LE_random_orientations.png'
filename2 = os.path.join(my_experiment_ro.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 10: example images with random orientations
fig, ax = plt.subplots(nrows, n_le, figsize=(30,10))
for row in range(nrows):
    for i, lef in enumerate(labelling_efficiencies):
        ax[row,i].imshow(counting_measures_per_le_ro[lef]["images"][row][0], cmap = "grey")
        ax[row,i].set_title(f"LE = {lef}")
        ax[row,i].set_xticks([])
        ax[row,i].set_yticks([])
filename = my_experiment_ro.date_as_string + 'vlab4mic_PCNAs_across_LE_random_orientations.png'
filename2 = os.path.join(my_experiment_ro.output_directory, filename)
fig.suptitle("Random orientations. Precision: " + str(PRECISION_NM) + " nm")
plt.tight_layout()
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 11: one example particle per number of detected corners, random orientations
lab_eff_reference = 0.3
_df = counting_measures_per_le_ro[lab_eff_reference]["dataframe"]
for ncorners in threshold_for_plot:
    ids = list(_df.loc[_df["n_corners"] == ncorners].replicate)
    if not ids:
        continue
    fig = plot_ring_corners(counting_measures_per_le_ro[lab_eff_reference]["rows"][ids[0]])
    filename = my_experiment_ro.date_as_string + 'vlab4mic_pcna_corners_example_random_orientations_' + str(ncorners) + '.png'
    filename2 = os.path.join(my_experiment_ro.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()


########### SMLM modalities
# Top-view PCNA (random in-plane rotations, no random orientations) imaged
# with the localisation precision reported for two SMLM modalities. Only the localisation precision changes; corners are counted on
# the localisations, so the rendering kernel stays at RENDERING_KERNEL_NM.

PUBLISHED = {"DNA-PAINT": 3.0, "dSTORM": 6.5}
probe = dict(
    probe_template="Linker",
    probe_name="pcna_ser186",
    probe_target_type="Atom_residue",
    probe_target_value={"atoms": ["CA"], "residues": ["SER"], "position": 186},
    probe_distance_to_epitope=0,
)

counting_measures_per_modality = dict()
for modality_label, precision in PUBLISHED.items():
    counting_measures_per_modality[modality_label] = dict()
    for labeff in labelling_efficiencies:
        _, _, my_experiment_mod = experiments.image_vsample(
            structure="4ZTD",
            **probe,
            labelling_efficiency=labeff,
            random_orientations=False,  # always seen from above
            random_rotations=True,  # random in-plane rotation only
            multimodal=["SMLM"],
            number_of_particles=1,
            sample_dimensions=[100, 100, 50],
            sample_inital_orientation=[0, 0, 1],  # trimer seen from above
            run_simulation=False,
            clear_experiment=True,
            random_seed=RANDOM_SEED,
        )
        my_experiment_mod.update_modality(
            modality_name="SMLM",
            pixelsize_nm=PIXEL_NM,
            psf_voxel_nm=PIXEL_NM,
            simulate_localistations=True,
            nlocalisations=nlocalisations,
            lateral_precision=precision,
            axial_precision=precision,
            rendering_kernel_nm=RENDERING_KERNEL_NM,
        )
        replicates_mod = my_experiment_mod.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
        results, images_list, rows = ring_corners_from_replicates(
            replicates=replicates_mod,
            **sectorcounting_params_le
        )
        counting_measures_per_modality[modality_label][labeff] = dict(
            dataframe=results, images=images_list, rows=rows
        )

merged_modalities = pd.concat(
    [
        counting_measures_per_modality[m][lef]["dataframe"].assign(
            modality=f"{m} ({PUBLISHED[m]} nm)", source=lef
        )
        for m in PUBLISHED
        for lef in labelling_efficiencies
    ],
    ignore_index=True,
)
correct_counts = (
    merged_modalities.assign(correct=merged_modalities["n_corners"] == N_SEGMENTS)
    .groupby(["modality", "source"])["correct"]
    .sum()
    .unstack("modality")
)
print(correct_counts)
# same count from occupied segments only, i.e. without requiring separate spots
print(
    merged_modalities.assign(correct=merged_modalities["n_segments_occupied"] == N_SEGMENTS)
    .groupby(["modality", "source"])["correct"]
    .sum()
    .unstack("modality")
)

# plot 12: PCNAs detected with 3 corners vs LE, per SMLM modality
fig, ax = plt.subplots(figsize=(5, 3.5))
for modality_label in correct_counts.columns:
    ax.plot(correct_counts.index, correct_counts[modality_label], marker="o", lw=2, label=modality_label)
# PCNAs with all three sites labelled: the best any modality can do
ax.plot(
    labelling_efficiencies,
    n_replicates_le * np.array(labelling_efficiencies) ** N_SEGMENTS,
    color="grey", ls="--", lw=1, label="all 3 sites labelled (LE³)",
)
ax.set_xlabel("Labelling efficiency")
ax.set_ylabel(f"PCNAs with {N_SEGMENTS} corners detected")
ax.set_ylim(0, n_replicates_le)
ax.set_xticks(labelling_efficiencies)
ax.legend()
ax.set_title(
    f"PCNA correctly detected (n = {n_replicates_le} per LE)" + "\n"
    + "Min locs per sector = " + str(min_locs_per_segment)
)
filename = my_experiment_mod.date_as_string + 'vlab4mic_pcna_3corners_detected_SMLM_modalities.png'
filename2 = os.path.join(my_experiment_mod.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 13: example images per modality across LE
fig, ax = plt.subplots(len(PUBLISHED), n_le, figsize=(30, 6))
for row, modality_label in enumerate(PUBLISHED):
    for i, lef in enumerate(labelling_efficiencies):
        ax[row, i].imshow(counting_measures_per_modality[modality_label][lef]["images"][0][0], cmap="grey")
        ax[row, i].set_title(f"{modality_label}, LE = {lef}")
        ax[row, i].set_xticks([])
        ax[row, i].set_yticks([])
filename = my_experiment_mod.date_as_string + 'vlab4mic_PCNAs_across_LE_SMLM_modalities.png'
filename2 = os.path.join(my_experiment_mod.output_directory, filename)
plt.tight_layout()
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 14: gallery of rendered PCNAs per modality at reference LEs.
# DNA-PAINT particles where three separate spots were detected are framed.
GALLERY_ROWS, GALLERY_COLS = 5, 10
HIGHLIGHT_MODALITY = "DNA-PAINT"
HIGHLIGHT_COLOR = "tab:orange"
for lab_eff_reference in [0.7, 1]:
    fig = plt.figure(figsize=(2 * GALLERY_COLS, 2.1 * GALLERY_ROWS * len(PUBLISHED)))
    subfigs = fig.subfigures(len(PUBLISHED), 1, hspace=0.05)
    for subfig, modality_label in zip(subfigs, PUBLISHED):
        measures = counting_measures_per_modality[modality_label][lab_eff_reference]
        n_corners = measures["dataframe"]["n_corners"].to_numpy()
        n_shown = min(GALLERY_ROWS * GALLERY_COLS, len(measures["images"]))
        highlight = modality_label == HIGHLIGHT_MODALITY
        ax = subfig.subplots(GALLERY_ROWS, GALLERY_COLS)
        for k, a in enumerate(ax.flat):
            a.set_xticks([])
            a.set_yticks([])
            if k >= n_shown:
                a.axis("off")
                continue
            a.imshow(measures["images"][k][0], cmap="grey")
            if highlight and n_corners[k] == N_SEGMENTS:
                for spine in a.spines.values():
                    spine.set_color(HIGHLIGHT_COLOR)
                    spine.set_linewidth(4)
        title = f"{modality_label} ({PUBLISHED[modality_label]} nm precision)"
        if highlight:
            n_resolved = int((n_corners[:n_shown] == N_SEGMENTS).sum())
            title += f": {n_resolved}/{n_shown} with {N_SEGMENTS} separate spots (framed)"
        subfig.suptitle(title, fontsize=16)
    fig.suptitle(f"PCNA, LE = {lab_eff_reference}", fontsize=18, y=1.02)
    filename = my_experiment_mod.date_as_string + 'vlab4mic_PCNA_gallery_SMLM_modalities_' + str(lab_eff_reference) + 'LE.png'
    filename2 = os.path.join(my_experiment_mod.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()
