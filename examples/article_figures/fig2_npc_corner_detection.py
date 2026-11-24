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

RANDOM_SEED = 7
PRECISION_NM = 10
MIN_SEPARATION_NM = 2 * PRECISION_NM  # peak separation and DBSCAN radius
DBSCAN_MIN_SAMPLES = 3
PEAK_THRESHOLD_REL = 0.2  # default of count_resolved_sites
N_SEGMENTS = 8  # maximum number of segments of the ring in 2D
nlocalisations=10
EXAMPLE_REPLICATE = 0
RENDERING_KERNEL = 3
depth_of_field_nm = 20
depth_of_field_nm = 1000

labelling_efficiencies = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
replicates_per_le = dict()
parameters_le = dict()
n_replicates_le = 100
probe_DoL = None
import copy

for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="GFP_w_nanobody",
        probe_name="GFP_w_nanobody",
        probe_target_type="Sequence",
        probe_target_value="ELAVGSL",
        probe_distance_to_epitope = 0,
        labelling_efficiency=labeff,
        probe_DoL=probe_DoL,
    )
    # GFP_w_nanobody has no default target: Nup96 C-terminus
    _, _, my_experiment_le = experiments.image_vsample(
        structure="7R5K",
        **primary_le,
        random_rotations=True,
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[500, 500, 1],
        sample_inital_orientation=[0, 0, 1],  # pore seen from above
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED,
    )
    my_experiment_le.update_modality(
        modality_name="SMLM",
        simulate_localistations=True,
        nlocalisations=nlocalisations,
        depth_of_field_nm=depth_of_field_nm,
        lateral_precision=PRECISION_NM,
        axial_precision=PRECISION_NM,
        precision_model="fixed",
        lateral_resolution_nm=RENDERING_KERNEL,
    )
    
    replicates_le = my_experiment_le.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le[labeff] = copy.copy(replicates_le)
    parameters_le[labeff] = my_experiment_le.get_parameters()

# all outputs of this script go in one folder named after it
output_dir = os.path.join(my_experiment_le.output_directory,
                          os.path.splitext(os.path.basename(__file__))[0])
os.makedirs(output_dir, exist_ok=True)


#######
# split_z: False counts corners on both rings of the pore together,
# True keeps only the ring above the mean z (one side of the pore)
SPLIT_Z_CONDITIONS = {"Both sides of NPC": False, "One side of NPC": True}
min_locs_per_segment = 4
n_corners = 8
sectorcounting_params_le = dict(
    modality="SMLM",
    channel="ch0",
    n_corners=n_corners,
    radius=52,
    r_min=30,
    r_max=70,
    min_per_section=min_locs_per_segment,
    min_locs_to_circlefit=nlocalisations,
    max_radius=80,
)
# scalar columns of the ring_corners_from_replicates table, saved as csv
CSV_COLUMNS = ["replicate", "n_localisations", "radius_fit", "n_corners", "n_segments_occupied"]

counting_measures = dict()  # sides -> labelling efficiency -> dataframe, images, rows
for sides, split_z in SPLIT_Z_CONDITIONS.items():
    counting_measures[sides] = dict()
    for lef in labelling_efficiencies:
        results, images_list, rows = ring_corners_from_replicates(
            replicates=replicates_per_le[lef],
            split_z=split_z,
            **sectorcounting_params_le
        )
        counting_measures[sides][lef] = dict(dataframe=results, images=images_list, rows=rows)

merged = pd.concat(
    [counting_measures[sides][lef]["dataframe"][CSV_COLUMNS].assign(
        sides=sides, split_z=split_z, labelling_efficiency=lef)
     for sides, split_z in SPLIT_Z_CONDITIONS.items() for lef in labelling_efficiencies],
    ignore_index=True,
)
merged["radius_fit"] = merged["radius_fit"].astype(float)
summary = merged.groupby(["sides", "labelling_efficiency"], sort=False).agg(
    n_corners_mean=("n_corners", "mean"), n_corners_sd=("n_corners", "std"),
    radius_fit_mean=("radius_fit", "mean"), radius_fit_sd=("radius_fit", "std"),
    n_radius_fit=("radius_fit", "count"), n_total=("n_corners", "size"),
)
prefix = os.path.join(output_dir,
                      my_experiment_le.date_as_string + 'vlab4mic_npc_corners_')
merged.to_csv(prefix + "particles.csv", index=False)
summary.to_csv(prefix + "summary.csv")
print(summary)


######## plot 1, both split_z conditions in one boxplot
fig, ax = plt.subplots(figsize=(7, 3.5))
sns.boxplot(data=merged, x="labelling_efficiency", y="n_corners", hue="sides", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Number of Corners detected")
plt.title("Corners detected on NPC simuations" + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
ax.legend(title=None)
filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_NPC_simuations.png'
filename2 = os.path.join(output_dir, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

# plot 6, both split_z conditions in one boxplot
fig, ax = plt.subplots(figsize=(7, 3.5))
sns.boxplot(data=merged, x="labelling_efficiency", y="radius_fit", hue="sides", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Radius fit")
ax.legend(title=None)
filename = my_experiment_le.date_as_string + 'vlab4mic_npc_radius_fit_LE.png'
filename2 = os.path.join(output_dir, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

for sides, split_z in SPLIT_Z_CONDITIONS.items():
    # one output per split_z condition, one side marked in the file name
    suffix = '_one_side' if split_z else ''

    #### plots 2 and 3, mean and std shown for LE 0.3
    for lab_eff_reference in [0.3, 0.4]:
        counts = counting_measures[sides][lab_eff_reference]["dataframe"].n_corners
        hist = np.bincount(counts, minlength=n_corners + 1)
        label = "detected"
        if lab_eff_reference == 0.3:
            label += f"\nmean ± std = {counts.mean():.2f} ± {counts.std():.2f}"
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.bar(np.arange(n_corners + 1), hist, label=label)
        ax.set_xlabel("corners per NPC")
        ax.set_ylabel("NPCs")
        if lab_eff_reference == 0.3:
            ax.set_ylim(top=hist.max() * 1.35)  # room for the mean and std in the legend
        ax.legend()
        ax.set_title("Corners detected on NPC simuations. LE = " + str(lab_eff_reference) + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
        filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_NPC_simuations_' + str(lab_eff_reference) + 'LE' + suffix + '.png'
        filename2 = os.path.join(output_dir, filename)
        fig.savefig(filename2, dpi=300, bbox_inches='tight')
        plt.close()



### plot 4
# the images do not depend on split_z
images_per_le = counting_measures["Both sides of NPC"]
n_le = len(labelling_efficiencies)
nrows = 5
fig, ax = plt.subplots(nrows, n_le, figsize=(30,10))
for row in range(nrows):
    for i, lef in enumerate(labelling_efficiencies):
        ax[row,i].imshow(images_per_le[lef]["images"][row][0], cmap = "grey")
        ax[row,i].set_title(f"LE = {lef}")
        ax[row,i].set_xticks([])
        ax[row,i].set_yticks([])
filename = my_experiment_le.date_as_string + 'vlab4mic_NPCs_across_LE.png'
filename2 = os.path.join(output_dir, filename)
if probe_DoL is None:
    maintitle = "probe_DoL: max" + ". DoF: " + str(depth_of_field_nm)
else:
    maintitle = "probe_DoL: " + str(probe_DoL) + ". DoF: " + str(depth_of_field_nm)
fig.suptitle(maintitle)
plt.tight_layout()
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()


## plot 5
# e.g. pores where only 3 or fewer corners were detected
#threshold_for_plot = 4
threshold_for_plot = [4, 6, 8]
lab_eff_reference = 0.3
examples = []  # (file name suffix, row)
for sides, split_z in SPLIT_Z_CONDITIONS.items():
    _df = counting_measures[sides][lab_eff_reference]["dataframe"]
    for ncorners in threshold_for_plot:
        ids = list(_df.loc[_df["n_corners"] == ncorners].replicate)
        if not ids:
            continue
        suffix = '_one_side' if split_z else ''
        examples.append((str(ncorners) + suffix, counting_measures[sides][lab_eff_reference]["rows"][ids[0]]))
# same axis limits on every example, wide enough to show all their localisations
example_lim = 1.3 * sectorcounting_params_le["r_max"]
for _, row in examples:
    if row["centre"] is not None and len(row["points"]):
        dist = np.abs(row["points"][:, :2] - row["centre"]).max()
        example_lim = max(example_lim, 1.05 * dist)
for name, row in examples:
    fig = plot_ring_corners(row, lim=example_lim)
    filename = my_experiment_le.date_as_string + 'vlab4mic_npc_corners_example_' + name + '.png'
    filename2 = os.path.join(output_dir, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()