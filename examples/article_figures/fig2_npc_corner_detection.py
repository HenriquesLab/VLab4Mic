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



####### 
counting_measures_per_le = dict()
split_z = False
min_locs_per_segment = 4
sectorcounting_params_le = dict(
    modality="SMLM",
    channel="ch0",
    split_z = split_z,
    n_corners=8,
    radius=52,
    r_min=30,
    r_max=70,
    min_per_section=min_locs_per_segment,
    min_locs_to_circlefit=nlocalisations,
    max_radius=80,
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
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
plt.title("Corners detected on NPC simuations" + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_NPC_simuations.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

#### plot 2
n_corners = 8
lab_eff_reference = 0.3
counts = counting_measures_per_le[lab_eff_reference]["dataframe"].n_corners
hist =  np.bincount(counts, minlength=n_corners + 1)
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar(np.arange(9), hist, label="detected")
ax.set_xlabel("corners per NPC")
ax.set_ylabel("NPCs")
ax.legend()
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
ax.set_title("Corners detected on NPC simuations. LE = " + str(lab_eff_reference) + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_NPC_simuations_0.3LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()



### plot 3
from scipy.stats import binom

n_corners = 8
lab_eff_reference = 0.4
counts = counting_measures_per_le[lab_eff_reference]["dataframe"].n_corners
hist =  np.bincount(counts, minlength=n_corners + 1)

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar(np.arange(9), hist, label="detected")
ax.set_xlabel("corners per NPC")
ax.set_ylabel("NPCs")
ax.legend()
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
ax.set_title("Corners detected on NPC simuations. LE = " + str(lab_eff_reference) + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))
filename = my_experiment_le.date_as_string + 'vlab4mic_corners_detected_on_NPC_simuations_0.4LE.png'
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
filename = my_experiment_le.date_as_string + 'vlab4mic_NPCs_across_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
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
_df = counting_measures_per_le[lab_eff_reference]["dataframe"]
for ncorners in threshold_for_plot:
    corners = results[_df.n_corners == ncorners]
    ids = list(_df.loc[_df["n_corners"] == ncorners].replicate)
    for example in ids[:1]:
        fig = plot_ring_corners(counting_measures_per_le[lab_eff_reference]["rows"][example])
    filename = my_experiment_le.date_as_string + 'vlab4mic_npc_corners_example_' + str(ncorners) + '.png'
    filename2 = os.path.join(my_experiment_le.output_directory, filename)
    fig.savefig(filename2, dpi=300, bbox_inches='tight')
    plt.close()

# plot 6
fig, ax = plt.subplots(figsize=(5, 3.5))
sns.boxplot(data=merged, x="source", y="radius_fit", ax=ax)
plt.xlabel("Labelling efficiency")
plt.ylabel("Radius fit")
filename = my_experiment_le.date_as_string + 'vlab4mic_npc_radius_fit_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()