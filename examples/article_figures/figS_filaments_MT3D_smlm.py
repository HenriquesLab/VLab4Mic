import numpy as np
import matplotlib.pyplot as plt
from vlab4mic import experiments
import copy
import seaborn as sns
import os
import io
import json
import zipfile
import requests
from vlab4mic.analysis import particle_measures
import pandas as pd

# Experimental 3D SMLM localisations of microtubules (ShareLoc.xyz, Zenodo 5541063)
smlm_experimental_path = "https://zenodo.org/records/5541063/files/MT3D_160528_C1C2_1424K/data.smlm"

structure = "5SYF"
probe_model = "1IGT"
random_seed = 24
sample_dimensions = [1000, 1000, 10]
common = dict(
    structure=structure,
    sample_dimensions=sample_dimensions,                    # nm
    number_of_particles=1,                               # a single labelled particle
    modality="SMLM",
    SMLM=dict(
        lateral_precision=5,        # nm
        axial_precision=5,          # nm
        n_localisations=30,
        precision_model="fixed",    # not photon-limited
    ),
    clear_experiment=True,
    random_seed=random_seed,
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

def crop_filament_segment(localisations, center, length=500, width=150, search_radius=400):
    """
    Crop a straight filament segment and align it to the x axis.

    The filament direction is the principal axis of the localisations within
    search_radius of center (xy, nm). Localisations within length/2 along and
    width/2 across this axis are kept, rotated so the filament is horizontal
    and centred at the origin (z centred on its median).
    """
    xyz = localisations[["x", "y", "z"]].to_numpy(dtype=float)
    near = np.hypot(xyz[:, 0] - center[0], xyz[:, 1] - center[1]) < search_radius
    segment_center = xyz[near, :2].mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov((xyz[near, :2] - segment_center).T))
    direction = eigenvectors[:, -1]
    rotation = np.array([[direction[0], direction[1]], [-direction[1], direction[0]]])
    along_across = (xyz[:, :2] - segment_center) @ rotation.T
    keep = (np.abs(along_across[:, 0]) <= length / 2) & (np.abs(along_across[:, 1]) <= width / 2)
    segment = np.zeros((keep.sum(), 3))
    segment[:, :2] = along_across[keep]
    segment[:, 2] = xyz[keep, 2] - np.median(xyz[keep, 2])
    return segment, keep, segment_center, direction

localisations_all = load_smlm_table(smlm_experimental_path)

# isolated, straight microtubule; 500 nm long as the simulated filament in fig2_filaments.py
segment_center = (26500, 15200)    # nm
segment_length = 500               # nm
segment_width = 150                # nm
fluorophore_positions, in_segment, center_xy, direction_xy = crop_filament_segment(
    localisations_all, center=segment_center, length=segment_length, width=segment_width
)
print(f"{fluorophore_positions.shape[0]} localisations used as fluorophore positions")

### Fluorophore positions from the localisations
primary_1 = dict(
    probe_template="Antibody",
    probe_name="ATubulin_epitopes",
    probe_distance_to_epitope = 0,
    probe_model=probe_model,
    probe_target_type="Sequence",
    probe_target_value="EEVGV",
    probe_steric_hindrance=0,     # localisations are not excluded by probe size
    labelling_efficiency=1,
)
# without normals the targets are labelled directly: each selected
# localisation becomes a fluorophore at its own position
particle_generator_data = dict(
    targets={
        "ATubulin_epitopes":dict(
            coordinates = fluorophore_positions.copy(),
            normals=None
        )
    },
    scale=1e-9,
    reference_point = np.array([0.0, 0.0, 0.0]),
    axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0]))
    )

__, noiseless_, my_experiment = experiments.image_vsample(
    probe_list=[primary_1],
    **common
)
my_experiment.particle.source = dict()
my_experiment.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False
)
my_experiment.particle.generate_instance()
my_experiment.build(modules=["coordinate_field", "imager"])

# overview of the dataset and the cropped segment
fig, axs = plt.subplots(1, 2, figsize=(20, 10))
overview, xedges, yedges = np.histogram2d(
    localisations_all.x, localisations_all.y, bins=int(localisations_all.x.max() / 40)
)
axs[0].imshow(overview.T, cmap="grey", vmax=np.percentile(overview, 99.5), origin="upper",
              extent=[xedges[0], xedges[-1], yedges[-1], yedges[0]])
axs[0].scatter(localisations_all.x[in_segment], localisations_all.y[in_segment], s=0.5, c="#F58B19")
axs[0].set_title("Experimental localisations (segment in orange)")
zoom, xedges, yedges = np.histogram2d(
    fluorophore_positions[:, 0], fluorophore_positions[:, 1],
    bins=[np.arange(-segment_length / 2, segment_length / 2 + 5, 5), np.arange(-segment_width / 2, segment_width / 2 + 5, 5)]
)
axs[1].imshow(zoom.T, cmap="grey", origin="lower", extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
axs[1].set_title("Cropped segment, aligned to x (nm)")
filename = my_experiment.date_as_string + 'vlab4mic_MT3D_segment_overview.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], emitters_plotsize=1, axis_object=ax1)
ax2 = fig.add_subplot(1,2,2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[0,0,0], emitters_plotsize=1, axis_object=ax2)
filename = my_experiment.date_as_string + 'vlab4mic_MT3D_filament.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

def measure_filaments(replicates,
                      max_length_nm=500,
                      pixel_size=1,
                      modality = "SMLM",
                      channel="ch0",
                      fill_holes = True,
                      smooth_sigma=1,
                      row_binning=1,
                      min_snr=2,
                      mode="single",
                      threshold_method = "otsu",
                      threshold_method_kwargs = None,):
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
            apparent_continuity = (major_lenght / max_length_nm)/ n_segments

        rows.append(dict(
            replicate=result["replicate"],
            mean_filament_width = filament_metrics["mean"],
            std_filament_width = filament_metrics["std"],
            median_filament_width = filament_metrics["median"],
            mask_continuity=apparent_continuity,
            mask_major_lenght = major_lenght,
            mask_number_of_segments = n_segments,
            mask_prefilter_smooth_sigma=smooth_sigma,
            pixelsize = pixel_size,
            row_binning=row_binning,
            s_width=s_width,
            min_eccentricity=min_eccentricity,
            max_eccentricity=max_eccentricity,
            mean_eccentricity=mean_eccentricity
        ))
    df = pd.DataFrame(rows)
    return df, filament_metrics_list, masks, blurred_list, tables

### Apparent continuity

labelling_efficiencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
replicates_per_le = dict()
n_replicates_le = 100
for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="Linker",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope = 0,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        probe_steric_hindrance=0,
        labelling_efficiency=labeff,
    )
    __, noiseless_, my_experiment_le = experiments.image_vsample(
        probe_list=[primary_le],
        **common
    )
    my_experiment_le.particle.source = dict()
    my_experiment_le.particle.load_source(**copy.deepcopy(particle_generator_data))
    my_experiment_le.set_virtualsample_params(
        sample_dimensions=sample_dimensions,
        update_mode = False
    )
    my_experiment_le.particle.generate_instance()
    my_experiment_le.build(modules=["coordinate_field", "imager"])
    replicates_le = my_experiment_le.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le[labeff] = copy.copy(replicates_le)

filament_measures_per_le = dict()
smooth_sigma = 0.0

filaments_params_le = dict(
    pixel_size=2,
    max_length_nm=segment_length,
    modality="SMLM",
    channel="ch0",
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=smooth_sigma,
    threshold_method = "otsu"
)

for i in labelling_efficiencies:
    df_, filament_metrics_list_, masks_, blurred_, tables_ = measure_filaments(
    replicates=replicates_per_le[i],
    **filaments_params_le
    )
    filament_measures_per_le[i] = dict()
    filament_measures_per_le[i]["dataframe"] = df_
    filament_measures_per_le[i]["filament_metrics_list"] = filament_metrics_list_
    filament_measures_per_le[i]["masks"] = masks_
    filament_measures_per_le[i]["blurred"] = blurred_
    filament_measures_per_le[i]["tables"] = tables_

dfs = dict()
for lef in labelling_efficiencies:
    dfs[lef] = filament_measures_per_le[lef]["dataframe"]

merged = pd.concat(
    [d.assign(source=name) for name, d in dfs.items()],
    ignore_index=True,
)
fig = plt.figure()
ax = fig.add_subplot(111)
sns.swarmplot(data=merged, x="source", y="mean_eccentricity", ax=ax, s=1)
filename = my_experiment_le.date_as_string + 'vlab4mic_MT3D_mean_eccentricities_per_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
plt.xlabel("Labelling efficiency")
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

n_le = len(labelling_efficiencies)
fig, ax = plt.subplots(1, n_le, figsize=(30,10))
for i, lef in enumerate(labelling_efficiencies):
    ax[i].imshow(filament_measures_per_le[lef]["blurred"][0], cmap = "grey")
    ax[i].set_title(f"LE = {lef}")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
filename = my_experiment_le.date_as_string + 'vlab4mic_MT3D_microtubules_across_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()
