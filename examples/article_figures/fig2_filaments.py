import numpy as np
import matplotlib.pyplot as plt
from vlab4mic import experiments
import copy
import seaborn as sns
import os
from vlab4mic.analysis import particle_measures
import pandas as pd

structure = "5SYF"
probe_model = "1IGT"          # used for every primary and secondary antibody
labelling_efficiency = 1
random_seed = 24  
primary_1 = dict(
    probe_template="Antibody",
    probe_name="ATubulin_epitopes",
    probe_distance_to_epitope = 0,
    probe_model=probe_model,
    probe_target_type="Sequence",
    probe_target_value="EEVGV",
    labelling_efficiency=1,
)


sample_dimensions = [1000, 1000, 10]
common = dict(
    structure=structure,
    structure_global_normal_orientation="local_plane",   # normals from a local plane fit
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

pts = line_of_points2(sample_dimensions[0], sample_dimensions[1], spacing=7, line_position=(0.5, 0.5), orientation=0, length=500)

__, noiseless_, my_experiment = experiments.image_vsample(
    probe_list=[primary_1],
    **common,
    random_orientations=True,
    xz_orientations=[90,],
)
my_experiment.clear_virtual_sample()
my_experiment.set_virtualsample_params(
    particle_positions=pts)
my_experiment.build(modules=["coordinate_field","imager"])

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], emitters_plotsize=1, axis_object=ax1)
ax2 = fig.add_subplot(1,2,2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[0,0,0], emitters_plotsize=1, axis_object=ax2)
filename = my_experiment.date_as_string + 'vlab4mic_microtubule_filament.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 


n_molecules = len(my_experiment.coordinate_field.molecules)
epitopes_new = np.empty(shape=(0,3))
normals_new = np.empty(shape=(0,3))
for particle_number in range(n_molecules):
    epitopes_new = np.vstack((epitopes_new, my_experiment.coordinate_field.molecules[particle_number].source["targets"]["ATubulin_epitopes"]["coordinates"] ))
    normals_new = np.vstack((normals_new, my_experiment.coordinate_field.molecules[particle_number].source["targets"]["ATubulin_epitopes"]["normals"] ))

assert epitopes_new.shape == (1872,3)

absolute_positions = epitopes_new.copy()
normals = normals_new.copy()
particle_generator_data = dict(
    targets={
        "ATubulin_epitopes":dict(
            coordinates = absolute_positions,
            normals=normals
        )
    },
    scale=1e-9,
    reference_point = my_experiment.coordinate_field.get_field_param("absolute_reference_point").copy(),
    axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0]))
    )
probe_parameters_custom = copy.deepcopy(my_experiment.particle.labels["ATubulin_epitopes"])

### Filament diameter
n_replicates = 100
labelling_efficiency_primary = 1
labelling_efficiency_secondary = 1

primary_1 = dict(
    probe_template="Antibody",
    probe_name="ATubulin_epitopes",
    probe_distance_to_epitope = 0,
    probe_model=probe_model,
    probe_target_type="Sequence",
    probe_target_value="EEVGV",
    labelling_efficiency=1,
    probe_paratope= "TITCH",
    probe_secondary_epitope = "EKKNW"
)
secondary_1 = dict(
    probe_template="Antibody",
    probe_name="Secondary_L1",
    probe_model=probe_model,
    probe_target_type="Primary",
    probe_target_value="ATubulin_epitopes",
    probe_conjugation_target_info=dict(
        type="Atom_residue",
        value=dict(atoms=["CA"], residues=["LYS"]),
    ),
    probe_paratope= "TITCH", 
    #probe_DoL=2,
    labelling_efficiency=labelling_efficiency_secondary,
)
secondary_2 = dict(
    probe_template="Antibody",
    probe_name="Secondary_L2",
    probe_model=probe_model,
    probe_target_type="Primary",
    probe_target_value="ATubulin_epitopes",
    probe_conjugation_target_info=dict(
        type="Atom_residue",
        value=dict(
            atoms= ["CA", ],
            residues=["GLN", ],
            position=305
            ),
    ),
    probe_DoL=None,
    labelling_efficiency=labelling_efficiency,
)
### First strategy
__, noiseless_, my_experiment_probe1 = experiments.image_vsample(
    probe_list=[primary_1],
    **common
)
my_experiment_probe1.particle.source = dict()
my_experiment_probe1.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe1.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe1.particle.generate_instance()
my_experiment_probe1.build(modules=["coordinate_field", "imager"])
replicates_1 = my_experiment_probe1.run_replicates(n_replicates=n_replicates, modality="SMLM")
## Second strategy
__, noiseless_, my_experiment_probe2 = experiments.image_vsample(
    primary_probe=primary_1,
    secondary_probe=secondary_1,
    **common
)
my_experiment_probe2.particle.source = dict()
my_experiment_probe2.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe2.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe2.particle.generate_instance()
my_experiment_probe2.build(modules=["coordinate_field", "imager"])
replicates_2 = my_experiment_probe2.run_replicates(n_replicates=n_replicates, modality="SMLM")
## Third
__, noiseless_, my_experiment_probe3 = experiments.image_vsample(
    primary_probe=primary_1,
    secondary_probe=secondary_2,
    **common
)
my_experiment_probe3.particle.source = dict()
my_experiment_probe3.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe3.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe3.particle.generate_instance()
my_experiment_probe3.build(modules=["coordinate_field", "imager"])
replicates_3 = my_experiment_probe3.run_replicates(n_replicates=n_replicates, modality="SMLM")

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,3,1)
ax1.imshow(replicates_1[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax1.set_title("Primary Antibody")
ax2 = fig.add_subplot(1,3,2)
ax2.imshow(replicates_2[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax2.set_title("Primary Antibody + Secondary Antibody (NHS)")
ax3 = fig.add_subplot(1,3,3)
ax3.imshow(replicates_3[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax3.set_title("Primary Antibody + Secondary Antibody (site-specific)")
filename = my_experiment_probe3.date_as_string + 'vlab4mic_microtubules_simulation_examples.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 

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

filaments_params = dict(
    pixel_size=2, 
    max_length_nm=500, 
    modality="SMLM", 
    channel="ch0", 
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=2,
    threshold_method = "otsu"
)

df1, filament_metrics_list1, masks1, blurred1, tables1 = measure_filaments(
    replicates=replicates_1, 
    **filaments_params
)
df2, filament_metrics_list2, masks2, blurred2, tables2 = measure_filaments(
    replicates=replicates_2, 
    **filaments_params
    )
df3, filament_metrics_list3, masks3, blurred3, tables3 = measure_filaments(
    replicates=replicates_3, 
    **filaments_params
    )

combined_df_sw = pd.DataFrame(data={
    "Primary": df1["s_width"],
    "Primary_Secondary_SiteSpec": df3["s_width"],
    "Primary_Secondary_NHS": df2["s_width"],
    }
)
custom_palette = {
    "Primary": "#F58B19",
    "Primary_Secondary_SiteSpec": "green",
    "Primary_Secondary_NHS": "#a34bb0",
}

fig = plt.figure()
ax = fig.add_subplot(111)
filament_diameters = sns.histplot(data=combined_df_sw, palette=custom_palette, bins=30, ax=ax)
filename = my_experiment_probe3.date_as_string + 'vlab4mic_apparent_diameter_microtubule.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 
fig.close()

### Apparent continuity

labelling_efficiencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1]
replicates_per_le = dict()
n_replicates_le = 100
for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="Antibody",
        probe_name="ATubulin_epitopes",
        probe_distance_to_epitope = 0,
        probe_model=probe_model,
        probe_target_type="Sequence",
        probe_target_value="EEVGV",
        labelling_efficiency=labeff,
        probe_paratope= "TITCH", # sequence in 1IGT
        probe_secondary_epitope = "EKKNW"
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
    max_length_nm=500, 
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
sns.swarmplot(data=merged, x="source", y="mean_eccentricity", ax=ax)
filename = my_experiment_le.date_as_string + 'vlab4mic_mean_eccentricities_per_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
plt.xlabel("Labelling efficiency")
fig.savefig(filename2, dpi=300, bbox_inches='tight')
fig.close()