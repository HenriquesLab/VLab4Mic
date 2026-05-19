from vlab4mic import experiments
import matplotlib.pyplot as plt
import numpy as np
import os
random_seed = 24

modalities = ["SMLM",]
lat_axial_resolution_nm = 3 # the lateral and axial standard deviation of the PSF at rendering only
render_pixelsize_nm = 2
voxel_size = 2
number_of_locs_per_emitter = 30
depth_of_field_nm=500
sample_dimensions=[250,250,10]
axial_to_lateral_ratio = 3.2
secondary_dol = 1.5

# Primary and secondary
primary = dict(
    probe_template = "Antibody",
    probe_name="Primary-clathrin",
    probe_target_type = "Sequence",
    probe_target_value = "EQATETQ",
    labelling_efficiency = 1
)

secondary = dict(
    probe_template = "Antibody",
    probe_name="Secondary-clathrin",
    probe_target_type = "Primary",
    probe_target_value = "Primary-clathrin",
    probe_DoL = secondary_dol
)

# create the experiment
_1, _2, my_experiment = experiments.image_vsample(
    structure = "1XI5",
    primary_probe = primary,
    secondary_probe = secondary,
    multimodal=modalities,
    SMLM = dict(
        lateral_resolution_nm = lat_axial_resolution_nm,
        axial_resolution_nm = lat_axial_resolution_nm,
        pixelsize_nm=render_pixelsize_nm,
        psf_voxel_nm=voxel_size,
        depth_of_field_nm=depth_of_field_nm,
    ),
    sample_dimensions=sample_dimensions,
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)

# parameterise cases to show
locs_precision = [0,1,5,10,15, 20]
render_psf_sizes = [10,3]
number_of_locs_per_emitter = 30
images_per_rendersize_per_precision = dict()
for render_size in render_psf_sizes:
    render_size_str = str(render_size)
    images_per_rendersize_per_precision[render_size_str] = dict()
    for loc_p in locs_precision:
        if loc_p > 0:
            my_experiment.update_modality(
                modality_name="SMLM", 
                simulate_localistations=True, 
                lateral_precision=loc_p, 
                axial_precision=axial_to_lateral_ratio*loc_p, 
                lateral_resolution_nm=render_size,
                nlocalisations=number_of_locs_per_emitter)
            my_experiment.build(modules=["imager"])
            smlm_locs, noiseless = my_experiment.run_simulation()
            images_per_rendersize_per_precision[render_size_str][str(loc_p)] = smlm_locs
        else:
            my_experiment.update_modality(modality_name="SMLM", simulate_localistations=False)
            my_experiment.build(modules=["imager"])
            smlm_no_locs, noiseless = my_experiment.run_simulation()
            images_per_rendersize_per_precision[render_size_str][str(loc_p)] = smlm_no_locs


# plot in a grid
render_psf_sizes = [10,3]
n_psf_sizes = len(render_psf_sizes)
n_precision = len(locs_precision)
fig, axs = plt.subplots(nrows=n_psf_sizes, ncols=n_precision, squeeze=False, figsize=[10,5])
cmap = "grey" #"hot"
row = 0
for render_size in images_per_rendersize_per_precision.keys():
    col = 0
    set_ylabel = True
    for loc_p in images_per_rendersize_per_precision[render_size].keys():
        axs[row, col].imshow(images_per_rendersize_per_precision[render_size][loc_p]["SMLM"]["ch0"][0], cmap=cmap)
        if loc_p == "0":
            if set_ylabel:
                axs[row, col].set_title("No locs")
                axs[row, col].set_ylabel(f"Render size: {render_size}")
                axs[row, col].set_xticks([])
                axs[row, col].set_yticks([])
                set_ylabel = False
            else:
                axs[row, col].set_title("No locs")
                axs[row, col].set_axis_off()
        else:
            axs[row, col].set_title(f"Loc precision: {loc_p}")
            axs[row, col].set_axis_off()
        col+=1
    row+=1
plt.tight_layout()

filename = my_experiment.date_as_string + 'figS_localisations_simulation_ccs.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches="tight")