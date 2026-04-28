from vlab4mic import experiments
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import os
########################## 
random_seed = 24
dome_model = "YOUR/PATH/TO/CIF"
flat_model = "YOUR/PATH/TO/CIF"
primary = dict(
    probe_template = "Antibody",
    probe_name="custom",
    probe_target_type = "Sequence",
    probe_target_value = "NNRIA",
    probe_distante_to_epitope = 0,
    probe_DoL=4,  
)
############################ DOME #######################################
images, noiseless, my_experiment = experiments.image_vsample(
    structure=dome_model,
    structure_is_path=True,
    probe_list=[primary,],
    clear_experiment=True,
    multimodal=["Widefield", "Confocal", "AiryScan", "STED", "SMLM"],
    run_simulation=False,
    sample_dimensions=[500,500,10],
    random_seed=random_seed
)
my_experiment.update_modality(modality_name="STED", depth_of_field_nm=1000)
my_experiment.update_modality(modality_name="SMLM", depth_of_field_nm=1000)
images, noiseless = my_experiment.run_simulation()
fig = plt.figure(figsize=[10,20])
ax = fig.add_subplot(1, 3, 1, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 50000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[90,0,0],
    target_size = 0,
    atoms_size = 1,
    atoms_alpha = 0.01,
    reference_point = False,
    axesoff=True,
    )
ax = fig.add_subplot(1, 3, 2, projection="3d")
my_experiment.structure.show_target_labels(
    with_assembly_atoms=False, 
    axis_object=ax, 
    show_axis=False, 
    view_init=[90,0,0],
    target_size = 10,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 3, projection="3d")
my_experiment.particle.gen_axis_plot(
 axis_object=ax,
 view_init=[90,0,0],
 emitter_plotsize=5,
 xlim=[-1250,1250], ylim=[-1250,1250], zlim=[0,1000]
)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_DOME_labelled_structure.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=[20,30])
for i, mod_name in enumerate(images.keys()):
    ax = fig.add_subplot(1,5,i+1)
    ax.imshow(images[mod_name]["ch0"][0], cmap="grey")
    ax.set_title(mod_name)
    ax.set_axis_off()
length_nm = 200
nm = 1e-9
pixelsize = (my_experiment.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * my_experiment.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = 5 
scalebar = AnchoredSizeBar(ax.transData,
                           length_px, 
                           "",
                           "lower right",
                           borderpad=0.5,
                           color='white',
                           frameon=False,
                           size_vertical=hight_px)
ax.add_artist(scalebar)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_DOME_image_simulation.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

########################## FLAT MODEL #########################

my_experiment.select_structure(
    structure_path=flat_model,
    build=False
)
my_experiment.set_virtualsample_params(
    random_rotations=True,
    rotation_angles=[85],
)
my_experiment.build()
images, noiseless = my_experiment.run_simulation()

fig = plt.figure(figsize=[10,20])
ax = fig.add_subplot(1, 3, 1, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 50000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[90,0,0],
    target_size = 0,
    atoms_size = 1,
    atoms_alpha = 0.01,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 2, projection="3d")
my_experiment.structure.show_target_labels(
    with_assembly_atoms=False, 
    axis_object=ax, 
    show_axis=False, 
    view_init=[90,0,0],
    target_size = 10,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 3, projection="3d")
my_experiment.particle.gen_axis_plot(
 axis_object=ax,
 view_init=[90,0,0],
 emitter_plotsize=5,
 xlim=[-1250,1250], ylim=[-1250,1250], zlim=[0,1000],
)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_FLAT_labelled_structure.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=[20,30])
for i, mod_name in enumerate(images.keys()):
    ax = fig.add_subplot(1,5,i+1)
    ax.imshow(images[mod_name]["ch0"][0], cmap="grey")
    ax.set_title(mod_name)
    ax.set_axis_off()
length_nm = 200
nm = 1e-9
pixelsize = (my_experiment.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * my_experiment.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = 5 
scalebar = AnchoredSizeBar(ax.transData,
                           length_px, 
                           "",
                           "lower right",
                           borderpad=0.5,
                           color='white',
                           frameon=False,
                           size_vertical=hight_px)
ax.add_artist(scalebar)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_FLAT_image_simulation.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()