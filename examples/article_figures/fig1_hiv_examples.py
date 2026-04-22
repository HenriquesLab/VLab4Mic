from vlab4mic import experiments
import matplotlib.pyplot as plt
import os
seed = 1


images, noisless, my_experiment = experiments.image_vsample(
    structure="3J3Y",
    probe_template="anti-p24_primary_antibody_HIV",
    probe_DoL=3,
    structural_integrity_small_cluster=20,
    structural_integrity_large_cluster=100,
    sample_dimensions=[1000,1000,100],
    random_orientations=True,
    structural_integrity=1,
    clear_experiment=True,
    number_of_particles=30,
    multimodal=["Widefield", "Confocal", "AiryScan", "STED", "SMLM", "Reference"],
    random_seed=seed,
    run_simulation=False,
)
my_experiment.set_virtualsample_params(minimal_distance=110)
my_experiment.build(modules=["coordinate_field"])
my_experiment.update_modality(modality_name="SMLM",
                              pixelsize_nm=2,
                              psf_voxel_nm=2,
                              depth_of_field_nm=1000,
                              lateral_resolution_nm=5,
                              simulate_localistations=True,
                              lateral_precision=5,
                              axial_precision=5,
                              nlocalisations=20)
my_experiment.update_modality(modality_name="STED", depth_of_field_nm=1000, lateral_resolution_nm=30)
my_experiment.set_modality_acq(modality_name="STED", exp_time=0.001)
my_experiment.update_modality(modality_name="Reference",
                              pixelsize_nm=2,
                              psf_voxel_nm=2,
                              depth_of_field_nm=1000,
                              lateral_resolution_nm=2,
                              simulate_localistations=True,
                              lateral_precision=0.001,
                              axial_precision=0.001,
                              nlocalisations=30)
my_experiment.set_modality_acq(modality_name="Reference", exp_time=1)
my_experiment.build(modules=["imager"])
images, noiseless = my_experiment.run_simulation()
modalities = list(images.keys())
n_modalities = len(modalities)
fig = plt.figure(figsize=[40,10])

for i, name in enumerate(modalities):
    ax = fig.add_subplot(1, n_modalities, i+1)
    ax.imshow(images[name]["ch0"][0], cmap="grey")
    ax.set_axis_off()

filename = my_experiment.date_as_string + 'vlab4mic_hiv_modalities.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()


target_colour="#01579D"
fig = plt.figure(figsize=[20,10])
ax = fig.add_subplot(1, 2, 1, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 100000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[20,0,0],
    target_size = 0,
    atoms_size = 10,
    atoms_alpha = 0.01,
    reference_point = False,
    target_plotcolour=target_colour
    )
ax = fig.add_subplot(1, 2, 2, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 100000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[20,0,0],
    target_size = 15,
    atoms_size = 10,
    atoms_alpha = 0.01,
    reference_point = False,
    target_plotcolour=target_colour
    )
probe_template = list(my_experiment.probe_parameters.keys())[0]
title = "Structure: 3J3Y \n Probe: " + probe_template
title = title + "\n" + "Sequence: " + my_experiment.probe_parameters[probe_template]["target"]["value"]
ax.set_title(title)
filename = my_experiment.date_as_string + 'vlab4mic_hiv_structure_epitopes.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=[5,5])
ax = fig.add_subplot(1, 1, 1, projection="3d")
my_experiment.particle.show_probe(
    plotcolour="#984ea3", axesoff=True,
    with_structural_atoms=True,
    emitter_plotsize=70,
    atoms_plotalpha=0.05,
    atoms_plotsize=10,
    atoms_fraction=0.3,
    view_init=[0,40,100],
    central_axis=False,
    use_dol=True,
    axis_object = ax,
    )

filename = my_experiment.date_as_string + 'vlab4mic_hiv_antibody.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()


target_colour="#01579D"
plotmarker="o"
emitter_plotmarker = "o"
hiv_antibody=dict(
        source_plotsize=3,
        emitter_plotsize=50,
        source_plotcolour=target_colour,
        source_plotalpha=1,
        source_plotmarker=plotmarker,
        emitter_plotmarker=emitter_plotmarker
    )
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(1, 1, 1, projection="3d")
my_experiment.particle.gen_axis_plot(
        reference_point=False,
        view_init=[20, 0, 0],
        axesoff=True,
        show_axis=False,
        with_sources=True,
        axis_object=ax,
        with_normals=False,
        **hiv_antibody)
filename = my_experiment.date_as_string + 'vlab4mic_hiv_labelled_particle.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

n_examples = 8
fig = plt.figure(figsize=[20, 10])
row=2
my_experiment.set_structural_integrity(
            structural_integrity=1,
            structural_integrity_small_cluster=20,
            structural_integrity_large_cluster=100,
        )
my_experiment.build(modules=["particle",])
for i in range(n_examples):
    ax = fig.add_subplot(row, int(n_examples/2), i+1, projection="3d")
    if i == 4:
        my_experiment.set_structural_integrity(
            structural_integrity=0.5
        )
        my_experiment.build(modules=["particle",])

    my_experiment.particle.generate_instance()
    my_experiment.particle.gen_axis_plot(
        reference_point=False,
        view_init=[20, 0, 0],
        axesoff=True,
        show_axis=False,
        with_sources=True,
        axis_object=ax,
        with_normals=False,
        source_plotsize=1,
        source_plotalpha=0.4,
        emitter_plotsize=5,
        )

filename = my_experiment.date_as_string + 'vlab4mic_hiv_labelled_particle_examples.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()