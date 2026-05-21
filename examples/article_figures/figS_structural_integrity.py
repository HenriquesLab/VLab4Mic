from vlab4mic import experiments
import matplotlib.pyplot as plt
random_seed = 24
import os

structure_probe = dict(
    HIV={"structure_id":"3J3Y", 
         "probe_parameters": dict(
             probe_template="Linker",
             probe_target_type="Sequence",
             probe_target_value= "SPRTLNA",
             probe_distance_to_epitope=0,
         ),
         "structural_integrity_params": dict(
             structural_integrity = 1,
             structural_integrity_small_cluster = 20,
            structural_integrity_large_cluster = 100
         )},
    NPC={"structure_id":"7R5K",
         "probe_parameters": dict(
             probe_template="Linker",
             probe_target_type="Sequence",
             probe_target_value= "ELAVGSL",
             probe_distance_to_epitope=0,
         ),
         "structural_integrity_params": dict(
             structural_integrity = 1,
             structural_integrity_small_cluster = 300,
            structural_integrity_large_cluster = 600
         )},
    CCP={"structure_id":"1XI5",
         "probe_parameters": dict(
             probe_template="Linker",
             probe_target_type="Sequence",
             probe_target_value= "EQATETQ",
             probe_distance_to_epitope=0,
         ),
         "structural_integrity_params": dict(
             structural_integrity = 1,
             structural_integrity_small_cluster = 100,
            structural_integrity_large_cluster = 200
         )},
)
n_structures = len(list(structure_probe.keys()))
experiment_examples = dict()

for structure_name, structure_parameters in structure_probe.items():
    imgs, noiseless, experiment_examples[structure_name] = experiments.image_vsample(
        structure=structure_parameters["structure_id"],
        **structure_parameters["probe_parameters"],
        **structure_parameters["structural_integrity_params"],
        clear_experiment=True,
        sample_dimensions=[400,400,10],
        random_seed=random_seed,
        modality="SMLM",
        run_simulation=False,
    )

images_per_structure = dict()
fig = plt.figure(figsize=[10,10])
fig2 = plt.figure(figsize=[10,10])
fig_number = 0
for structure_name, structure_experiment in experiment_examples.items():
    images_per_structure[structure_name] = dict()
    for fraction in [1.0, 0.5, 0.3]:
        structure_experiment.set_structural_integrity(
            structural_integrity = fraction,
        )
        structure_experiment.build(modules=["particle","coordinate_field", "imager"])
        ax = fig.add_subplot(n_structures, 3, fig_number+1, projection="3d")
        ax2 = fig2.add_subplot(n_structures, 3, fig_number+1)
        structure_experiment.coordinate_field.molecules[0].gen_axis_plot(
            axis_object=ax, 
            show_axis=False, 
            view_init=[90,0,0],
            with_sources=False,
            emitter_plotsize=5,
        )
        img, noisless = structure_experiment.run_simulation()
        ax2.imshow(img["SMLM"]["ch0"][0], cmap="grey")
        ax2.set_axis_off()
        fig_number+=1
plt.tight_layout()
    
filename_1 = experiment_examples["HIV"].date_as_string + 'figS_structural_integrity_fluorophores.png'
filename_2 = experiment_examples["HIV"].date_as_string + 'figS_structural_integrity_images.png'
filename1 = os.path.join(experiment_examples["HIV"].output_directory, filename_1)
filename2 = os.path.join(experiment_examples["HIV"].output_directory, filename_2)


fig.savefig(filename1,dpi=300, bbox_inches='tight')
plt.close()

fig2.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()