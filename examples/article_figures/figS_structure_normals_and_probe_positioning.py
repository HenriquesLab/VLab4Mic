from vlab4mic import experiments
import matplotlib.pyplot as plt
import os
random_seed = 24

structure_probe = dict(
    HIV={"structure_id":"3J3Y", 
         "probe_parameters": dict(
             probe_template="Antibody",
             probe_target_type="Sequence",
             probe_target_value= "SPRTLNA"
         )},
    NPC={"structure_id":"7R5K",
         "probe_parameters": dict(
             probe_template="Antibody",
             probe_target_type="Sequence",
             probe_target_value= "ELAVGSL"
         )},
    CCP={"structure_id":"1XI5",
         "probe_parameters": dict(
             probe_template="Antibody",
             probe_target_type="Sequence",
             probe_target_value= "EQATETQ"
         )},
)
experiment_examples = dict()

for structure_name, structure_parameters in structure_probe.items():
    vsample, experiment_examples[structure_name] = experiments.generate_virtual_sample(
        structure=structure_parameters["structure_id"],
        **structure_parameters["probe_parameters"],
        clear_experiment=True,
        random_seed=random_seed
    )


fig = plt.figure(figsize=[20,10])
n = len(list(experiment_examples.keys()))
for i, experiment_name in enumerate(experiment_examples.keys()):
    ax = fig.add_subplot(1, n+1, i+1, projection="3d")
    total = experiment_examples[experiment_name].structure.num_assembly_atoms
    fraction = 10000/total
    experiment_examples[experiment_name].structure.show_target_labels(
        with_assembly_atoms=True, 
        axis_object=ax, 
        show_axis=False, 
        assembly_fraction = fraction,
        view_init=[90,0,0],
        target_size = 5,
        target_plotcolour="#377eb8",
        atoms_size = 5,
        atoms_alpha = 0.01,
        reference_point = False,
        with_normals=True
        )
    probe_template = list(experiment_examples[experiment_name].probe_parameters.keys())[0]
    title = "Structure: " + experiment_name + " \n Probe: " + probe_template
    title = title + "\n" + "Sequence: " + experiment_examples[experiment_name].probe_parameters[probe_template]["target"]["value"]
    ax.set_title(title)

filename = experiment_examples["HIV"].date_as_string + 'figS_Structure_normals.png'
filename2 = os.path.join(experiment_examples["HIV"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=[20,10])
n = len(list(experiment_examples.keys()))
for i, experiment_name in enumerate(experiment_examples.keys()):
    ax = fig.add_subplot(1, n+1, i+1, projection="3d")
    total = experiment_examples[experiment_name].structure.num_assembly_atoms
    fraction = 10000/total
    experiment_examples[experiment_name].particle.gen_axis_plot(
        axis_object=ax, 
        show_axis=False, 
        view_init=[90,0,0],
        with_sources=True,
        emitter_plotsize=3,
        source_plotsize=1
        )
    probe_template = list(experiment_examples[experiment_name].probe_parameters.keys())[0]
    title = "Structure: " + experiment_name + " \n Probe: " + probe_template
    title = title + "\n" + "Sequence: " + experiment_examples[experiment_name].probe_parameters[probe_template]["target"]["value"]
    ax.set_title(title)
plt.tight_layout()

filename = experiment_examples["HIV"].date_as_string + 'figS_probe_positioning_to_normals.png'
filename2 = os.path.join(experiment_examples["HIV"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()