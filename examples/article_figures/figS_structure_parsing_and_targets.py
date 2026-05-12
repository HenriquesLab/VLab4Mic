from vlab4mic import experiments
import matplotlib.pyplot as plt
from vlab4mic.utils.sample import arrays
import os
random_seed = 24

experiment_examples = dict()
structure_probe = dict(
    HIV={"structure_id":"3J3Y",
         "probe_template":"HIV_capsid_p24_direct"},
    NPC={"structure_id":"7R5K",
         "probe_template":"NPC_Nup96_Cterminal_direct"},
)
for structure_name, structure_parameters in structure_probe.items():
    vsample, experiment_examples[structure_name] = experiments.generate_virtual_sample(
        structure=structure_parameters["structure_id"],
        probe_template=structure_parameters["probe_template"],
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
        target_size = 0,
        atoms_size = 20,
        atoms_alpha = 0.03,
        reference_point = False,
        with_normals=False,
        )
    
plt.tight_layout()
filename = experiment_examples["HIV"].date_as_string + 'figS_structure_parsing_asymmetricUnits.png'
filename2 = os.path.join(experiment_examples["HIV"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()



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
        target_size = 20,
        target_plotcolour="#377eb8",
        atoms_size = 10,
        atoms_alpha = 0.02,
        reference_point = False,
        )
    probe_template = list(experiment_examples[experiment_name].probe_parameters.keys())[0]
    title = "Structure: " + experiment_name + " \n Probe: " + probe_template
    ax.set_title(title)

filename = experiment_examples["HIV"].date_as_string + 'figS_structure_specific_targets.png'
filename2 = os.path.join(experiment_examples["HIV"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

################### 
experiment_examples["NPC"].remove_probes()
experiment_examples["NPC"].add_probe(
    probe_template="NPC_Nup96_Cterminal_direct", 
)
experiment_examples["NPC"].add_probe(
    probe_template="NPC_Nup96_Cterminal_direct", 
    probe_name="NPC_Nup107_Cterminal_direct",
    probe_target_type="Sequence", 
    probe_target_value="GYEIQ")
experiment_examples["NPC"].build(modules=["particle"])
experiment_examples["NPC"].structure.plotting_params["NPC_Nup96_Cterminal_direct"]["plotcolour"] = "#377eb8"
experiment_examples["NPC"].structure.plotting_params["NPC_Nup107_Cterminal_direct"]["plotcolour"] = "#ff7f00"

fig = plt.figure(figsize=[20,10])
ax = fig.add_subplot(121, projection="3d")
experiment_examples["NPC"].structure.show_target_labels(
    with_assembly_atoms=True,
    assembly_fraction=0.005,
    show_axis=False,
    reference_point=False,
    target_size = 50,
    atoms_size = 20,
    atoms_alpha = 0.01,
    view_init=[90,0,0],
    axis_object=ax,
)
ax.set_title("Top view")
ax = fig.add_subplot(122, projection="3d")
experiment_examples["NPC"].structure.show_target_labels(
    with_assembly_atoms=True,
    assembly_fraction=0.005,
    show_axis=False,
    reference_point=False,
    target_size = 50,
    atoms_size = 20,
    atoms_alpha = 0.01,
    view_init=[0,0,0],
    axis_object=ax,
)
ax.set_title("Lateral view")
filename = experiment_examples["NPC"].date_as_string + 'figS_NPC_two_targets.png'
filename2 = os.path.join(experiment_examples["NPC"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

# Relabel with NHS
for structure_common_name, experiment in experiment_examples.items():
    experiment.remove_probes()
    experiment.add_probe(
        probe_template="NHS_ester",
        labelling_efficiency=1
    )
    experiment.build(modules=["particle",])

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
        target_plotcolour="#377eb8",
        target_size = 1,
        atoms_size = 0,
        atoms_alpha = 0.01,
        reference_point = False,
        )

filename = experiment_examples["NPC"].date_as_string + 'figS_structure_NHS.png'
filename2 = os.path.join(experiment_examples["NPC"].output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()