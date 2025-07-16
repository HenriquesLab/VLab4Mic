import supramolsim.generate.labelled_instance as lp
from supramolsim import workflows
from supramolsim.utils import data_format
import pytest
import copy
import numpy as np

def test_empty_particle():
    particle = lp.LabeledInstance()


def test_create_particle():
    particle = lp.create_particle()


# test with source and label params from fixture

structure_list = [
    "7R5K",
]
labels_per_structure = {
    "7R5K": [
        "NPC_Nup96_Cterminal_direct",
    ]
}
generic_labels = [
    "NHS_ester",
]


@pytest.mark.parametrize(
    "structure_id, structure_label",
    [(key, val) for key, key_vals in labels_per_structure.items() for val in key_vals],
)
def test_add_specific_labels(structure_id, structure_label, configuration_directory):
    configuration_path = configuration_directory
    fluorophore_id = "AF647"
    # loading structure
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_path
    )
    labels_list = []
    labels_list.append(
        data_format.structural_format.label_builder_format(
            structure_label, fluorophore_id
        )
    )
    particle, label_params_list = workflows.particle_from_structure(
        structure, labels_list, configuration_path
    )
    assert particle.get_ref_point().shape == (3,)


@pytest.mark.parametrize("structure_id", structure_list)
@pytest.mark.parametrize("generic_label", generic_labels)
def test_add_generic_labels(structure_id, generic_label, configuration_directory):
    configuration_path = configuration_directory
    fluorophore_id = "AF647"
    # loading structure
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_path
    )
    labels_list = []
    labels_list.append(
        data_format.structural_format.label_builder_format(
            generic_label, fluorophore_id
        )
    )
    particle, label_params_list = workflows.particle_from_structure(
        structure, labels_list, configuration_path
    )
    assert particle.get_ref_point().shape == (3,)


def test_add_probe_with_peptide_motif(experiment_7r5k_base):
    test_experiment = copy.copy(experiment_7r5k_base)
    test_experiment.remove_probes()
    list_of_proteins = test_experiment.structure.list_protein_names()
    chain_name = list_of_proteins[np.random.randint(0, len(list_of_proteins))]
    probe_name = "Antibody"
    test_experiment.add_probe(
        probe_name=probe_name,
        peptide_motif={
            "chain_name": chain_name,
            "position": "cterminal",
        },
    )
    assert len(test_experiment.probe_parameters) == 1
    assert test_experiment.probe_parameters[probe_name]["label_name"] == probe_name
    assert test_experiment.probe_parameters[probe_name]["target_info"]["type"] == "Sequence"
    assert type(test_experiment.probe_parameters[probe_name]["target_info"]["value"]) == str