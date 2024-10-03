import supramolsim.generate.labelled_instance as lp
from supramolsim import workflows, data_format
import pytest


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
        "7R5K_Nup96_Cterminal_direct",
    ]
}
generic_labels = [
    "Generic_NHS_ester",
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
    particle = workflows.particle_from_structure(
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
    particle = workflows.particle_from_structure(
        structure, labels_list, configuration_path
    )
    assert particle.get_ref_point().shape == (3,)
