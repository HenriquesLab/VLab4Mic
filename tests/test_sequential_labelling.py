import pytest
from supramolsim import workflows, data_format




def test_linker_secondary(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id1 = "Mock_linker"
    target_info1 = dict(
        type="Sequence",
        value="EQATETQ"
    )
    fluorophore_id1 = "AF647"
    lab_eff1 = 1
    tmp_label1 = data_format.structural_format.label_builder_format(
        label_id1, fluorophore_id1, lab_eff1, target_info1
    )

    label_id2 = "Mock_primary_antibody"
    target_info2 = dict(
        type="Primary",
        value="Mock_linker"
    )
    fluorophore_id2 = "AF647"
    lab_eff2 = 1
    tmp_label2 = data_format.structural_format.label_builder_format(
        label_id2, fluorophore_id2, lab_eff2, target_info2
    )
    particle2, label_params_list = workflows.particle_from_structure(
        structure, [tmp_label1, tmp_label2], configuration_directory
    )
    assert len(particle2.secondary.keys()) > 0
    first_labels = particle2.emitters[target_info2["value"]].shape
    particle2.sequential_labelling=True
    particle2.generate_instance()
    second_labels = particle2.emitters[target_info2["value"]].shape
    assert first_labels[0] != second_labels[0]



def test_primary_secondary(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id1 = "Mock_primary_antibody"
    target_info1 = dict(
        type="Sequence",
        value="EQATETQ"
    )
    fluorophore_id1 = "AF647"
    lab_eff1 = 1
    tmp_label1 = data_format.structural_format.label_builder_format(
        label_id1, fluorophore_id1, lab_eff1, target_info1
    )

    label_id2 = "Mock_primary_antibody"
    target_info2 = dict(
        type="Primary",
        value="Mock_antibody"
    )
    fluorophore_id2 = "AF647"
    lab_eff2 = 1
    tmp_label2 = data_format.structural_format.label_builder_format(
        label_id2, fluorophore_id2, lab_eff2, target_info2
    )
    particle2, label_params_list = workflows.particle_from_structure(
        structure, [tmp_label1, tmp_label2], configuration_directory
    )
    assert len(particle2.secondary.keys()) > 0
    first_labels = particle2.emitters[target_info2["value"]].shape
    particle2.sequential_labelling=True
    particle2.generate_instance()
    second_labels = particle2.emitters[target_info2["value"]].shape
    assert first_labels[0] != second_labels[0]