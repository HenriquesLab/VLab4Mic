from supramolsim import workflows
from supramolsim.utils import data_format
import pytest


def test_simple_imaging_system(configuration_directory):
    configuration_path = configuration_directory
    selected_mods = [
        "STED",
    ]
    fluorophore_id = "AF647"
    structure_id = "2RCJ"
    generic_label = "NHS_ester"
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
    exported_field, coordinates_field = workflows.field_from_particle(particle)
    imaging_system, modality_parameters = workflows.create_imaging_system(
        exported_field, selected_mods, configuration_path
    )
    assert imaging_system.get_absoulte_reference_point().shape == (1, 3)


def test_multi_imaging_system(configuration_directory):
    configuration_path = configuration_directory
    selected_mods = ["STED", "Confocal"]
    fluorophore_id = "AF647"
    structure_id = "2RCJ"
    generic_label = "NHS_ester"
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
    exported_field, coordinates_field = workflows.field_from_particle(particle)
    imaging_system, modality_parameters = workflows.create_imaging_system(
        exported_field, selected_mods, configuration_path
    )
    assert imaging_system.get_absoulte_reference_point().shape == (1, 3)


def test_image_from_field(configuration_directory, gt_structural_model_field):
    configuration_path = configuration_directory
    selected_mods = [
        "STED",
    ]
    imaging_system, modality_parameters = workflows.create_imaging_system(
        gt_structural_model_field, selected_mods, configuration_path
    )
    assert imaging_system.get_absoulte_reference_point().shape == (1, 3)
