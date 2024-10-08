from supramolsim import workflows, data_format
import pytest


def test_simple_imaging_system(configuration_directory):
    configuration_path = configuration_directory
    selected_mods = [
        "STED_demo",
    ]
    fluorophore_id = "AF647"
    structure_id = "2RCJ"
    generic_label = "Generic_NHS_ester"
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
    exported_field, coordinates_field = workflows.field_from_particle(particle)
    imaging_system = workflows.create_imaging_system(
        exported_field, selected_mods, configuration_path
    )
    assert imaging_system.get_absoulte_reference_point().shape == (1, 3)
