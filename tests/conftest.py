import supramolsim
from supramolsim import workflows, data_format
import pytest
import os


@pytest.fixture(scope="module")
def configuration_directory():
    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    conf_dif = os.path.join(pck_dir, "configuration")
    return conf_dif


@pytest.fixture(scope="module")
def structure_test(configuration_directory):
    structure_id = "7R5K"
    configuration_path = configuration_directory
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_path
    )
    return structure


@pytest.fixture(scope="module")
def gt_structural_model_field(configuration_directory):
    configuration_path = configuration_directory
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
    return exported_field
