from supramolsim import workflows
import pytest

structure_list = [
    "7R5K",
    "2RCJ",
]


@pytest.mark.parametrize("structure_id", structure_list)
def test_load_structure(structure_id, configuration_directory):
    configuration_path = configuration_directory
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_path
    )
    assert structure.assembly_refpt.shape == (3,)


# test user-input file
