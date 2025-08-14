from vlab4mic import workflows
import pytest
import copy

structure_list = [
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


def test_structure_normals(experiment_7r5k_base):
    structure = copy.deepcopy(experiment_7r5k_base.structure)
    structure.assign_normals2targets()
    structure.show_target_labels(with_normals=True, show_axis=True)
