from supramolsim.utils.data_format.structural_format import label_builder_format
import pytest
import numpy as np
from supramolsim import workflows

defects = np.linspace(0, 1, 27)


@pytest.fixture(scope="module", autouse=True)
def labelled_particle(configuration_directory):
    structure_id = "7R5K"
    structure_label = "7R5K_Nup96_Cterminal_direct"
    configuration_path = configuration_directory
    fluorophore_id = "AF647"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_path
    )
    labels_list = []
    labels_list.append(label_builder_format(structure_label, fluorophore_id))
    particle, label_params_list = workflows.particle_from_structure(
        structure, labels_list, configuration_path
    )
    return particle


@pytest.mark.parametrize("defect", defects)
def test_labeff_on_specific_labels(defect, labelled_particle):
    assert labelled_particle.get_ref_point().shape == (3,)
    labelled_particle.add_defects(
        eps1=300,
        xmer_neigh_distance=600,
        deg_dissasembly=defect,
    )
    assert labelled_particle.defects_target_normals is not None
