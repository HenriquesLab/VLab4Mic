import supramolsim
from supramolsim import workflows
from supramolsim.utils import data_format
from supramolsim import experiments, sweep_generator
import pytest
import os


@pytest.fixture(scope="module")
def configuration_directory():
    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    conf_dif = os.path.join(pck_dir, "configs")
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
    return exported_field


@pytest.fixture(scope="module")
def experiment_7r5k_base():
    structure_id = "7R5K"
    probe = "NPC_Nup96_Cterminal_direct"
    fluorophore_id = "AF647"
    virtual_sample = "square1x1um_randomised" 
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    imaging_output7r5k, exp7r5k = experiments.image_vsample(
        structure=structure_id,
        probe_name=probe,
        virtual_sample_template=virtual_sample,
        multimodal=modalities,
        run_simulation=False 
    )
    return exp7r5k


@pytest.fixture(scope="module")
def sweep_gen():
    sweep_object = sweep_generator.sweep_generator()
    return sweep_object