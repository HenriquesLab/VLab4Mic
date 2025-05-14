import supramolsim
from supramolsim import workflows
from supramolsim.utils import data_format
from supramolsim import experiments 
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
    modalities = ["Widefield", "Confocal", "SMLM", "STED", "AiryScan"]
    selected_mods = {}
    default_aqc = dict(
        nframes=2,
        exp_time=0.005
    )
    for mod in modalities:
        selected_mods[mod] = default_aqc
    myexperiment = experiments.ExperimentParametrisation()
    myexperiment.structure_id = structure_id
    myexperiment.structure_label = probe
    myexperiment.fluorophore_id = fluorophore_id
    myexperiment.coordinate_field_id = virtual_sample
    myexperiment.selected_mods = selected_mods
    myexperiment.build(use_locals=True)
    return myexperiment
