import supramolsim
from supramolsim import workflows, data_format
from supramolsim.experiments import create_experiment_parametrisation
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
    selected_mods = dict(
        STED=None,
        Confocal=None,
    )
    structure_and_labels = dict(
        structure_id="7R5K",
        structure_label="NPC_Nup96_Cterminal_direct",
        fluorophore_id="AF647",
    )
    defects_eps_d = dict(eps1=300, eps2=600)
    savging = dict(
        experiment_id="SupraMolSim_experiment", output_directory="", save=True
    )
    Experiment_generator = create_experiment_parametrisation(
        structure_and_labels=structure_and_labels,
        modalities_acquisition=selected_mods,
        savging=savging,
        defects_params=defects_eps_d,
        use_locals=True,
    )
    return Experiment_generator
