import vlab4mic
from vlab4mic import workflows
from vlab4mic.utils import data_format
from vlab4mic import experiments, sweep_generator
import pytest
import os


def pytest_addoption(parser):
    parser.addoption(
        "--no-network",
        action="store_true",
        default=False,
        help="Skip all tests marked with @pytest.mark.network (for offline / sandboxed environments).",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--no-network"):
        skip = pytest.mark.skip(
            reason="requires network access to download .cif files (--no-network)"
        )
        for item in items:
            if "network" in item.keywords:
                item.add_marker(skip)



@pytest.fixture(scope="module")
def configuration_directory():
    pck_dir = os.path.dirname(os.path.abspath(vlab4mic.__file__))
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
    coordinate_exp = experiments.ExperimentParametrisation()
    coordinate_exp.exported_coordinate_field
    return coordinate_exp.exported_coordinate_field


@pytest.fixture(scope="module")
def experiment_7r5k_base():
    structure_id = "7R5K"
    probe = "NPC_Nup96_Cterminal_direct"
    fluorophore_id = "AF647"
    virtual_sample = "square1x1um_randomised"
    modalities = ["Widefield", "Confocal", "SMLM", "STED"]
    imaging_output7r5k, images_noiseless, exp7r5k = experiments.image_vsample(
        structure=structure_id,
        probe_template=probe,
        virtual_sample_template=virtual_sample,
        multimodal=modalities,
        run_simulation=False,
        clear_experiment=True,
    )
    return exp7r5k


@pytest.fixture(scope="module")
def sweep_gen():
    sweep_object = sweep_generator.sweep_generator()
    return sweep_object
