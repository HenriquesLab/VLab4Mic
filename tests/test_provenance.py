import os

import yaml
import pytest

import vlab4mic
from vlab4mic import experiments


def test_version_is_exposed():
    assert isinstance(vlab4mic.__version__, str) and vlab4mic.__version__


@pytest.mark.network
def test_parameters_recorded_with_outputs(tmp_path):
    images, _, experiment = experiments.image_vsample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        multimodal=["SMLM"], random_seed=8,
        structure_global_normal_orientation="local_plane",
    )
    params = experiment.get_parameters()
    assert params["random_seed"] == 8
    assert params["normals"]["mode"] == "local_plane"
    assert "SMLM" in params["imaging_modalities"]
    written = experiment.export_positions(str(tmp_path))
    yml = [p for p in written if p.endswith("_parameters.yml")]
    assert len(yml) == 1
    with open(yml[0]) as f:
        saved = yaml.safe_load(f)
    assert saved["normals"]["mode"] == "local_plane"
    assert saved["vlab4mic_version"] == vlab4mic.__version__


@pytest.mark.network
def test_single_modality_save_writes_files(tmp_path):
    _, _, experiment = experiments.image_vsample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        multimodal=["Widefield"], run_simulation=False, random_seed=9,
    )
    experiment.output_directory = str(tmp_path)
    experiment.run_simulation(modality="Widefield", save=True)
    files = os.listdir(tmp_path)
    assert any(f.endswith("_parameters.yml") for f in files)
    assert any(f.endswith("_emitters.csv") for f in files)
