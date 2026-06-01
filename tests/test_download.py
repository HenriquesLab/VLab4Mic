import os
from pathlib import Path

import pytest

import vlab4mic
from vlab4mic import download


@pytest.fixture
def configs_dir():
    """Installed package configs dir (mirrors conftest's configuration_directory,
    but defined locally so this module needs no heavy imports to collect)."""
    return os.path.join(os.path.dirname(vlab4mic.__file__), "configs")


def test_default_path_resolves_to_package_configs():
    """download_suggested_structures() must default to the installed package
    configs dir, which actually contains structures/*.yaml."""
    expected = os.path.join(os.path.dirname(vlab4mic.__file__), "configs")
    structures = Path(expected) / "structures"
    assert structures.is_dir()
    assert list(structures.glob("*.yaml"))


def test_structure_yaml_keys_are_nested_under_model(configs_dir):
    """Regression: metadata lives under the `model:` mapping with a capital
    `ID`, so download must read data["model"]["ID"]/["title"], not data["id"]/
    data["title"]."""
    import yaml

    structures = Path(configs_dir) / "structures"
    yaml_files = [
        f for f in structures.glob("*.yaml") if not f.stem.startswith("_template")
    ]
    assert yaml_files, "no real structure yaml files found"

    for f in yaml_files:
        with open(f, "r") as fh:
            data = yaml.load(fh, Loader=yaml.FullLoader)
        # The keys the fixed code relies on must exist (ID drives the download;
        # title is display-only and may be blank for some structures).
        assert data["model"]["ID"]
        assert "title" in data["model"]
        # ...and the old flat keys must NOT (proves the original code was broken).
        assert "id" not in data
        assert "title" not in data


def test_download_skips_existing_without_network(
    configs_dir, monkeypatch, tmp_path
):
    """End-to-end-ish: with every .cif already present, the function walks all
    real structure yamls, parses their metadata, and never hits the network."""
    src = Path(configs_dir) / "structures"
    data_root = tmp_path / "configs"
    structures = data_root / "structures"
    structures.mkdir(parents=True)

    real_ids = []
    for yaml_file in src.glob("*.yaml"):
        (structures / yaml_file.name).write_text(yaml_file.read_text())
        if not yaml_file.stem.startswith("_template"):
            # pre-create the .cif so the function skips downloading
            (structures / yaml_file.with_suffix(".cif").name).write_text("stub")
            real_ids.append(yaml_file.stem)

    assert real_ids, "expected at least one real structure"

    def _boom(*args, **kwargs):
        raise AssertionError("network download attempted for an existing .cif")

    monkeypatch.setattr(download, "download_file", _boom)

    # Must not raise (no KeyError on parsing, no network call on skip).
    download.download_suggested_structures(data_path=str(data_root))
