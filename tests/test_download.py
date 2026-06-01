import os
from pathlib import Path

import pytest

import vlab4mic
from vlab4mic import workflows
from vlab4mic import download


@pytest.fixture
def configs_dir():
    """Installed package configs dir (mirrors conftest's configuration_directory,
    but defined locally so this module needs no heavy imports to collect)."""
    return os.path.join(os.path.dirname(vlab4mic.__file__), "configs")


def test_default_yaml_source_resolves_to_package_configs():
    """The default YAML source remains the installed package configs dir."""
    expected = os.path.join(os.path.dirname(vlab4mic.__file__), "configs")
    structures = Path(expected) / "structures"
    assert structures.is_dir()
    assert list(structures.glob("*.yaml"))


def test_default_structure_download_dir_uses_user_storage(monkeypatch, tmp_path):
    monkeypatch.delenv("VLAB4MIC_STRUCTURE_DIR", raising=False)
    monkeypatch.setattr(download.Path, "home", lambda: tmp_path)

    assert download.get_structure_download_dir() == (
        tmp_path / ".vlab4mic" / "structures"
    )


def test_structure_download_dir_can_be_overridden(monkeypatch, tmp_path):
    custom_dir = tmp_path / "custom-structures"
    monkeypatch.setenv("VLAB4MIC_STRUCTURE_DIR", str(custom_dir))

    assert download.get_structure_download_dir() == custom_dir


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


def test_download_suggested_structures_writes_to_user_cache(
    configs_dir, monkeypatch, tmp_path
):
    """YAML configs are read from data_path, but CIFs are written to cache."""
    src = Path(configs_dir) / "structures"
    data_root = tmp_path / "configs"
    structures = data_root / "structures"
    structures.mkdir(parents=True)
    cache_dir = tmp_path / "cache" / "structures"

    real_ids = []
    for yaml_file in src.glob("*.yaml"):
        (structures / yaml_file.name).write_text(yaml_file.read_text())
        if not yaml_file.stem.startswith("_template"):
            real_ids.append(yaml_file.stem)

    assert real_ids, "expected at least one real structure"
    downloaded = []

    def _fake_download(url, fname, **kwargs):
        downloaded.append((url, Path(fname)))
        Path(fname).write_text("stub")

    monkeypatch.setattr(download, "download_file", _fake_download)

    download.download_suggested_structures(
        data_path=str(data_root), download_dir=str(cache_dir)
    )

    assert downloaded
    assert all(fname.parent == cache_dir for _url, fname in downloaded)
    assert list(cache_dir.glob("*.cif"))
    assert not list(structures.glob("*.cif"))


def test_download_skips_existing_cache_without_network(
    configs_dir, monkeypatch, tmp_path
):
    """With every .cif already cached, metadata parsing does not hit network."""
    src = Path(configs_dir) / "structures"
    data_root = tmp_path / "configs"
    structures = data_root / "structures"
    structures.mkdir(parents=True)
    cache_dir = tmp_path / "cache" / "structures"
    cache_dir.mkdir(parents=True)

    real_ids = []
    for yaml_file in src.glob("*.yaml"):
        (structures / yaml_file.name).write_text(yaml_file.read_text())
        if not yaml_file.stem.startswith("_template"):
            (cache_dir / yaml_file.with_suffix(".cif").name).write_text("stub")
            real_ids.append(yaml_file.stem)

    assert real_ids, "expected at least one real structure"

    def _boom(*args, **kwargs):
        raise AssertionError("network download attempted for an existing .cif")

    monkeypatch.setattr(download, "download_file", _boom)

    download.download_suggested_structures(
        data_path=str(data_root), download_dir=str(cache_dir)
    )


def test_verify_structure_uses_user_cache(monkeypatch, tmp_path):
    cache_dir = tmp_path / "vlab-cache" / "structures"
    monkeypatch.setenv("VLAB4MIC_STRUCTURE_DIR", str(cache_dir))

    def _fake_download(url, fname, **kwargs):
        Path(fname).write_text("stub")

    monkeypatch.setattr(download, "download_file", _fake_download)

    cif_path = download.verify_structure("9I0K")

    assert cif_path == cache_dir / "9I0K.cif"
    assert cif_path.read_text() == "stub"


def test_verify_structure_keeps_explicit_directory(monkeypatch, tmp_path):
    explicit_dir = tmp_path / "explicit"

    def _fake_download(url, fname, **kwargs):
        Path(fname).write_text("stub")

    monkeypatch.setattr(download, "download_file", _fake_download)

    cif_path = download.verify_structure("1XI5", structure_dir=str(explicit_dir))

    assert cif_path == explicit_dir / "1XI5.cif"
    assert cif_path.exists()


def test_load_structure_downloads_to_default_cache(configs_dir, monkeypatch):
    calls = []

    def _fake_verify_structure(structure_id):
        calls.append(structure_id)
        return Path("/tmp/vlab4mic-cache") / f"{structure_id}.cif"

    def _fake_build_structure_cif(**kwargs):
        return object()

    monkeypatch.setattr(workflows, "verify_structure", _fake_verify_structure)
    monkeypatch.setattr(workflows, "build_structure_cif", _fake_build_structure_cif)

    structure, structure_params = workflows.load_structure("1XI5", configs_dir)

    assert structure is not None
    assert structure_params["model"]["ID"] == "1XI5"
    assert calls == ["1XI5"]


def test_load_structure_from_local_file_does_not_verify(configs_dir, monkeypatch):
    def _boom(*args, **kwargs):
        raise AssertionError("local structure paths must not call verify_structure")

    build_calls = []

    def _fake_build_structure_cif(**kwargs):
        build_calls.append(kwargs)
        return object()

    monkeypatch.setattr(workflows, "verify_structure", _boom)
    monkeypatch.setattr(workflows, "build_structure_cif", _fake_build_structure_cif)

    structure, structure_params = workflows.load_structure(
        "LOCAL",
        configs_dir,
        structure_path="/tmp/local-model.pdb",
        structure_format="PDB",
    )

    assert structure is not None
    assert structure_params["model"]["format"] == "PDB"
    assert build_calls[0]["cif_file"] == "/tmp/local-model.pdb"
