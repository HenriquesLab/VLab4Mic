import os

import numpy as np
import pytest

from vlab4mic.generate import labels
from vlab4mic.utils.io.yaml_functions import load_yaml
from vlab4mic.workflows import estimate_probe_size, probe_model


def test_size_from_all_atoms():
    atoms = np.array([[0, 0, 0], [10, 0, 0], [0, 4, 0], [3, 3, 3], [5, 1, 1]], float)
    params = {"structural_atoms": {"coordinates": atoms}, "coordinates": atoms[:2]}
    assert estimate_probe_size(params) == pytest.approx(np.hypot(10, 4))


def test_size_without_model_uses_entity_points():
    entity = np.array([[0, 0, 0], [0, 0, 7], [0, 2, 7]], float)
    assert estimate_probe_size({"coordinates": entity}) == pytest.approx(np.hypot(7, 2))


@pytest.mark.network
def test_antibody_size_is_largest_dimension(configuration_directory):
    # the steric size is the largest dimension of the probe model, not the
    # spread of its conjugation sites
    params = load_yaml(os.path.join(configuration_directory, "probes", "Antibody.yaml"))
    label_object, _ = labels.construct_label(label_config_dictionary=params)
    probe, emitters, *_ = probe_model(
        model=label_object.model, binding=label_object.binding,
        conjugation_sites=label_object.conjugation, epitope=label_object.epitope,
        config_dir=configuration_directory,
    )
    size = estimate_probe_size({"structural_atoms": probe.assembly_atoms,
                                "coordinates": emitters})
    # 1HZH spans about 17 nm (coordinates in Angstrom)
    assert 160 < size < 180
