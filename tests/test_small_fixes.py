import numpy as np
import pytest

from vlab4mic import sweep_generator
from vlab4mic.utils.transform.normals import normals_by_scaling


def test_single_site_normal_is_defined():
    # regression: one site gave a zero normal (division by zero in wobble)
    normals = normals_by_scaling(np.array([[1.0, 2.0, 3.0]]))
    np.testing.assert_allclose(normals, [[0, 0, 1]])


def test_load_reference_image_with_nothing_does_not_crash():
    gen = sweep_generator.sweep_generator()
    gen.reference_image = None
    gen.load_reference_image()
    assert gen.reference_image is None


def test_last_epitope_method():
    from vlab4mic.utils.transform.cif_builder import summarize_epitope_atoms

    atoms = [np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]], float),
             np.array([[5, 5, 5], [6, 6, 6]], float)]
    summary = summarize_epitope_atoms(atoms, method="last")
    np.testing.assert_allclose(summary, [[2, 2, 2], [6, 6, 6]])


def test_secondary_epitope_formats():
    from vlab4mic.generate import labels
    from vlab4mic.utils.io.yaml_functions import load_yaml
    import os
    from vlab4mic.experiments import ExperimentParametrisation

    config = ExperimentParametrisation().configuration_path
    base = load_yaml(os.path.join(config, "probes", "Antibody.yaml"))
    for info in ["LSPGK", {"type": "Sequence", "value": "LSPGK"},
                 {"target": {"type": "Sequence", "value": "LSPGK"}}]:
        params = dict(base, epitope_target_info=info)
        _, label_params = labels.construct_label(label_config_dictionary=params)
        assert label_params["epitope"]["target"] == {"type": "Sequence", "value": "LSPGK"}
