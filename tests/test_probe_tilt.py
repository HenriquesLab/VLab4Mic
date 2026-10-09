import os

import numpy as np
import pytest

from vlab4mic.generate import labels
from vlab4mic.utils.io.yaml_functions import load_yaml
from vlab4mic.utils.transform.cif_builder import indirect_labelling
from vlab4mic.utils.transform.normals import tilt_direction


@pytest.mark.parametrize("tilt", [0, 15, 30, 60])
def test_tilt_angle_is_exact(tilt):
    rng = np.random.default_rng(0)
    normal = np.array([0.0, 0.0, 2.0])
    for _ in range(50):
        d = tilt_direction(normal, tilt, rng=rng)
        assert np.linalg.norm(d) == pytest.approx(1)
        assert np.degrees(np.arccos(np.clip(d[2], -1, 1))) == pytest.approx(tilt, abs=1e-6)


def test_tilt_azimuth_is_uniform():
    rng = np.random.default_rng(1)
    dirs = np.array([tilt_direction([0, 0, 1], 30, rng=rng) for _ in range(4000)])
    # the mean in-plane component vanishes for a uniform azimuth
    assert np.abs(dirs[:, :2].mean(axis=0)).max() < 0.02


def test_template_and_override(configuration_directory):
    params = load_yaml(os.path.join(configuration_directory, "probes", "Antibody.yaml"))
    _, label_params = labels.construct_label(label_config_dictionary=params)
    assert label_params["binding"]["tilt"] == 0
    params["tilt_theta"] = 20
    _, label_params = labels.construct_label(label_config_dictionary=params)
    assert label_params["binding"]["tilt"] == 20


def test_probes_are_placed_at_the_tilt_angle():
    # labelling entity: anchor, axis point, one emitter 5 units along the axis
    entity = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 5]], dtype=float)
    epitopes = np.array([[10.0 * i, 0, 0] for i in range(20)])
    normals = np.tile([0, 0, 1.0], (20, 1))
    label_data = {
        "emitters": entity,
        "labelling_efficiency": 1,
        "minimal_distance": 0,
        "conjugation_sites": {"DoL": None},
        "binding": {"tilt": 25, "wobble_range": {"theta": None}},
    }
    np.random.seed(0)
    emitters, n, _ = indirect_labelling(
        {"coordinates": epitopes, "normals": normals}, label_data
    )
    # each emitter sits 5 units from its epitope at 25 degrees from the normal
    heights = emitters[:, 2]
    np.testing.assert_allclose(heights, 5 * np.cos(np.radians(25)), atol=1e-6)
