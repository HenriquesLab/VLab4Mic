import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.coordinates_field import Field


def _field(**kwargs):
    field = Field()
    field.create_minimal_field(nmolecules=500, random_seed=2, **kwargs)
    field.generate_random_orientations()
    return np.array(field.molecules_params["orientations"])


@pytest.mark.parametrize("tilt", [5, 10, 30])
def test_tilt_mode_stays_within_cone(tilt):
    axes = _field(orientation_tilt_max=tilt)
    angles = np.degrees(np.arccos(np.clip(axes @ [0, 0, 1], -1, 1)))
    assert angles.max() <= tilt + 1e-6
    # uniform over the cap: mean cos matches (1 + cos max) / 2
    assert (axes @ [0, 0, 1]).mean() == pytest.approx(
        (1 + np.cos(np.radians(tilt))) / 2, abs=0.01
    )


def test_tilt_mode_about_set_orientation():
    axes = _field(orientation_tilt_max=10, sample_inital_orientation=[1, 0, 0])
    angles = np.degrees(np.arccos(np.clip(axes @ [1, 0, 0], -1, 1)))
    assert angles.max() <= 10 + 1e-6


@pytest.mark.network
def test_ventral_geometry_in_virtual_sample():
    _, experiment = experiments.generate_virtual_sample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        number_of_particles=6, clear_experiment=True, random_seed=5,
    )
    experiment.set_virtualsample_params(
        orientation_tilt_max=10, random_rotations=True
    )
    experiment.build(modules=["coordinate_field"], use_self_particle=True)
    axes = [m.get_axis()["direction"] for m in experiment.coordinate_field.molecules]
    for a in axes:
        a = np.asarray(a, float) / np.linalg.norm(a)
        assert np.degrees(np.arccos(np.clip(abs(a[2]), -1, 1))) <= 10 + 1e-3
