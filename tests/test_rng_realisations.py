import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.coordinates_field import Field


def _sample(seed, **kwargs):
    sample, experiment = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        number_of_particles=4,
        random_orientations=True,
        random_rotations=True,
        clear_experiment=True,
        random_seed=seed,
        **kwargs,
    )
    return sample, experiment


@pytest.mark.network
def test_random_placement_with_minimal_distance_is_reproducible():
    # regression: positions drawn with a minimal distance ignored the seed
    s1, _ = _sample(11)
    s2, _ = _sample(11)
    np.testing.assert_array_equal(
        s1["field_emitters"]["AF647"], s2["field_emitters"]["AF647"]
    )


@pytest.mark.network
def test_rebuilt_samples_are_independent_realisations():
    # regression: every build reused the same seed, so orientations,
    # rotations and positions were identical across realisations
    _, experiment = _sample(11)
    first = np.array(experiment.coordinate_field.molecules_params["absolute_positions"])
    experiment.build(modules=["coordinate_field"], use_self_particle=True)
    second = np.array(experiment.coordinate_field.molecules_params["absolute_positions"])
    assert not np.allclose(first, second)


def test_field_draws_are_reproducible_and_not_correlated():
    def build(seed):
        field = Field()
        field.create_minimal_field(
            nmolecules=5, random_orientations=True, random_rotations=True,
            random_seed=seed,
        )
        field.generate_random_orientations()
        field.initialise_random_rotations()
        return field
    a, b = build(3), build(3)
    np.testing.assert_array_equal(
        np.array(a.molecules_params["relative_positions"]),
        np.array(b.molecules_params["relative_positions"]),
    )
    np.testing.assert_array_equal(
        a.molecules_params["rotations"], b.molecules_params["rotations"]
    )
    # positions and rotations come from one stream, not two identical ones
    x_positions = np.array(a.molecules_params["relative_positions"])[:, 0]
    assert not np.allclose(np.sort(x_positions) * 360, np.sort(a.molecules_params["rotations"]))


def test_axial_offsets_follow_seed():
    # regression for issue #135
    def build(seed):
        field = Field()
        field.create_minimal_field(
            nmolecules=6, axial_offset=list(range(0, 600, 50)), random_seed=seed
        )
        field.calculate_absolute_reference()
        field._gen_abs_from_rel_positions()
        return np.array(field.molecules_params["absolute_positions"])[:, 2]
    np.testing.assert_array_equal(build(5), build(5))
    assert not np.array_equal(build(5), build(6))
