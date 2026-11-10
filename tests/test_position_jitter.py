"""Random particle placement (position_jitter_nm)."""
import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.coordinates_field import Field


def _field_positions(seed, **kwargs):
    field = Field()
    field.create_minimal_field(random_seed=seed, **kwargs)
    return np.array(field.get_molecule_param("absolute_positions")).reshape(-1, 3)


def test_no_jitter_by_default():
    np.testing.assert_array_equal(_field_positions(1), _field_positions(2))


def test_jitter_offsets_x_and_y_within_amplitude():
    centre = _field_positions(1)
    offsets = np.vstack([_field_positions(seed, position_jitter_nm=50) - centre
                         for seed in range(200)])
    # absolute positions are in nm (field scale 1e-9 m)
    assert np.all(np.abs(offsets[:, :2]) <= 50)
    assert offsets[:, :2].std(axis=0).min() > 20  # uniform on +-50: s.d. ~29
    np.testing.assert_array_equal(offsets[:, 2], 0)  # z unchanged


def test_jitter_per_axis_and_reproducible():
    centre = _field_positions(1)
    offsets = np.vstack([_field_positions(seed, position_jitter_nm=[10, 20, 30]) - centre
                         for seed in range(200)])
    assert np.all(np.abs(offsets) <= [10, 20, 30])
    assert offsets[:, 2].std() > 10
    np.testing.assert_array_equal(_field_positions(7, position_jitter_nm=50),
                                  _field_positions(7, position_jitter_nm=50))


def test_jitter_needs_one_or_three_values():
    with pytest.raises(ValueError):
        _field_positions(1, position_jitter_nm=[10, 20])


def _centroids(results):
    return np.array([r["field_emitters"]["AF647"].mean(axis=0) for r in results])


@pytest.mark.network
def test_replicates_are_placed_at_different_positions():
    def replicates(**kwargs):
        _, _, experiment = experiments.image_vsample(
            structure="1XI5", probe_template="Antibody", probe_target_type="Sequence",
            probe_target_value="EQATETQ", multimodal=["STED"], run_simulation=False,
            clear_experiment=True, random_seed=3, **kwargs,
        )
        return experiment.run_replicates(6, modality="STED")

    fixed = _centroids(replicates())
    jittered = _centroids(replicates(position_jitter_nm=100))
    assert fixed[:, :2].std(axis=0).max() < 5  # same position (labelling varies the centroid)
    assert jittered[:, :2].std(axis=0).min() > 20


@pytest.mark.network
def test_sweep_repetitions_are_placed_at_different_positions():
    from vlab4mic.sweep_generator import run_parameter_sweep

    sweep = run_parameter_sweep(
        structures=["1XI5"], probe_templates=["Antibody"], probe_target_type=["Sequence"],
        probe_target_value=["EQATETQ"], modalities=["STED"], sweep_repetitions=6,
        labelling_efficiency=[1.0], position_jitter_nm=100, use_reference=False,
        run_analysis=False, return_generator=True, save_sweep_images=False,
        save_analysis_results=False, analysis_plots=False, random_seed=3,
    )
    sweep.generate_virtual_samples()
    (realisations,) = sweep.virtual_samples.values()
    centroids = _centroids(realisations)
    assert centroids[:, :2].std(axis=0).min() > 20
    # offsets do not accumulate across repetitions
    assert np.abs(centroids[:, :2] - centroids[:, :2].mean(axis=0)).max() < 150
