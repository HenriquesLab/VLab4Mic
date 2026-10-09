import os

import numpy as np
import pytest

from vlab4mic import sweep_generator
from vlab4mic.sweep_generator import values_from_range


@pytest.mark.parametrize(
    "rng, expected",
    [
        ((0, 0.3, 0.1), [0, 0.1, 0.2, 0.3]),
        ((0, 0.7, 0.1), [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]),
        ((0.1, 0.7, 0.2), [0.1, 0.3, 0.5, 0.7]),
        ((0, 1, 0.25), [0, 0.25, 0.5, 0.75, 1]),
        ((80, 120, 20), [80, 100, 120]),
    ],
)
def test_values_from_range_includes_stop(rng, expected):
    # regression: (0, 0.3, 0.1) gave 3 values with a changed step
    np.testing.assert_allclose(values_from_range(*rng), expected)


def test_integer_ranges_stay_integer():
    assert values_from_range(80, 120, 20) == [80, 100, 120]
    assert all(isinstance(v, int) for v in values_from_range(1, 5, 2))


def test_tuples_apply_to_parameters_without_settings():
    # regression: tuples were ignored for parameters not listed in
    # parameter_settings.yaml (all optical parameters)
    gen = sweep_generator.sweep_generator()
    gen.set_sweep_parameters(pixelsize_nm=(80, 120, 20))
    assert gen.params_by_group["modality"]["pixelsize_nm"] == [80, 100, 120]


@pytest.fixture(scope="module")
def small_sweep(tmp_path_factory):
    img_mask = np.zeros((8, 8))
    img_mask[2, 2] = img_mask[5, 5] = 1
    return sweep_generator.run_parameter_sweep(
        sweep_repetitions=1,
        exp_time=[0.001, 0.01],
        modalities=["Widefield"],
        image4vsample=img_mask,
        image4vsample_parameters=dict(pixelsize=100, mode="mask", npositions=2),
        reference_image=img_mask,
        reference_image_parameters={"ref_pixelsize": 100},
        output_name="test",
        return_generator=True,
        analysis_plots=False,
        save_sweep_images=False,
        save_analysis_results=False,
        run_analysis=True,
        random_seed=1,
    )


def test_swept_exposure_time_is_used(small_sweep):
    # regression: swept exp_time values were replaced by {"channels": [...]}
    used = set()
    for parameters in small_sweep.acquisition_outputs_parameters.values():
        for p in parameters:
            if isinstance(p, dict) and "exp_time" in p:
                used.add(float(p["exp_time"]))
    assert used == {0.001, 0.01}


def test_reference_saved_without_overwriting(small_sweep, tmp_path):
    # regression: the reference was saved under the last condition's name
    small_sweep.save_images(output_directory=str(tmp_path))
    files = sorted(os.listdir(tmp_path))
    assert "reference.tiff" in files
    n_conditions = len(small_sweep.acquisition_outputs)
    assert sum(f.endswith(".tiff") for f in files) == n_conditions + 1
