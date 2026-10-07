import numpy as np
import pytest

from vlab4mic import experiments


@pytest.fixture(scope="module")
def smlm_experiment():
    _, _, experiment = experiments.image_vsample(
        multimodal=["SMLM"], run_simulation=False, random_seed=3
    )
    return experiment


def _smlm_total_intensity(experiment, exp_time, photons_per_second):
    experiment.imager.set_fluorophore_photons_per_second(
        "AF647", photons_per_second
    )
    np.random.seed(3)
    _, _, noiseless, _ = experiment.imager.generate_imaging(
        modality="SMLM", exp_time=exp_time
    )
    return np.asarray(noiseless["ch0"]).sum()


@pytest.mark.network
@pytest.mark.parametrize(
    "exp_time, photons_per_second",
    [(0.01, 100000), (0.001, 1000000)],
)
def test_smlm_intensity_follows_photon_budget(
    smlm_experiment, exp_time, photons_per_second
):
    # regression: SMLM localisations used a fixed 100 photons, ignoring
    # fluorophore brightness and exposure time
    base = _smlm_total_intensity(smlm_experiment, 0.001, 100000)
    scaled = _smlm_total_intensity(smlm_experiment, exp_time, photons_per_second)
    assert scaled / base == pytest.approx(10, rel=0.01)
