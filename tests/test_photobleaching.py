import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.generate.imaging import Imager


def test_constant_emission_decays_with_photobleaching():
    np.random.seed(0)
    photons = Imager._generate_constant_emission_modality(
        None, nframes=10, nemitters=20000, photons_per_second=100,
        photobleaching_rate=0.3,
    )
    mean_per_frame = photons.mean(axis=0)
    expected = 100 * np.exp(-0.3 * np.arange(10))
    np.testing.assert_allclose(mean_per_frame, expected, rtol=0.05)


def test_no_photobleaching_by_default():
    np.random.seed(1)
    photons = Imager._generate_constant_emission_modality(
        None, nframes=10, nemitters=20000, photons_per_second=100
    )
    np.testing.assert_allclose(photons.mean(axis=0), 100, rtol=0.02)


def test_blinking_model_raises_clear_error():
    imager, _ = experiments.build_virtual_microscope()
    modality = list(imager.modalities.keys())[0]
    imager.modalities[modality]["emission"] = "blinking"
    imager.fluorophore_params["f"] = {"photons_per_second": 1, "blinking": {}}
    with pytest.raises(NotImplementedError):
        imager.calculate_photons_per_frame(modality, "f", 3, 2)


@pytest.mark.network
def test_frame_series_loses_intensity():
    # reviewer 1: a ten-frame series showed no loss of intensity
    _, _, experiment = experiments.image_vsample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        multimodal=["Widefield"], run_simulation=False, random_seed=3,
    )
    experiment.set_photobleaching_rate("AF647", 0.2)
    np.random.seed(0)
    _, _, noiseless, _ = experiment.imager.generate_imaging(
        modality="Widefield", nframes=10, exp_time=0.01
    )
    totals = np.asarray(noiseless["ch0"]).sum(axis=(1, 2))
    assert totals[-1] < 0.5 * totals[0]
