import os

import numpy as np
import pytest

from vlab4mic import experiments
from vlab4mic.utils.io.localisation_table import (
    COLUMNS,
    read_localisation_table,
    write_localisation_table,
)


def test_table_round_trip(tmp_path):
    coords = np.array([[1.5, 2.25, 3.0], [10, 20, 30]])
    path = write_localisation_table(
        str(tmp_path / "t.csv"), coords, photons=[100, 200], uncertainty=5
    )
    with open(path) as f:
        header = f.readline().strip()
    assert header == ",".join(f'"{c}"' for c in COLUMNS)
    table = read_localisation_table(path)
    np.testing.assert_allclose(table["x [nm]"], [1.5, 10])
    np.testing.assert_allclose(table["intensity [photon]"], [100, 200])
    np.testing.assert_allclose(table["uncertainty [nm]"], [5, 5])
    np.testing.assert_allclose(table["frame"], [1, 1])


def test_emitter_table_has_empty_uncertainty(tmp_path):
    path = write_localisation_table(str(tmp_path / "e.csv"), np.zeros((2, 3)))
    assert np.isnan(read_localisation_table(path)["uncertainty [nm]"]).all()


@pytest.mark.network
def test_export_emitters_and_localisations(tmp_path):
    images, _, experiment = experiments.image_vsample(
        structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
        multimodal=["SMLM"], random_seed=4,
    )
    positions = experiment.imager.get_positions("SMLM")
    record = positions["ch0"][list(positions["ch0"])[0]]
    assert record["emitters"].shape[1] == 3 and len(record["emitters"]) > 0
    assert record["localisations"] is not None
    written = experiment.export_positions(str(tmp_path))
    names = sorted(os.path.basename(p) for p in written)
    assert any(n.endswith("_emitters.csv") for n in names)
    assert any(n.endswith("_localisations.csv") for n in names)
    emitters = read_localisation_table(
        [p for p in written if p.endswith("_emitters.csv")][0]
    )
    np.testing.assert_allclose(emitters["x [nm]"], record["emitters"][:, 0], atol=1e-3)
