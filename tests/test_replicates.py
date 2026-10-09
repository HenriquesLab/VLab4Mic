import os

import numpy as np
import pytest

from vlab4mic import experiments

KWARGS = dict(
    structure="7R5K",
    probe_template="NPC_Nup96_Cterminal_direct",
    labelling_efficiency=0.5,
    multimodal=["Widefield"],
    number_of_particles=2,
)


def _emitters(result):
    return np.concatenate(list(result["field_emitters"].values()))


@pytest.mark.network
def test_replicates_are_independent_and_reproducible(tmp_path):
    results, _ = experiments.run_replicates(
        n_replicates=3, random_seed=7, export_directory=str(tmp_path), **KWARGS
    )
    assert [r["replicate"] for r in results] == [0, 1, 2]
    # each realisation has its own labelling and placement
    assert not np.array_equal(_emitters(results[0]), _emitters(results[1]))
    assert "Widefield" in results[0]["images"]
    record = results[0]["positions"]["Widefield"]["ch0"]
    assert len(list(record.values())[0]["emitters"]) > 0
    files = os.listdir(tmp_path)
    assert any("rep000" in f and f.endswith("_emitters.csv") for f in files)

    again, _ = experiments.run_replicates(n_replicates=3, random_seed=7, **KWARGS)
    for a, b in zip(results, again):
        np.testing.assert_array_equal(_emitters(a), _emitters(b))
