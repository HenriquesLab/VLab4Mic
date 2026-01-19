from vlab4mic.utils.data_format.structural_format import label_builder_format
import pytest
import numpy as np
from vlab4mic import workflows

incomplete_labelling_values = np.linspace(0, 1, 3)


@pytest.mark.parametrize("incomplete_labelling", incomplete_labelling_values)
def test_add_incomplete_labelling(incomplete_labelling, experiment_7r5k_base):
    assert experiment_7r5k_base.generators_status("particle")
    experiment_7r5k_base.particle.add_incomplete_labelling(
        eps1=300,
        xmer_neigh_distance=600,
        deg_dissasembly=incomplete_labelling,
    )
    if incomplete_labelling == 0:
        assert experiment_7r5k_base.particle.incomplete_labelling_target_normals is None
    else:
        assert experiment_7r5k_base.particle.incomplete_labelling_target_normals is not None
