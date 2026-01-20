from vlab4mic.utils.data_format.structural_format import label_builder_format
import pytest
import numpy as np
from vlab4mic import workflows

structural_integrity_values = np.linspace(0, 1, 3)


@pytest.mark.parametrize("structural_integrity", structural_integrity_values)
def test_add_structural_integrity(structural_integrity, experiment_7r5k_base):
    assert experiment_7r5k_base.generators_status("particle")
    experiment_7r5k_base.particle.add_structural_integrity(
        eps1=300,
        xmer_neigh_distance=600,
        integrity=structural_integrity,
    )
    if structural_integrity == 1:
        assert experiment_7r5k_base.particle.structural_integrity_target_normals is None
    else:
        assert experiment_7r5k_base.particle.structural_integrity_target_normals is not None
