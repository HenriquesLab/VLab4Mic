from supramolsim.utils.data_format.structural_format import label_builder_format
import pytest
import numpy as np
from supramolsim import workflows

defects = np.linspace(0, 1, 3)

@pytest.mark.parametrize("defect", defects)
def test_add_defects(defect, experiment_7r5k_base):
    assert experiment_7r5k_base.generators_status("particle")
    experiment_7r5k_base.particle.add_defects(
        eps1=300,
        xmer_neigh_distance=600,
        deg_dissasembly=defect,
    )
    if defect == 0:
        assert experiment_7r5k_base.particle.defects_target_normals is None
    else:
        assert experiment_7r5k_base.particle.defects_target_normals is not None
