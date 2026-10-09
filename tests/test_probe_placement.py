import numpy as np
import pytest

from vlab4mic.utils.transform.points_transforms import (
    decorate_epitopes_normals,
    labeling_reorient_set,
)

# labelling entity: pivot (anchor), axis point, then 3 emitters
ENTITY = np.array(
    [[0, 0, 0], [0, 0, 1], [0, 0, 5], [1, 0, 5], [0, 1, 5]], dtype=float
)


@pytest.mark.parametrize(
    "direction", [[0, 0, -1], [0, 0, 1], [1, 0, 0], [0, 1e-9, -1]]
)
def test_probe_axis_follows_direction_including_antiparallel(direction):
    # regression: an axis exactly opposite the target direction was not flipped
    direction = np.array(direction, dtype=float)
    points, _ = labeling_reorient_set(ENTITY, 0, 1, direction, np.zeros(3))
    axis = points[1] - points[0]
    np.testing.assert_allclose(axis / np.linalg.norm(axis),
                               direction / np.linalg.norm(direction), atol=1e-6)


def _decorate(dol, n=3):
    normals = np.tile([0, 0, 1.0], (n, 1))
    epitopes = np.arange(3 * n, dtype=float).reshape(n, 3) * 10
    return decorate_epitopes_normals([normals, epitopes], ENTITY, dol=dol)


def test_zero_dol_gives_no_emitters_without_crash():
    # regression: IndexError when every probe drew 0 fluorophores
    emitters, normals = _decorate(dol=0)
    assert emitters.shape == (0, 3)
    assert normals == []


def test_normals_follow_probes_with_dol():
    # regression: normals were not returned when DoL was set
    np.random.seed(0)
    emitters, normals = _decorate(dol=100)
    assert len(normals) == 3
    assert emitters.shape == (9, 3)


def test_no_dol_uses_every_site():
    emitters, normals = _decorate(dol=None)
    assert emitters.shape == (9, 3)
    assert len(normals) == 3
