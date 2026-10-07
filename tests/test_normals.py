import numpy as np
import pytest

from vlab4mic.utils.transform.normals import (
    normals_by_local_plane,
    normals_by_scaling,
)


def flat_lattice(n_side=15, spacing=10.0, z=50.0):
    g = np.stack(
        np.meshgrid(np.arange(n_side) * spacing, np.arange(n_side) * spacing), -1
    ).reshape(-1, 2)
    return np.c_[g, np.full(len(g), z)]


def fibonacci_sphere(n=400, radius=50.0):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return radius * np.c_[
        np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    ]


def cylinder(radius=12.5, length=200.0, n_around=26, n_along=40):
    angles = np.linspace(0, 2 * np.pi, n_around, endpoint=False)
    zs = np.linspace(-length / 2, length / 2, n_along)
    a, z = np.meshgrid(angles, zs)
    return np.c_[radius * np.cos(a.ravel()), radius * np.sin(a.ravel()), z.ravel()]


def unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def test_local_plane_normals_on_flat_lattice_follow_reference():
    points = flat_lattice()
    up = normals_by_local_plane(points, reference_vector=np.array([0, 0, 1]))
    np.testing.assert_allclose(up, np.tile([0, 0, 1.0], (len(points), 1)), atol=1e-8)
    down = normals_by_local_plane(points, reference_vector=np.array([0, 0, -1]))
    np.testing.assert_allclose(down, -up, atol=1e-8)


def test_scaling_normals_fail_on_flat_lattice():
    # documents the limitation that local_plane addresses: on a flat surface
    # the scaling normals lie in the plane instead of along its normal
    points = flat_lattice()
    normals = unit(normals_by_scaling(points)[np.any(points[:, :2] != 70, axis=1)])
    assert np.allclose(normals[:, 2], 0)


def test_local_plane_normals_on_sphere_point_outward():
    points = fibonacci_sphere()
    normals = normals_by_local_plane(points)
    radial = unit(points)
    cos = np.einsum("ni,ni->n", normals, radial)
    assert np.all(cos > 0.98)


def test_local_plane_normals_on_cylinder_are_radial():
    points = cylinder()
    normals = normals_by_local_plane(points, n_neighbours=8)
    radial = unit(np.c_[points[:, :2], np.zeros(len(points))])
    cos = np.einsum("ni,ni->n", normals, radial)
    # sites away from the two open ends
    body = np.abs(points[:, 2]) < 80
    assert np.all(cos[body] > 0.98)


def test_local_plane_normals_are_unit_length():
    normals = normals_by_local_plane(fibonacci_sphere())
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1)


@pytest.mark.parametrize("n_points", [0, 1, 2])
def test_local_plane_normals_with_too_few_points(n_points):
    points = fibonacci_sphere()[:n_points]
    normals = normals_by_local_plane(points, reference_vector=np.array([0, 0, 2]))
    assert normals.shape == (n_points, 3)
    np.testing.assert_allclose(normals, np.tile([0, 0, 1.0], (n_points, 1)))


class _StructureStub:
    """Minimal object to call MolecularStructure normals methods on."""

    def __init__(self, mode, normal_vector=None):
        self.label_targets = {
            "a": {"coordinates": flat_lattice(), "normals": None},
            "b": {"coordinates": flat_lattice(z=0), "normals": None},
        }
        self.normals_params = {"mode": mode, "normal_vector": normal_vector}
        self.axis = {"pivot": np.zeros(3), "direction": np.array([0, 0, 1])}


@pytest.mark.parametrize("mode", ["scaling", "local_plane", "global", "structure_axis"])
def test_assign_normals_to_single_target(mode):
    # regression: the "global" and "structure_axis" branches used an
    # undefined name when a single target was given
    from vlab4mic.generate.molecular_structure import MolecularReplicates as MolecularStructure

    stub = _StructureStub(mode, normal_vector=np.array([0, 0, 1]))
    stub._compute_target_normals = (
        lambda coords: MolecularStructure._compute_target_normals(stub, coords)
    )
    MolecularStructure.assign_normals2targets(stub, target="a")
    assert stub.label_targets["a"]["normals"].shape == (225, 3)
    assert stub.label_targets["b"]["normals"] is None


def test_local_plane_mode_uses_structure_axis_on_flat_targets():
    from vlab4mic.generate.molecular_structure import MolecularReplicates as MolecularStructure

    stub = _StructureStub("local_plane")
    stub.axis["direction"] = np.array([0, 0, -1])
    stub._compute_target_normals = (
        lambda coords: MolecularStructure._compute_target_normals(stub, coords)
    )
    MolecularStructure.assign_normals2targets(stub)
    for target in stub.label_targets.values():
        np.testing.assert_allclose(target["normals"][:, 2], -1, atol=1e-8)
