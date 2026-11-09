"""
The faster implementations give the same output as the ones they replace.
"""
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from vlab4mic import workflows
from vlab4mic.generate import molecular_structure
from vlab4mic.utils.transform.cif_builder import do_sym_operation
from vlab4mic.utils.transform.image_convolution import (
    frame_by_volume_convolution,
    render_truncate,
    voxelize_active_fluorophores_withrange,
)


def _loop_sym_operation(pts_set, sym_operation):
    # previous implementation: one matrix-vector product per point
    rot, translation = sym_operation
    return np.array([np.matmul(rot, p) + translation for p in pts_set])


def test_sym_operation_matches_per_point_product():
    rng = np.random.default_rng(0)
    pts = rng.normal(0, 100, (5000, 3))
    rot = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    translation = rng.normal(size=3)
    expected = _loop_sym_operation(pts, [rot, translation])
    np.testing.assert_array_equal(do_sym_operation(pts, [rot, translation]), expected)


def test_sym_operation_of_no_points_is_empty():
    result = do_sym_operation([], [np.eye(3), np.zeros(3)])
    assert result.size == 0


@pytest.mark.parametrize("projection_depth", [4, 10, 200])
def test_rendering_projects_then_filters_like_filtering_the_volume(projection_depth):
    rng = np.random.default_rng(1)
    coordinates = rng.uniform(0, 60, (300, 3))
    photons = rng.poisson(50, 300)
    ranges = [(0, 60), (0, 60), (0, 40)]
    sigma, zfocus = 2.5, 20
    frame = frame_by_volume_convolution(
        coordinates, photons, ranges, 1, None, zfocus,
        projection_depth=projection_depth, render_sigma_xy=sigma,
    )
    # previous implementation: filter the volume, then project
    volume = voxelize_active_fluorophores_withrange(coordinates, photons, ranges, 1)
    filtered = gaussian_filter(
        volume, sigma=(sigma, sigma, 0),
        truncate=render_truncate(volume.shape, sigma), mode="constant",
    )
    half = min(int(projection_depth / 2), int(volume.shape[2] / 2))
    expected = filtered[:, :, zfocus - half : zfocus + half].sum(axis=2)
    np.testing.assert_allclose(frame, expected, rtol=1e-12, atol=1e-12 * expected.max())


@pytest.mark.network
def test_cif_file_is_read_once(configuration_directory, monkeypatch):
    reads = []

    def counting_reader(path):
        reads.append(path)
        return MMCIF2Dict(path)

    monkeypatch.setattr(molecular_structure, "MMCIF2Dict", counting_reader)
    structure, _ = workflows.load_structure("1XI5", configuration_directory)
    assert reads == []
    assert dict(structure.CIFdictionary) == dict(MMCIF2Dict(structure.source_file))
    assert structure._parsed_mmcif_dict is None
