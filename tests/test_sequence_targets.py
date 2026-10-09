import numpy as np
import pytest

from vlab4mic import workflows


@pytest.fixture(scope="module")
def structure_1xi5(configuration_directory):
    structure, _ = workflows.load_structure("1XI5", configuration_directory)
    return structure


def _expected_sites(structure, motif):
    sites = []
    for chain in structure.struct[0]:
        for peptide in structure.ppgen.build_peptides(chain):
            sequence = str(peptide.get_sequence())
            start = sequence.find(motif)
            while start != -1:
                atoms = [a.coord for r in peptide[start:start + len(motif)] for a in r]
                sites.append(np.mean(atoms, axis=0))
                start = sequence.find(motif, start + 1)
    return np.array(sites)


@pytest.mark.network
@pytest.mark.parametrize("motif", ["EQATETQ", "EHLQLQN"])
def test_sequence_targets_average_exactly_the_motif(structure_1xi5, motif):
    # regression: the average included one residue after the motif and
    # residues were mapped through the whole chain, not the fragment
    structure = structure_1xi5
    asymmetric = structure.assymetric_defined
    structure.assymetric_defined = False  # compare asymmetric-unit sites
    try:
        structure.gen_targets_by_sequence("motif_test", motif)
    finally:
        structure.assymetric_defined = asymmetric
    found = np.array(structure.label_targets["motif_test"]["coordinates"])
    expected = _expected_sites(structure, motif)
    assert found.shape == expected.shape
    np.testing.assert_allclose(
        np.sort(found, axis=0), np.sort(expected, axis=0), atol=1e-3
    )


@pytest.mark.network
def test_absent_motif_gives_no_targets(structure_1xi5):
    structure_1xi5.gen_targets_by_sequence("absent", "WWWWWWWW")
    assert len(structure_1xi5.label_targets["absent"]["coordinates"]) == 0
