import vlab4mic.generate.labelled_instance as lp
from vlab4mic import workflows
from vlab4mic.utils import data_format
import pytest
import copy
import numpy as np


def test_empty_particle():
    particle = lp.LabeledInstance()


def test_create_particle():
    particle = lp.create_particle()


# test with source and label params from fixture

structure_list = [
    "7R5K",
]
labels_per_structure = {
    "7R5K": [
        "NPC_Nup96_Cterminal_direct",
    ]
}
generic_labels = [
    "NHS_ester",
]


def test_add_probe_with_peptide_motif(experiment_7r5k_base):
    test_experiment = copy.copy(experiment_7r5k_base)
    test_experiment.remove_probes()
    list_of_proteins = test_experiment.structure.list_protein_names()
    chain_name = list_of_proteins[np.random.randint(0, len(list_of_proteins))]
    probe_name = "Antibody"
    test_experiment.add_probe(
        probe_template=probe_name,
        peptide_motif={
            "chain_name": chain_name,
            "position": "cterminal",
        },
    )
    assert len(test_experiment.probe_parameters) == 1
    assert (
        test_experiment.probe_parameters[probe_name]["label_name"]
        == probe_name
    )
    assert (
        test_experiment.probe_parameters[probe_name]["target"]["type"]
        == "Sequence"
    )
    assert (
        type(test_experiment.probe_parameters[probe_name]["target"]["value"])
        == str
    )
