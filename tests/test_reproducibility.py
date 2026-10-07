from vlab4mic import experiments
import pytest
import numpy as np
import math


def test_exact_structure_labelling():
    random_seed = 25

    sample1, experiment1 = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        probe_DoL=5,
        number_of_particles=1,
        clear_experiment=True,
        random_seed=random_seed
    )
    sample2, experiment2 = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        probe_DoL=5,
        number_of_particles=1,
        clear_experiment=True,
        random_seed=random_seed
    )
    sample3, experiment3 = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        probe_DoL=5,
        number_of_particles=1,
        clear_experiment=True,
        random_seed=1
    )
    assert (experiment1.particle.emitters["Antibody"] == experiment2.particle.emitters["Antibody"]).all()
    assert (sample1["field_emitters"]["AF647"] == sample2["field_emitters"]["AF647"]).all()
    assert experiment1.particle.emitters["Antibody"].shape != experiment3.particle.emitters["Antibody"].shape


def test_virtual_sample_random_rotations():
    random_seed = 25

    sample1, experiment1 = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        probe_DoL=5,
        number_of_particles=5,
        clear_experiment=True,
        random_rotations=True,
        random_seed=random_seed
    )
    sample2, experiment2 = experiments.generate_virtual_sample(
        structure="1XI5",
        probe_template="Antibody",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
        probe_DoL=5,
        number_of_particles=5,
        random_rotations=True,
        clear_experiment=True,
        random_seed=random_seed
    )
    