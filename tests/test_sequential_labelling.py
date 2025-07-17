import pytest
from supramolsim import workflows, experiments
from supramolsim.utils import data_format


def test_primary_secondary_fromExperiment():
    sample, test_experiment = experiments.generate_virtual_sample(clear_probes=True)
    test_experiment.remove_probes()
    test_experiment.add_probe(
        probe_template = "Linker",
        probe_target_type = "Sequence",
        probe_target_value = "EQATETQ",
        )
    test_experiment.add_probe(
        probe_template = "Antibody",
        probe_target_type = "Primary",
        probe_target_value = "Linker")
    test_experiment.build(modules=["particle"])
    assert test_experiment.particle is not None
    assert test_experiment.exported_coordinate_field['field_emitters'] != {}