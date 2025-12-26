import pytest
from vlab4mic import workflows, experiments
from vlab4mic.utils import data_format


def test_primary_secondary_fromExperiment():
    sample, test_experiment = experiments.generate_virtual_sample(
        clear_probes=True
    )
    test_experiment.remove_probes()
    test_experiment.add_probe(
        probe_template="Linker",
        probe_target_type="Sequence",
        probe_target_value="EQATETQ",
    )
    test_experiment.add_probe(
        probe_template="Antibody",
        probe_target_type="Primary",
        probe_target_value="Linker",
    )
    test_experiment.build(modules=["particle"])
    assert test_experiment.particle is not None
    assert test_experiment.exported_coordinate_field["field_emitters"] != {}


def test_primaryAB_secondaryAB_fromExperiment():
    modalities = ["STED", "SMLM",]
    primary = dict(
        probe_template = "Antibody",
        probe_name="Primary-clathrin",
        probe_target_type = "Sequence",
        probe_target_value = "EQATETQ",
    )
    secondary = dict(
        probe_template = "Antibody",
        probe_name="Secondary-clathrin",
        probe_target_type = "Primary",
        probe_target_value = "Primary-clathrin"
    )
    image_outputs4, image_outputs_noiseless4, experiment = experiments.image_vsample(
        structure = "1XI5",
        primary_probe = primary,
        secondary_probe = secondary,
        multimodal=modalities,
        clear_experiment=True,
        run_simulation=True
    )

    assert len(list(experiment.probe_parameters.keys())) == 2
    assert len(list(experiment.particle.secondary.keys())) == 1 
    assert experiment.particle.sequential_labelling == True

