import pytest
from supramolsim import workflows, data_format


probes_1XI5_primaries = []
target_info1 = dict(type="Sequence", value="EQATETQ")
fluorophore_id1 = "AF647"
lab_eff1 = 1
linker1 = data_format.structural_format.label_builder_format(
    "Mock_linker", fluorophore_id1, lab_eff1, target_info1
)
antibody = data_format.structural_format.label_builder_format(
    "Mock_antibody", fluorophore_id1, lab_eff1, target_info1, as_linker=False
)
antibody_as_linker = data_format.structural_format.label_builder_format(
    "Mock_antibody", fluorophore_id1, lab_eff1, target_info1, as_linker=True
)
probes_1XI5_primaries = [linker1, antibody, antibody_as_linker]

# parameters for sequential labelling
probes_1XI5_primaries = [linker1, antibody_as_linker]
secondaries_ids = ["Mock_linker", "Mock_antibody"]


@pytest.mark.parametrize("probe_primary", probes_1XI5_primaries)
def test_primaries_only(configuration_directory, probe_primary):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    labels_list = [
        probe_primary,
    ]
    particle, label_params_list = workflows.particle_from_structure(
        structure, labels_list, configuration_directory
    )
    assert particle.get_ref_point().shape == (3,)


@pytest.mark.parametrize("primary", probes_1XI5_primaries)
@pytest.mark.parametrize("secondary_id", secondaries_ids)
def test_primary_with_secondary(configuration_directory, primary, secondary_id):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    target_info_secondary = dict(type="Primary", value=primary["label_id"])
    fluorophore_id1 = "AF488"
    lab_eff1 = 1
    secondary_probe = data_format.structural_format.label_builder_format(
        secondary_id, fluorophore_id1, lab_eff1, target_info_secondary
    )
    labels_list = [primary, secondary_probe]
    particle, label_params_list = workflows.particle_from_structure(
        structure, labels_list, configuration_directory
    )
    assert len(particle.secondary.keys()) > 0


def test_linker_secondary(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id1 = "Mock_linker"
    target_info1 = dict(type="Sequence", value="EQATETQ")
    fluorophore_id1 = "AF647"
    lab_eff1 = 1
    tmp_label1 = data_format.structural_format.label_builder_format(
        label_id1, fluorophore_id1, lab_eff1, target_info1
    )

    label_id2 = "Mock_antibody"
    target_info2 = dict(type="Primary", value="Mock_linker")
    fluorophore_id2 = "AF488"
    lab_eff2 = 1
    tmp_label2 = data_format.structural_format.label_builder_format(
        label_id2, fluorophore_id2, lab_eff2, target_info2, as_linker=False
    )
    particle2, label_params_list = workflows.particle_from_structure(
        structure, [tmp_label1, tmp_label2], configuration_directory
    )
    assert len(particle2.secondary.keys()) > 0
    epitopes = particle2.source["targets"][label_id1]["coordinates"].shape
    emitters = particle2.emitters[label_id1].shape
    assert epitopes != emitters


def test_primary_secondary(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id1 = "Mock_antibody"
    target_info1 = dict(type="Sequence", value="EQATETQ")
    fluorophore_id1 = "AF647"
    lab_eff1 = 1
    tmp_label1 = data_format.structural_format.label_builder_format(
        label_id1, fluorophore_id1, lab_eff1, target_info1, as_linker=True
    )

    label_id2 = "Mock_antibody"
    target_info2 = dict(type="Primary", value="Mock_antibody")
    fluorophore_id2 = "AF488"
    lab_eff2 = 1
    tmp_label2 = data_format.structural_format.label_builder_format(
        label_id2, fluorophore_id2, lab_eff2, target_info2
    )
    particle2, label_params_list = workflows.particle_from_structure(
        structure, [tmp_label1, tmp_label2], configuration_directory
    )
    assert len(particle2.secondary.keys()) > 0
    epitopes = particle2.source["targets"][label_id1]["coordinates"].shape
    emitters = particle2.emitters[label_id1].shape
    assert epitopes != emitters
