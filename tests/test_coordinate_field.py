import supramolsim.generate.coordinates_field as field
from supramolsim import workflows, data_format


def test_create_minimal_field():
    nparticles = 24
    test_field = field.create_min_field(nparticles=nparticles)
    test_field.molecules_params["nMolecules"] == nparticles


def test_minfield_with_particles(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id = "Mock_linker"
    target_info = dict(
        type="Sequence",
        value="EQATETQ"
    )
    fluorophore_id = "AF647"
    lab_eff = 1
    tmp_label1 = data_format.structural_format.label_builder_format(
        label_id, fluorophore_id, lab_eff, target_info
    )
    particle, label_params_list = workflows.particle_from_structure(
        structure, [tmp_label1], configuration_directory
    )
    nparticles = 24
    random_placing = True
    # minimal distance should be calculated automatically from particle
    coordinates_field = field.create_min_field(nparticles=nparticles,
                                                random_placing=random_placing)
    coordinates_field.create_molecules_from_InstanceObject(particle)
    coordinates_field.construct_static_field()
    exported_field = coordinates_field.export_field()
    assert coordinates_field.molecules_params["minimal_distance"] is not None