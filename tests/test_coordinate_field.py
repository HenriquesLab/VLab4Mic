import supramolsim.generate.coordinates_field as field
from supramolsim import workflows
from supramolsim.utils import data_format
import numpy as np
import copy


def test_create_minimal_field():
    nparticles = 24
    test_field = field.create_min_field(nparticles=nparticles)
    test_field.molecules_params["nMolecules"] == nparticles
    test_field.change_number_of_molecules(25)
    test_field.show_field()
    test_field.expand_isotropically(factor=2)

def test_gen_positions_from_image():
    img_mask = np.random.rand(24,24)
    p = 0.9
    img_mask[img_mask >= p] = 1
    img_mask[img_mask < p] = 0
    xyz_relative, image_physical_size = field.gen_positions_from_image(img_mask, pixelsize = 100)


def test_minfield_with_particles(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id = "Linker"
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



def test_nparticles_constraints(configuration_directory):
    structure_id = "1XI5"
    structure, structure_param = workflows.load_structure(
        structure_id, configuration_directory
    )
    label_id = "Linker"
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
    nparticles = 1000
    random_placing = True
    # minimal distance should be calculated automatically from particle
    coordinates_field = field.create_min_field(nparticles=nparticles,
                                                random_placing=random_placing)
    coordinates_field.create_molecules_from_InstanceObject(particle)
    coordinates_field.construct_static_field()
    exported_field = coordinates_field.export_field()
    assert coordinates_field.molecules_params["minimal_distance"] is not None
    assert len(coordinates_field.molecules) < nparticles
    #
    minimal_distance=1000
    coordinates_field2 = field.create_min_field(nparticles=nparticles,
                                                random_placing=random_placing,
                                                minimal_distance=minimal_distance)
    coordinates_field2.create_molecules_from_InstanceObject(particle)
    coordinates_field2.construct_static_field()
    exported_field2 = coordinates_field2.export_field()
    assert coordinates_field2.molecules_params["minimal_distance"] is not None
    assert len(coordinates_field2.molecules) < nparticles