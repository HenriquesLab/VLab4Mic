import supramolsim.generate.coordinates_field as field
from supramolsim import workflows, experiments
from supramolsim.utils import data_format
import numpy as np
import copy

def test_virtual_sample_params(experiment_7r5k_base):
    copy_exp = copy.deepcopy(experiment_7r5k_base)
    copy_exp.set_virtualsample_params(number_of_particles=2)
    copy_exp.build(modules=["virtualsample"])


def test_create_minimal_field():
    number_of_particles = 24
    test_field = field.create_min_field(number_of_particles=number_of_particles)
    test_field.molecules_params["nMolecules"] == number_of_particles
    test_field.change_number_of_molecules(25)
    test_field.show_field()
    test_field.expand_isotropically(factor=2)

def test_gen_positions_from_image(experiment_7r5k_base):
    copy_exp = copy.deepcopy(experiment_7r5k_base)
    img_mask = np.random.rand(24,24)
    p = 0.9
    img_mask[img_mask >= p] = 1
    img_mask[img_mask < p] = 0
    copy_exp.use_image_for_positioning(img_mask, mode="mask", pixelsize=100)



def test_nparticles_constraints():
    number_of_particles = 200
    sample, expr = experiments.generate_virtual_sample(number_of_particles=number_of_particles)
    assert expr.coordinate_field.molecules_params["minimal_distance"] is not None
    assert len(expr.coordinate_field.molecules) < number_of_particles
    assert expr.coordinate_field.molecules_params["nMolecules"] < number_of_particles