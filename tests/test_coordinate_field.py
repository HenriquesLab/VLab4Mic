import supramolsim.generate.coordinates_field as field
import pytest


def test_create_minimal_field():
    nparticles = 24
    test_field = field.create_min_field(nparticles=nparticles)
    test_field.molecules_params["nMolecules"] == nparticles
