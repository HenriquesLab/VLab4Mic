from .generate.molecular_structure import build_structure_cif
from .utils.io.yaml_functions import load_yaml
import os


def load_structure(structure_id: str = None, config_dir=None):
    """
    Reads a structure configuration file for structure_id and creates its object
    according to the ID. Downloads the file if necesary

    Args:
        structure_config: path of structrue configuration file
    Returns:
        Molecularstructure object
    """
    if config_dir is not None:
        print("Loading structure.")
        structure_dir = os.path.join(config_dir, "structures")
        config_file = structure_id + ".yaml"
        structure_configuration = os.path.join(structure_dir, config_file)
        structure_params = load_yaml(structure_configuration)
        # get CIF path
        cif_name = structure_id + ".cif"
        cif_file = os.path.join(structure_dir, cif_name)
        # build structure
        structure = build_structure_cif(cif_file=cif_file,
                                        struct_title=structure_params["title"],
                                        cif_id=structure_id)
        print("Structure Loaded!")
        if "labels" not in structure_params.keys():
            print("Structure has no predefined labels")
            structure_params["labels"] = list(["<None>"])
        return structure, structure_params
    else:
        print("No configuration directory exists")
