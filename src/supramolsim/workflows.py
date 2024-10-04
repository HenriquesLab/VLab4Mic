from .generate.molecular_structure import build_structure_cif
from .utils.io.yaml_functions import load_yaml
from .generate.molecular_structure import MolecularReplicates
from .generate import labelled_instance as labinstance
from .generate.labels import construct_label
from .generate.coordinates_field import create_min_field
from .generate.imaging import Imager
from .utils.data_format.configuration_format import compile_modality_parameters
import os


def load_structure(structure_id: str = None, config_dir=None):
    """
    Initialise a BioPython object for the PDB/CIF ID
    and retreive available information about specific labelling.
    Assumes configuration files for structure and labels exist
    in the configuration directory provided.
    Downloads the file if necesary.

    Args:
        structure_id:  (string) 4 letter ID of structure
        config_dir: (string) absolute path for configuration files
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
        structure = build_structure_cif(
            cif_file=cif_file,
            struct_title=structure_params["title"],
            cif_id=structure_id,
        )
        print("Structure Loaded!")
        if "labels" not in structure_params.keys():
            print("Structure has no predefined labels")
            structure_params["labels"] = list(["<None>"])
        return structure, structure_params
    else:
        print("No configuration directory exists")


def particle_from_structure(
    structure: MolecularReplicates, labels=list, config_dir=None
):
    """
    Create a label object from each pair of label and fluorophore configuration files
    and this label it to create targets on the structure object.

    Args:
        structure: (MolecularReplicates) Object that represent the parsed CIF
        labels: (list) list of dictionaries. Each dictionary
        contains its label ID and Fluorophore ID
    Returns:
        particle: (LabelledInstance)
    """
    if config_dir is not None:
        label_params_list = []
        label_config_dir = os.path.join(config_dir, "labels")
        for label in labels:
            label_name = label["label_id"] + ".yaml"
            label_config_path = os.path.join(label_config_dir, label_name)
            label_object, label_params = construct_label(
                label_config_path,
                label["fluorophore_id"],
                lab_eff=label["labelling_efficiency"],
            )
            label_params_list.append(label_params)
            # print(f"Label type is: {label_params["label_type"]}")
            structure.add_label(label_object)
            # print(label_params)
            if label_params["label_type"] == "BindingLabel":
                print("Label is indirect label")
                structure.assign_normals2targets()  # default is with scaling
        inst_builder = structure.create_instance_builder()
        particle = labinstance.create_particle(
            source_builder=inst_builder, label_params_list=label_params_list
        )
        return particle


def field_from_particle(
    particle: labinstance.LabeledInstance, field_config: str = None
):
    coordinates_field = create_min_field()
    if field_config is not None:
        print("Creating field from parameter files")
        coordinates_field.init_from_file(field_config)
    coordinates_field.create_molecules_from_InstanceObject(particle)
    coordinates_field.construct_static_field()
    return coordinates_field.export_field(), coordinates_field


def create_imaging_system(
    exported_field=None, modalities_id_list: list = None, config_dir=None, **kwargs
):
    """
    Reads the configuration files for the specified imaging modalities and
    creates the imager object

    args:
        particle: LabeledInstance
        modalities_list: (list) List of modalities ID
    returns:
        imager (Optics Object) Object containing modalities and structural model
    """

    if config_dir is not None:
        if exported_field is None:
            # minimal field
            coordinates_field = create_min_field()
            exported_field = coordinates_field.export_field()
            image_generator = Imager()
            image_generator.import_field(**exported_field)
        else:
            image_generator = Imager()
            image_generator.import_field(**exported_field)
        fluo_emission = dict()
        for fluo in exported_field["field_emitters"].keys():
            fluo_dir = os.path.join(config_dir, "fluorophores", fluo)
            fluopath = fluo_dir + ".yaml"
            image_generator.set_fluorophores_from_file(fluopath)
            fluoprams = load_yaml(fluopath)
            fluo_emission[fluo] = fluoprams["emission"]
        for mod in modalities_id_list:
            modality = compile_modality_parameters(mod, config_dir, fluo_emission)
            image_generator.set_imaging_modality(**modality)
        return image_generator
