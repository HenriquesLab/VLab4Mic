from .generate.molecular_structure import build_structure_cif
from .utils.io.yaml_functions import load_yaml
from .generate.molecular_structure import MolecularReplicates
from .generate import labelled_instance as labinstance
from .generate.labels import construct_label
from .generate.coordinates_field import create_min_field
from .generate.imaging import Imager
from .utils.data_format.configuration_format import (
    compile_modality_parameters,
    format_modality_acquisition_params,
)
from .download import verify_structure
import os
import copy
import itertools
import matplotlib.pyplot as plt
from IPython.utils import io
from tqdm import tqdm
import numpy as np
from .utils.transform.datatype import truncate
from .utils.data_format.structural_format import label_builder_format


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
        cif_file = verify_structure(structure_id, structure_dir)
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


# generate several modalities results
def generate_multi_imaging_modalities(
    image_generator,
    experiment_name="multi_imaging_modalities",
    savingdir=None,
    acquisition_param: dict = None,
    write=False,
    **kwargs,
):
    # kwargs will contain as keys the names of the modalities,
    # which shall correspond to the modality name
    # in the yaml file. each modality will have specific parameters for aquisitions
    # there must be a default imaging parameter for all, like 10 frames for each
    image_generator.set_experiment_name(experiment_name)
    if acquisition_param is None:
        print("No acquisition parameters defined. Using default on all modalities")
        for mod in image_generator.modalities.keys():
            # should verify that the path exist
            if savingdir is not None:
                image_generator.set_writing_directory(savingdir)
                acq_params = format_modality_acquisition_params(save=write)
            timeseries, calibration_beads = image_generator.generate_imaging(
                modality=mod, **acq_params
            )
    else:
        acquisition_parameters = copy.copy(acquisition_param)
        for mod, acq_params in acquisition_parameters.items():
            print(mod, acq_params)
            if savingdir is not None:
                image_generator.set_writing_directory(savingdir)
            for chan in acq_params["channels"]:
                print(f"imaging channel: {chan}")
                timeseries, calibration_beads = image_generator.generate_imaging(
                    modality=mod, channel=chan, **acq_params
                )


def param_sweep_generator(
    linspaces_dict,
    structure,
    imager,
    structure_label,
    fluorophore_id,
    configuration_path,
    defects_eps,
    exp_name,
    output_dir,
    write=False,
):
    linspace_list = [linspaces_dict["labelling_efficiency"], linspaces_dict["defects"]]
    total_lengths = [len(i) for i in linspace_list]
    total_len = np.prod(np.array(total_lengths))
    for combination in tqdm(itertools.product(*linspace_list), total=total_len):
        with io.capture_output() as captured:
            labeff = combination[0]
            defect = combination[1]
            iteration_name = (
                exp_name
                + "LEff_"
                + str(truncate(labeff, 3))
                + "_Defect_"
                + str(truncate(defect, 3))
            )
            # print(f"eff: {labeff}, defect: {defect}")
            labels_list = []
            labels_list.append(
                label_builder_format(
                    label_id=structure_label,
                    fluorophore_id=fluorophore_id,
                    labelling_efficiency=labeff,
                )
            )
            particle = particle_from_structure(
                structure, labels_list, configuration_path
            )
            particle.add_defects(
                eps1=defects_eps["eps1"],
                xmer_neigh_distance=defects_eps["eps2"],
                deg_dissasembly=defect,
            )
            if write:
                particle.show_instance(with_sources=True)
                fig_name = iteration_name + ".png"
                name_path = os.path.join(output_dir, fig_name)
                plt.savefig(name_path)
                plt.close()
            exported_field, fieldobject = field_from_particle(
                particle, field_config=None
            )
            imager.import_field(**exported_field)
            # test_imager.show_field()

            generate_multi_imaging_modalities(
                image_generator=imager,
                experiment_name=iteration_name,
                savingdir=output_dir,
                write=write,
            )
    return fieldobject, imager
