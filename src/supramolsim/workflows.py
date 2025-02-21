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
    This function require your configuration directory
    to contain a structure configuration file (yaml).
    It automatically downloads the CIF file of the structure if
    it does not exist locally.

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
        title = ""
        if "title" in structure_params["model"]:
            title = structure_params["model"]["title"]
        structure = build_structure_cif(
            cif_file=cif_file,
            struct_title=title,
            cif_id=structure_id,
        )
        print("Structure Loaded!")
        if "labels" not in structure_params.keys():
            print("Structure has no predefined labels")
            structure_params["labels"] = list(["<None>"])
        return structure, structure_params
    else:
        print("No configuration directory exists")


def probe_model(
    model, binding, conjugation_sites, config_dir, probe_name="probe", **kwargs
):
    anchor_point = None
    probe_epitope = dict(coordinates=None, normals=None)
    structural_model, structure_model_prams = load_structure(model["ID"], config_dir)
    # generate conjugation sites
    # print(conjugation_sites["target"]["value"])
    structural_model.gen_Targets(
        target_name=probe_name,
        target_type=conjugation_sites["target"]["type"],
        target_value=conjugation_sites["target"]["value"],
        **kwargs,
    )
    target_sites = structural_model.label_targets[probe_name]["coordinates"]
    # calculate axis for antibody by getting center of mass and paratope sequence
    if binding["paratope"]:
        structural_model.gen_Targets(
            target_name="paratope",
            target_type="Sequence",
            target_value=binding["paratope"],
            **kwargs,
        )
        paratopes = structural_model.label_targets["paratope"]["coordinates"]
        if paratopes.shape[0] > 1:
            print("more than one paratope sites recovered. Calculating average")
            paratope_site = np.mean(
                structural_model.label_targets["paratope"]["coordinates"], axis=0
            )
        else:
            paratope_site = structural_model.label_targets["paratope"]["coordinates"]
    direction_point = structural_model.assembly_refpt
    # extend the anchor from paratope
    diff = np.array([0, 0, 0])
    if binding["distance"]["to_target"] is not None:
        difference = paratope_site - direction_point
        unit_vector = difference / np.linalg.norm(difference)
        diff = unit_vector * binding["distance"]["to_target"]
    else:
        binding["distance"]["to_target"] = 0
    print(f"dif:{diff}")
    anchor_point = paratope_site + diff
    # generate a site for binding
    if "epitope" in kwargs.keys():
        print("generating epitope")
        print(kwargs["epitope"])
        structural_model.gen_Targets(
            target_name="epitope",
            target_type=kwargs["epitope"]["target"]["type"],
            target_value=kwargs["epitope"]["target"]["value"],
            **kwargs,
        )
        epitope_sites = structural_model.label_targets["epitope"]["coordinates"]
        if epitope_sites.shape[0] > 1:
            print("more than one paratope sites recovered. Calculating average")
            epitope_sites = np.mean(
                structural_model.label_targets["epitope"]["coordinates"], axis=0
            )
        probe_epitope["coordinates"] = epitope_sites

    # print(structural_model.label_targets)
    print(f"epitope site: {probe_epitope}")

    return structural_model, target_sites, anchor_point, direction_point, probe_epitope


def particle_from_structure(
    structure: MolecularReplicates, labels=list, config_dir=None
):
    """
    Create a labelled particle.
    First, build each label as label objects from the configuration
    directory.
    Each label is added to the structure object; this action
    generates a set of targets corresponding to that label.
    A labelled particle is initialised from these targets and
    their associated labels.

    Args:
        structure: (MolecularReplicates) Object that represent the parsed CIF
        labels: (list) list of dictionaries. Each dictionary
        contains its label ID and Fluorophore ID
    Returns:
        particle object: (LabelledInstance)
    """
    if config_dir is not None:
        label_params_list = []
        label_config_dir = os.path.join(config_dir, "probes")
        for label in labels:
            label_name = label["label_id"] + ".yaml"
            label_config_path = os.path.join(label_config_dir, label_name)
            if "target_info" in label.keys():
                label_object, label_params = construct_label(
                    label_config_path=label_config_path,
                    fluorophore_id=label["fluorophore_id"],
                    lab_eff=label["labelling_efficiency"],
                    target_info=label["target_info"],
                )
            else:
                label_object, label_params = construct_label(
                    label_config_path=label_config_path,
                    fluorophore_id=label["fluorophore_id"],
                    lab_eff=label["labelling_efficiency"],
                )
            # get model for antibody and add to label params
            if label_object.model:
                print("Generating conjugation sites")
                (
                    probe,
                    probe_emitter_sites,
                    anchor_point,
                    direction_point,
                    probe_epitope,
                ) = probe_model(
                    model=label_object.model,
                    binding=label_object.binding,
                    conjugation_sites=label_object.conjugation,
                    epitope=label_object.epitope,
                    config_dir=config_dir,
                )
                if anchor_point.shape == (3,):
                    print("setting new axis")
                    label_object.set_axis(pivot=anchor_point, direction=direction_point)
                if (
                    probe_epitope["coordinates"] is not None
                    and label_params["target"]["type"] != "Primary"
                ):
                    print("Generating linker from epitope site")
                    label_object.set_emitters(probe_epitope["coordinates"])
                else:
                    label_object.set_emitters(probe_emitter_sites)
                label_params["coordinates"] = label_object.gen_labeling_entity()
                print(label_params["coordinates"])
            label_params_list.append(label_params)
            structure.add_label(label_object)
            if label_params["binding"]["distance"]["to_target"] is not None:
                print("Label is indirect label")
                structure.assign_normals2targets()  # default is with scaling
            else:
                print("Label is direct")
        inst_builder = structure.create_instance_builder()
        particle = labinstance.create_particle(
            source_builder=inst_builder, label_params_list=label_params_list
        )
        return particle, label_params_list


def field_from_particle(
    particle: labinstance.LabeledInstance, field_config: str = None
):
    """
    Create a particle field from input particle object.

    A minimal field is initialised, bu default it defines a single particle
    at the middle of a square area of 1x1 micrometers. If a field
    configuration file is provided, the initialised field is adjusted from
    this configuration file.
    According to the definition of the field, one or more particle copies
    are positioned within the field.
    The exact positions of the emitters within this field is exported.

    Args:
        particle: (MolecularReplicates) Object that represent the parsed CIF
        field_config: (list) list of dictionaries. Each dictionary
    Returns:
        exported_field:
        coordinates_field:
    """
    if field_config is not None:
        print("Creating field from parameter files")
        # coordinates_field.init_from_file(field_config)
        ## if not file, use as dictionary
        coordinates_field = create_min_field(**field_config)
        # coordinates_field.init_from_file(field_config)
    else:
        coordinates_field = create_min_field()
    coordinates_field.create_molecules_from_InstanceObject(particle)
    coordinates_field.construct_static_field()
    return coordinates_field.export_field(), coordinates_field


def create_imaging_system(
    exported_field=None, modalities_id_list: list = None, config_dir=None, **kwargs
):
    """
    Create an imaging system object.

    Reads the configuration files for the specified imaging modalities
    in modalities_id_list and initialises the imager object.

    args:
        exported_field: output dictionary from the export_field method
        if the coordinate field
        modalities_id_list: (list) List of modalities IDs. The IDs must match
        the name of their configuration file
        config_dir: Configuration directory
    returns:
        imager (Imaging Object): Instance of Imaging class initialised with
        the input modalities
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
            # image_generator.set_fluorophores_from_file(fluopath)
            fluo_params = load_yaml(fluopath)
            image_generator.set_fluorophores_params(
                identifier=fluo,
                photon_yield=fluo_params["emission"]["photon_yield"],
                emission=fluo_params["emission"]["type"],
                blinking_rates=fluo_params["blinking_rates"],
            )
            fluo_emission[fluo] = fluo_params["emission"]["type"]
        modality_parameters = []
        for mod in modalities_id_list:
            modality = compile_modality_parameters(mod, config_dir, fluo_emission)
            modality_parameters.append(modality)
            image_generator.set_imaging_modality(**modality)
        return image_generator, modality_parameters


# generate several modalities results
def generate_multi_imaging_modalities(
    image_generator,
    experiment_name="multi_imaging_modalities",
    savingdir=None,
    acquisition_param: dict = None,
    write=False,
    **kwargs,
):
    """
    Generate imaging ouputs for each modality from the image generator.

    args:
        image_generator: Instance of Imaging class
        experiment_name: Name to add to each output file
        savingdir: Output directory
        acquisition_param: dictionary conatining aquisition parameters for
        each modality.
        write: If True, all output images will be writen at the savingdir
    returns:
        outputs: dictionary of image outputs per modality. If write is True, the
    """
    # kwargs will contain as keys the names of the modalities,
    # which shall correspond to the modality name
    # in the yaml file. each modality will have specific parameters for aquisitions
    # there must be a default imaging parameter for all, like 10 frames for each
    image_generator.set_experiment_name(experiment_name)
    outputs = dict()
    if acquisition_param is None:
        print("No acquisition parameters defined. Using default on all modalities")
        for mod in image_generator.modalities.keys():
            # should verify that the path exist
            if savingdir is not None:
                savingdir = savingdir + os.sep
                image_generator.set_writing_directory(savingdir)
            acq_params = format_modality_acquisition_params(save=write)
            timeseries, calibration_beads = image_generator.generate_imaging(
                modality=mod, **acq_params
            )
            outputs[mod] = timeseries
    else:
        acquisition_parameters = copy.copy(acquisition_param)
        for mod, acq_params in acquisition_parameters.items():
            print(mod, acq_params)
            if acq_params is None:
                acq_params = format_modality_acquisition_params(save=write)
            if savingdir is not None:
                savingdir = savingdir + os.sep
                image_generator.set_writing_directory(savingdir)
            for chan in acq_params["channels"]:
                print(f"imaging channel: {chan}")
                timeseries, calibration_beads = image_generator.generate_imaging(
                    modality=mod, channel=chan, **acq_params
                )
            outputs[mod] = timeseries
    return outputs
