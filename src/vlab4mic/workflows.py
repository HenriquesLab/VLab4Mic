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
from scipy.spatial.distance import pdist
from .utils.transform.datatype import truncate
from .utils.data_format.structural_format import label_builder_format
from pathlib import Path

def load_structure(structure_id: str = None, config_dir=None, structure_path=None):
    """
    Initialise a BioPython object for the PDB/CIF ID and retrieve available information about specific labelling.

    This function requires your configuration directory to contain a structure configuration file (yaml).
    It automatically downloads the CIF file of the structure if it does not exist locally.

    Parameters
    ----------
    structure_id : str, optional
        4-letter ID of structure.
    config_dir : str, optional
        Absolute path for configuration files.

    Returns
    -------
    structure : Molecularstructure
        The loaded structure object.
    structure_params : dict
        Structure parameters loaded from configuration.
    """
    if config_dir is not None:
        if structure_path is not None:
            print("Loading structure from local file")
            structure = build_structure_cif(
                cif_file=structure_path, struct_title="", cif_id=structure_id
            )
            structure_params = dict(
                model = {"ID": structure_id, "format": "CIF", "title": "", "database": ""},
            )
            print("Structure Loaded!")
            return structure, structure_params
        else:
            print("Loading structure.")
            structure_dir = os.path.join(config_dir, "structures")
            config_file = structure_id + ".yaml"
            structure_configuration = os.path.join(structure_dir, config_file)
            structure_configuration_check = Path(structure_configuration)
            if structure_configuration_check.exists():
                structure_params = load_yaml(structure_configuration)
            else:
                template = os.path.join(structure_dir, "_template_.yaml")
                structure_params = load_yaml(template)
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
    """
    Generate probe model and anchor points for a given structure and binding configuration.

    Parameters
    ----------
    model : dict
        Model information dictionary.
    binding : dict
        Binding configuration dictionary.
    conjugation_sites : dict
        Conjugation site information.
    config_dir : str
        Path to configuration directory.
    probe_name : str, optional
        Name for the probe (default "probe").
    **kwargs
        Additional keyword arguments, e.g., epitope.

    Returns
    -------
    structural_model : MolecularReplicates
        The structure with generated targets.
    target_sites : numpy.ndarray
        Coordinates of target sites.
    anchor_point : numpy.ndarray
        Anchor point for the probe.
    direction_point : numpy.ndarray
        Direction point for the probe.
    probe_epitope : dict
        Dictionary with epitope coordinates and normals.
    """
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
    if binding["distance"]["to_target"] is not None and binding["paratope"]:
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
        try:
            #print(kwargs["epitope"])
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
            print(f"epitope site: {probe_epitope}")
        except:
            print("No epitope generated")
            print(kwargs["epitope"])

    # print(structural_model.label_targets)
    

    return structural_model, target_sites, np.squeeze(anchor_point, axis=None), np.squeeze(direction_point, axis=None), probe_epitope


def particle_from_structure(
    structure: MolecularReplicates, labels=list, config_dir=None
):
    """
    Create a labelled particle from a structure and label definitions.

    Parameters
    ----------
    structure : MolecularReplicates
        Object representing the parsed CIF structure.
    labels : list of dict
        List of label dictionaries, each containing label ID and fluorophore ID.
    config_dir : str, optional
        Path to configuration directory.

    Returns
    -------
    particle : LabelledInstance
        The created labelled particle.
    label_params_list : list of dict
        List of label parameter dictionaries.
    """
    if config_dir is not None:
        label_params_list = []
        label_config_dir = os.path.join(config_dir, "probes")
        for label in labels:
            #label_name = label["label_id"] + ".yaml"
            #label_config_path = os.path.join(label_config_dir, label_name)
            label_object, label_params = construct_label(
                label_config_dictionary=label
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
                    print(anchor_point, direction_point)
                    direction_vector = direction_point - anchor_point
                    label_object.set_axis(pivot=anchor_point, direction=direction_vector)
                    print(f'label object axis: {label_object.params["axis"]}')
                if (
                    probe_epitope["coordinates"] is not None
                    and label_params["as_linker"]
                ):
                    # as_linker is True when this probe has a secondary linker
                    print("Generating linker from epitope site")
                    label_object.set_emitters(probe_epitope["coordinates"])
                else:
                    print(probe_epitope["coordinates"], label_params["as_linker"])
                    print("Setting emitters from probe conjugation sites")
                    label_object.set_emitters(probe_emitter_sites)
                label_params["coordinates"] = label_object.gen_labeling_entity()
                print(label_params["coordinates"])
            
            structure.add_label(label_object)
            if label_params["binding"]["distance"]["to_target"] is not None:
                print("Label is indirect label")
                structure.assign_normals2targets()  # default is with scaling
            else:
                print("Label is direct")
            if label_params["binding"]["distance"]["between_targets"] == "estimate":
                distances = pdist(label_params["coordinates"])
                label_params["binding"]["distance"]["between_targets"] = np.max(distances)
            label_params_list.append(label_params)
        inst_builder = structure.create_instance_builder()
        particle = labinstance.create_particle(
            source_builder=inst_builder, label_params_list=label_params_list
        )
        return particle, label_params_list


def field_from_particle(
    particle: labinstance.LabeledInstance, field_config: str = None, **kwargs
):
    """
    Create a particle field from an input particle object.

    A minimal field is initialised; by default it defines a single particle at the middle of a square area of 1x1 micrometers.
    If a field configuration file is provided, the initialised field is adjusted from this configuration file.

    Parameters
    ----------
    particle : LabeledInstance
        Particle object to place in the field.
    field_config : dict or str, optional
        Field configuration dictionary or file path.
    **kwargs
        Additional keyword arguments for field creation.

    Returns
    -------
    exported_field : dict
        Exported field dictionary.
    coordinates_field : CoordinatesField
        The coordinates field object.
    """
    if field_config is not None:
        print("Creating field from parameter files")
        # coordinates_field.init_from_file(field_config)
        ## if not file, use as dictionary
        coordinates_field = create_min_field(**field_config, **kwargs)
        # coordinates_field.init_from_file(field_config)
    else:
        coordinates_field = create_min_field(**kwargs)
    coordinates_field.create_molecules_from_InstanceObject(particle)
    coordinates_field.construct_static_field()
    return coordinates_field.export_field(), coordinates_field


def create_imaging_system(
    exported_field=None, modalities_id_list: list = None, config_dir=None, mod_params=dict(), fluorophore_parameters=None, **kwargs
):
    """
    Create an imaging system object.

    Reads the configuration files for the specified imaging modalities in modalities_id_list and initialises the imager object.

    Parameters
    ----------
    exported_field : dict, optional
        Output dictionary from the export_field method of the coordinate field.
    modalities_id_list : list of str, optional
        List of modality IDs. The IDs must match the name of their configuration file.
    config_dir : str, optional
        Path to configuration directory.
    mod_params : dict, optional
        Dictionary of modality parameters.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    image_generator : Imager
        Instance of Imaging class initialised with the input modalities.
    modality_parameters : list of dict
        List of modality parameter dictionaries.
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
            #fluo_dir = os.path.join(config_dir, "fluorophores", fluo)
            #fluopath = fluo_dir + ".yaml"
            # image_generator.set_fluorophores_from_file(fluopath)
            #try:
            #    fluo_params = load_yaml(fluopath)
            #except:
            fluo_params = fluorophore_parameters[fluo]
            image_generator.set_fluorophores_params(
                identifier=fluo,
                photon_yield=fluo_params["emission"]["photon_yield"],
                emission=fluo_params["emission"]["type"],
                blinking_rates=fluo_params["blinking_rates"],
            )
            fluo_emission[fluo] = fluo_params["emission"]["type"]
        modality_parameters = []
        for mod in modalities_id_list:
            if mod in mod_params.keys():
                image_generator.set_imaging_modality(**mod_params[mod])
                modality_parameters.append(mod_params[mod])
            else:
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
    Generate imaging outputs for each modality from the image generator.

    Parameters
    ----------
    image_generator : Imager
        Instance of Imaging class.
    experiment_name : str, optional
        Name to add to each output file. Default is "multi_imaging_modalities".
    savingdir : str, optional
        Output directory for saving images.
    acquisition_param : dict, optional
        Dictionary containing acquisition parameters for each modality.
    write : bool, optional
        If True, all output images will be written to savingdir.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    outputs : dict
        Dictionary of image outputs per modality.
    """
    # kwargs will contain as keys the names of the modalities,
    # which shall correspond to the modality name
    # in the yaml file. each modality will have specific parameters for aquisitions
    # there must be a default imaging parameter for all, like 10 frames for each
    image_generator.set_experiment_name(experiment_name)
    outputs = dict()
    outputs_noiseless = dict()
    if acquisition_param is None:
        print("No acquisition parameters defined. Using default on all modalities")
        for mod in image_generator.modalities.keys():
            # should verify that the path exist
            if savingdir is not None:
                savingdir = savingdir + os.sep
                image_generator.set_writing_directory(savingdir)
            acq_params = format_modality_acquisition_params(save=write)
            timeseries, calibration_beads, timeseries_noiseless, calibration_beads_noiseless  = image_generator.generate_imaging(
                modality=mod, **acq_params
            )
            outputs[mod] = timeseries
            outputs_noiseless[mod] = timeseries_noiseless
    else:
        acquisition_parameters = copy.copy(acquisition_param)
        for mod, acq_param in acquisition_parameters.items():
            #print(mod, acq_params)
            if acq_param is None:
                acq_params = format_modality_acquisition_params()
            else:
                acq_params = format_modality_acquisition_params(**acq_param)
            acq_params["save"] = write
            if savingdir is not None:
                savingdir = savingdir + os.sep
                image_generator.set_writing_directory(savingdir)
            for chan in acq_params["channels"]:
                print(f"imaging channel: {chan}")
                timeseries, calibration_beads, timeseries_noiseless, calibration_beads_noiseless = image_generator.generate_imaging(
                    modality=mod, channel=chan, **acq_params
                )
            outputs[mod] = timeseries
            outputs_noiseless[mod] = timeseries_noiseless
    return outputs, outputs_noiseless
