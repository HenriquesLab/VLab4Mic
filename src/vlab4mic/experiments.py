from dataclasses import dataclass, field, fields
from typing import List, Dict
import numpy as np
import vlab4mic
from .workflows import (
    load_structure,
    create_imaging_system,
    particle_from_structure,
    field_from_particle,
    generate_multi_imaging_modalities,
)
from .generate import coordinates_field
from .utils.data_format.structural_format import label_builder_format
from .utils.data_format import configuration_format
import os
from vlab4mic.utils.io.yaml_functions import load_yaml
import numpy as np
import os
import copy
from vlab4mic.utils.io import yaml_functions
from IPython.utils import io
from pathlib import Path
import sys

IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    output_path = "/content/vlab4mic_outputs"
else:
    output_path = Path.home() / "vlab4mic_outputs"


if not os.path.exists(output_path):
    os.makedirs(output_path)


@dataclass
class ExperimentParametrisation:
    experiment_id: str = "vLab4mic_experiment"
    structure_id: str = "1XI5"
    structure_format: str = "CIF"
    configuration_path: str = ""
    structure_label: str = "NHS_ester"
    fluorophore_id: str = "AF647"
    coordinate_field_id: str = "square1x1um_randomised"
    selected_mods: Dict[str, int] = field(default_factory=dict)
    imaging_modalities: Dict[str, int] = field(default_factory=dict)
    probe_parameters: Dict[str, int] = field(default_factory=dict)
    incomplete_labelling_eps: Dict[str, int] = field(default_factory=dict)
    sweep_pars: Dict[str, int] = field(default_factory=dict)
    objects_created: Dict[str, int] = field(default_factory=dict)
    output_directory: str = None
    example_structures = ["3J3Y", "7R5K", "1XI5", "8GMO"]
    example_modalities = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

    def __post_init__(self):
        self.output_directory = str(output_path)
        pck_dir = os.path.dirname(os.path.abspath(vlab4mic.__file__))
        local_dir = os.path.join(pck_dir, "configs")
        self.configuration_path = local_dir
        self.structure_path = None
        # keep track of objects created
        self.objects_created = dict(
            structure=False,
            particle=False,
            exported_coordinate_field=False,
            coordinate_field=False,
            imager=False,
            output_reference=False,
        )
        self.virtualsample_params = dict()
        self.incomplete_labelling_eps["incomplete_labelling"] = None
        self.incomplete_labelling_eps["eps1"] = None
        self.incomplete_labelling_eps["eps2"] = None
        self.incomplete_labelling_eps["use_incomplete_labelling"] = False
        # read information of local modalities configuration
        modalities_dir = os.path.join(local_dir, "modalities")
        modalities_names_list = []
        modality_parameters = {}
        self.config_modalities = dict()
        for mods in os.listdir(modalities_dir):
            if (
                os.path.splitext(mods)[-1] == ".yaml"
                and "_template" not in mods
            ):
                modalities_names_list.append(os.path.splitext(mods)[0])
        for mod in modalities_names_list:
            mod_info, mod_configuration = (
                configuration_format.compile_modality_parameters(
                    mod, local_dir
                )
            )
            modality_parameters[mod] = mod_info
            self.config_modalities[mod] = mod_configuration
        self.local_modalities_names = modalities_names_list
        self.local_modalities_parameters = modality_parameters
        probes_dir = os.path.join(local_dir, "probes")
        structure_dir = os.path.join(local_dir, "structures")
        self.config_probe_params = {}
        self.config_probe_models_names = []
        self.config_global_probes_names = []
        self.config_probe_per_structure_names = {}
        for p_file in os.listdir(probes_dir):
            if (
                os.path.splitext(p_file)[-1] == ".yaml"
                and "_template" not in p_file
            ):
                label_config_path = os.path.join(probes_dir, p_file)
                label_parmeters = vlab4mic.load_yaml(label_config_path)
                # print(label_parmeters)
                lablname = os.path.splitext(p_file)[0]
                if "Mock" in label_parmeters["known_targets"]:
                    self.config_probe_models_names.append(lablname)
                    self.config_probe_params[lablname] = label_parmeters
                elif "Generic" in label_parmeters["known_targets"]:
                    self.config_global_probes_names.append(lablname)
                    self.config_probe_params[lablname] = label_parmeters
                else:
                    self.config_probe_params[lablname] = label_parmeters
                    for struct in label_parmeters["known_targets"]:
                        if (
                            struct
                            in self.config_probe_per_structure_names.keys()
                        ):
                            self.config_probe_per_structure_names[
                                struct
                            ].append(lablname)
                        else:
                            self.config_probe_per_structure_names[struct] = [
                                lablname,
                            ]
        self.demo_structures = []
        # get available structure IDs
        self.structures_info_list = dict()
        structure_dir = os.path.join(self.configuration_path, "structures")
        fluorophores_dir = os.path.join(
            self.configuration_path, "fluorophores"
        )
        self.fluorophore_parameters = dict()
        for file in os.listdir(fluorophores_dir):
            if (
                os.path.splitext(file)[-1] == ".yaml"
                and "_template" not in file
            ):
                fluorophore_id = os.path.splitext(file)[0]
                self.fluorophore_parameters[fluorophore_id] = load_yaml(
                    os.path.join(fluorophores_dir, file)
                )
        self.fluorophore_parameters_template = load_yaml(
            os.path.join(fluorophores_dir, "_template_.yaml")
        )
        probes_dir = os.path.join(self.configuration_path, "probes")
        modalities_dir = os.path.join(self.configuration_path, "modalities")
        for file in os.listdir(structure_dir):
            if (
                os.path.splitext(file)[-1] == ".yaml"
                and "_template" not in file
            ):
                structure_params = load_yaml(os.path.join(structure_dir, file))
                struct_id = structure_params["model"]["ID"]
                if struct_id in self.example_structures:
                    strict_title = structure_params["model"]["title"]
                    id_title = struct_id + ": " + strict_title
                    self.structures_info_list[id_title] = struct_id
                    self.demo_structures.append(id_title)
        self.config_directories = dict(
            structure=structure_dir,
            fluorophores=fluorophores_dir,
            probes=probes_dir,
            modalities=modalities_dir,
            base=self.configuration_path,
        )
        param_settings_file = os.path.join(
            self.config_directories["base"], "parameter_settings.yaml"
        )
        self.param_settings = yaml_functions.load_yaml(param_settings_file)
        self.results = dict()
        self.create_example_experiment()
        self.modality_noise_images = dict()

    def select_structure(self, structure_id="1XI5", build=True, structure_path:str = None):
        """
        Select a molecular structure by its identifier and optionally build the structure module.

        Parameters
        ----------
        :param structure_id : str, optional
            The identifier of the structure to select. Default is "1XI5".
        :param build : bool, optional
            If True, triggers the build process for the structure module. Default is True.

        Returns
        -------
        None
        """
        if structure_path is not None:
            self.structure_path = structure_path
            self.structure_id = structure_id
            self.structure_format = structure_path.split(".")[-1].upper()
        else:
            self.structure_id = structure_id
        if build:
            self.build(modules=["structure"])

    def clear_structure(self):
        """
        Clear the current structure by resetting the structure ID and related parameters.
        This method sets the structure ID to None and resets the structure object if it exists.
        Returns
        -------
        None
        """
        self.structure_id = None
        self.structure_path = None
        if self.generators_status("structure"):
            self.structure = None
        self.objects_created["structure"] = False
        print("Structure cleared")

    def add_modality(self, modality_name, save=False, modality_template=None, **kwargs):
        """
        Add a new imaging modality to the experiment.

        If the specified modality name exists in the local modalities, it copies the parameters from the local template.
        Otherwise, it uses the 'Widefield' modality parameters as a template to create a new modality.
        Additional parameters for the modality can be specified via keyword arguments and will override the defaults.
        Optionally, the modality output can be saved by setting `save` to True.

        Parameters
        ----------
        :param modality_name : str
            The name of the modality to add.
        :param save : bool, optional
            If True, saves the modality output. Default is False.
        :param **kwargs
            Arbitrary keyword arguments to override or set specific modality parameters.

        Returns
        -------
        None
        """
        if modality_name in self.local_modalities_names:
            self.imaging_modalities[modality_name] = copy.deepcopy(
                self.local_modalities_parameters[modality_name]
            )
            for param, value in kwargs.items():
                self.imaging_modalities[modality_name][param] = value
            self.set_modality_acq(modality_name, save=save)
        else:
            print(
                f"Modality {modality_name} not found in demo modalities. Using template to create new."
            )
            if modality_template is not None:
                template_modality = modality_template
            else:
                template_modality = "Widefield"
            print(f"Using {template_modality} as template.")
            self.imaging_modalities[modality_name] = copy.deepcopy(
                self.local_modalities_parameters[template_modality]
            )
            self.imaging_modalities[modality_name]["modality"] = modality_name
            for param, value in kwargs.items():
                self.imaging_modalities[modality_name][param] = value
            self.set_modality_acq(modality_name, save=save)
            self.local_modalities_names.append(modality_name)
            self.local_modalities_parameters[modality_name] = copy.deepcopy(
                self.imaging_modalities[modality_name]
            )

    def update_modality(
        self,
        modality_name,
        pixelsize_nm: int = None,
        lateral_resolution_nm: int = None,
        axial_resolution_nm: int = None,
        psf_voxel_nm: int = None,
        depth_of_field_nm: int = None,
        remove=False,
        **kwargs,
    ):
        """
        Update or remove an imaging modality's parameters.

        This method allows updating specific parameters of an imaging modality, such as pixel size, lateral and axial resolution, and PSF voxel size.
        If `remove` is True, the modality is removed from the internal dictionaries.

        Parameters
        ----------
        :param modality_name : str
            The name of the imaging modality to update or remove.
        :param pixelsize_nm : int, optional
            The new pixel size in nanometers. If provided, updates the detector pixel size.
        :param lateral_resolution_nm : int, optional
            The new lateral resolution in nanometers. If provided, updates the lateral standard deviations of the PSF.
        :param axial_resolution_nm : int, optional
            The new axial resolution in nanometers. If provided, updates the axial standard deviation of the PSF.
        :param psf_voxel_nm : int, optional
            The new PSF voxel size in nanometers. If provided, updates the PSF voxel size for all axes.
        :param remove : bool, optional
            If True, removes the specified modality from the internal dictionaries. Default is False.

        Notes
        -----
        After updating any parameters, the imaging modality is reconfigured in the imager.
        """
        if remove:
            self.imaging_modalities.pop(modality_name, None)
            self.selected_mods.pop(modality_name, None)
        else:
            changes = False
            if pixelsize_nm is not None:
                self.imaging_modalities[modality_name]["detector"][
                    "pixelsize"
                ] = (pixelsize_nm / 1000)
                changes = True
            if lateral_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name][
                    "psf_params"
                ]["voxelsize"][0]
                self.imaging_modalities[modality_name]["psf_params"][
                    "std_devs"
                ][0] = (lateral_resolution_nm / voxel_size)
                self.imaging_modalities[modality_name]["psf_params"][
                    "std_devs"
                ][1] = (lateral_resolution_nm / voxel_size)
                changes = True
            if axial_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name][
                    "psf_params"
                ]["voxelsize"][0]
                self.imaging_modalities[modality_name]["psf_params"][
                    "std_devs"
                ][2] = (axial_resolution_nm / voxel_size)
                changes = True
            if psf_voxel_nm is not None:
                self.imaging_modalities[modality_name]["psf_params"][
                    "voxelsize"
                ] = [
                    psf_voxel_nm,
                    psf_voxel_nm,
                    psf_voxel_nm,
                ]
                changes = True
            if depth_of_field_nm is not None:
                depth_in_slices = None
                voxel_size = self.local_modalities_parameters[modality_name][
                    "psf_params"
                ]["voxelsize"][2]
                depth = int(depth_of_field_nm / voxel_size)
                self.imaging_modalities[modality_name]["psf_params"][
                    "depth"
                ] = depth
                changes = True
            if changes:
                self.imager.set_imaging_modality(
                    **self.imaging_modalities[modality_name]
                )

    def set_modality_acq(
        self,
        modality_name,
        exp_time=0.001,
        noise=True,
        save=False,
        nframes=1,
        channels=[
            "ch0",
        ],
        **kwargs,
    ):
        """
        Configure and set acquisition parameters for a specified imaging modality.

        Parameters
        ----------
        :param modality_name : str
            The name of the imaging modality to configure.
        :param exp_time : float, optional
            Exposure time for the acquisition in seconds. Default is 0.001.
        :param noise : bool, optional
            Whether to include noise in the acquisition. Default is True.
        :param save : bool, optional
            Whether to save the acquired data. Default is False.
        :param nframes : int, optional
            Number of frames to acquire. Default is 1.
        :param channels : list of str, optional
            List of channel names to use. Default is ["ch0"].
        :param **kwargs
            Additional keyword arguments for modality-specific parameters.

        Notes
        -----
        If the specified modality is not available in `self.imaging_modalities`, a message is printed and no action is taken.
        Acquisition parameters are formatted using `configuration_format.format_modality_acquisition_params` and stored in `self.selected_mods`.
        """
        if modality_name in self.imaging_modalities.keys():
            self.selected_mods[modality_name] = (
                configuration_format.format_modality_acquisition_params(
                    exp_time, noise, save, nframes, channels
                )
            )
        else:
            print("Modality not selected")

    def clear_modalities(self):
        """
        Clear all selected imaging modalities and their acquisition parameters.
        This method resets the `imaging_modalities` and `selected_mods` dictionaries to empty states,
        effectively removing all configured modalities from the experiment.
        Returns
        -------
        None
        """
        modality_names = list(self.imaging_modalities.keys())
        for modality_name in modality_names:
            self.update_modality(modality_name, remove=True)
        if self.generators_status("imager"):
            self.imager = None
            self.objects_created["imager"] = False
        print("All modalities cleared")

    def reset_to_defaults(self, module="acquisitions", **kwargs):
        """
        Reset acquisition parameters or other module settings to their default values.

        Parameters
        ----------
        :param module : str, optional
            The module to reset. Default is "acquisitions".
        :param **kwargs
            Additional keyword arguments passed to the reset function.

        Returns
        -------
        None
        """
        if module == "acquisitions":
            for mod_name in self.selected_mods.keys():
                self.set_modality_acq(mod_name, **kwargs)

    def _build_structure(self, keep=True):
        """
        Build the structure object for the experiment.

        Parameters
        ----------
        :param keep : bool, optional
            If True, store the structure in the experiment. Default is True.

        Returns
        -------
        None
        """
        if self.structure_id and self.structure_path:
            struct, struct_param = load_structure(
                self.structure_id, self.configuration_path, self.structure_path, self.structure_format
            )
            self.structure = struct
            self.objects_created["structure"] = True
        elif self.structure_id:
            struct, struct_param = load_structure(
                self.structure_id, self.configuration_path
            )
            self.structure = struct
            self.objects_created["structure"] = True
        else:
            print("No structure ID")

    def _build_label(self, lab_eff=1, keep=False, **kwargs):
        """
        Create a list with one dictionary with the setup info for creating a label.

        Assumes there is only one label stored in the attribute "structure_label" and one in "fluorophore_id".

        Parameters
        ----------
        :param lab_eff : float, optional
            Labelling efficiency. Default is 1.
        :param keep : bool, optional
            If True, keep the label. Default is False.
        :param **kwargs
            Additional keyword arguments for label creation.

        Returns
        -------
        list
            List of label dictionaries.
        """
        labels_list = []
        for probe_name, probe_params in self.probe_parameters.items():
            labels_list.append(
                copy.deepcopy(self.probe_parameters[probe_name])
                #
            )
        return labels_list

    def _check_if_incomplete_labellings(self):
        """
        Check if incomplete_labelling parameters are set and update the use_incomplete_labelling flag.

        Returns
        -------
        None
        """
        if (
            self.incomplete_labelling_eps["incomplete_labelling"]
            and self.incomplete_labelling_eps["eps1"]
            and self.incomplete_labelling_eps["eps2"]
        ):
            self.incomplete_labelling_eps["use_incomplete_labelling"] = True
        else:
            self.incomplete_labelling_eps["use_incomplete_labelling"] = False

    def _build_particle(self, lab_eff=1.0, incomplete_labelling_build=None, keep=False):
        """
        Build the particle object for the experiment, optionally adding incomplete_labellings.

        Parameters
        ----------
        :param lab_eff : float, optional
            Labelling efficiency. Default is 1.0.
        :param incomplete_labelling_build : float or None, optional
            Defect parameter to use. Default is None.
        :param keep : bool, optional
            If True, store the particle in the experiment. Default is False.

        Returns
        -------
        particle or None
            The particle object if created, else None.
        """
        if self.generators_status("structure"):
            self._check_if_incomplete_labellings()
            labels_list = self._build_label(lab_eff=lab_eff)
            if len(labels_list) > 0:
                particle, label_params_list = particle_from_structure(
                    self.structure, labels_list, self.configuration_path
                )
                if self.incomplete_labelling_eps["use_incomplete_labelling"]:
                    if incomplete_labelling_build is not None:
                        incomplete_labelling = incomplete_labelling_build
                    else:
                        incomplete_labelling = self.incomplete_labelling_eps["incomplete_labelling"]
                    particle.add_incomplete_labellings(
                        eps1=self.incomplete_labelling_eps["eps1"],
                        xmer_neigh_distance=self.incomplete_labelling_eps["eps2"],
                        deg_dissasembly=incomplete_labelling,
                    )
                else:
                    particle.add_incomplete_labellings(
                        deg_dissasembly=0,
                    )
                if keep:
                    self.particle = particle
                    self.objects_created["particle"] = True
                return particle

    def _build_coordinate_field(
        self,
        use_self_particle=True,
        keep=False,
        coordinate_field_path=None,
        **kwargs,
    ):
        """
        Build the coordinate field for the experiment, either from the current particle or as a minimal field.

        Parameters
        ----------
        :param use_self_particle : bool, optional
            If True, use the current particle to build the field. Default is True.
        :param keep : bool, optional
            If True, store the field in the experiment. Default is False.
        :param coordinate_field_path : str or None, optional
            Path to a coordinate field file. Default is None.
        :param **kwargs
            Additional keyword arguments for field creation.

        Returns
        -------
        exported_field or None
            The exported field if created, else None.
        """
        self.exported_coordinate_field = None
        if use_self_particle and self.generators_status("particle"):
            print("creating field from existing particle")
            exported_field, fieldobject = field_from_particle(
                self.particle, **self.virtualsample_params, **kwargs
            )
            self.virtualsample_params["minimal_distance"] = (
                fieldobject.molecules_params["minimal_distance"]
            )
            if keep:
                self.exported_coordinate_field = exported_field
                self.objects_created["exported_coordinate_field"] = True
                self.coordinate_field = fieldobject
                self.objects_created["coordinate_field"] = True
            return exported_field
        else:
            # create minimal field
            fieldobject = coordinates_field.create_min_field(
                **self.virtualsample_params, **kwargs
            )
            exported_field = fieldobject.export_field()
            if keep:
                self.exported_coordinate_field = exported_field
                self.objects_created["exported_coordinate_field"] = True
                self.coordinate_field = fieldobject
                self.objects_created["coordinate_field"] = True

    def _build_imager(self, use_local_field=False, prints=True):
        """
        Build the imager object for the experiment using the selected imaging modalities.

        Parameters
        ----------
        :param use_local_field : bool, optional
            If True, use the local exported coordinate field. Default is False.
        :param prints : bool, optional
            If True, print status messages. Default is True.

        Returns
        -------
        None
        """
        if self.imaging_modalities:
            # print(f"Using selected mods: {self.imaging_modalities.keys()}")
            mods_list = list(self.imaging_modalities.keys())
            if use_local_field and self.generators_status(
                "exported_coordinate_field"
            ):
                self.imager, modality_parameters = create_imaging_system(
                    exported_field=self.exported_coordinate_field,
                    modalities_id_list=mods_list,
                    mod_params=self.imaging_modalities,
                    config_dir=self.configuration_path,
                    fluorophore_parameters=self.fluorophore_parameters,
                )
            else:
                if prints:
                    print(
                        "Local field missing or unused. Creating imager without particles"
                    )
                self.imager, modality_parameters = create_imaging_system(
                    modalities_id_list=mods_list,
                    mod_params=self.imaging_modalities,
                    config_dir=self.configuration_path,
                    fluorophore_parameters=self.fluorophore_parameters
                )
            # update the base acquisition parameters to consider all channels
            imager_channels = []
            anymod = list(self.imager.modalities.keys())[0]
            for chann in self.imager.modalities[anymod]["filters"].keys():
                imager_channels.append(chann)
            nchannels = len(imager_channels)
            self.reset_to_defaults(module="acquisitions", channels=imager_channels)
            self.objects_created["imager"] = True
        else:
            self.imager = None
            print("No modalities")

    def _gen_modality_noise_images(self):
        for modality_name in self.imager.modalities.keys():
            self.modality_noise_images[modality_name], _beads, _im_noiseless, _b_noiseless = self.imager.generate_imaging(
                modality=modality_name,
                noise=True, exp_time=0)

    def generators_status(self, generator_name):
        """
        Return the status of the specified generator.

        Parameters
        ----------
        :param generator_name : str
            The name of the generator whose status is to be retrieved.

        Returns
        -------
        object
            The status or object associated with the given generator name from the objects_created dictionary.

        Raises
        ------
        KeyError
            If the specified generator_name does not exist in objects_created.
        """
        return self.objects_created[generator_name]

    def build(
        self,
        use_locals=True,
        modules: list = [
            "all",
        ],
        **kwargs,
    ):
        """
        Build various simulation objects based on the specified modules.

        Parameters
        ----------
        :param use_locals : bool, optional
            Determines whether to use local variables or attributes when building objects. Default is True.
        :param modules : list, optional
            List of module names to build. If "all" is included, builds all available modules: "structure", "particle", "coordinate_field", and "imager". Default is ["all"].
        :param **kwargs
            Additional keyword arguments for building modules.

        Notes
        -----
        The order of building is: structure, particle, coordinate_field, imager.
        """
        print("Building objects")
        if "all" in modules:
            build_list = [
                "structure",
                "particle",
                "coordinate_field",
                "imager",
            ]
        else:
            build_list = modules
        if "structure" in build_list:
            self._build_structure()
        if "particle" in build_list:
            self._build_particle(keep=use_locals)
        if "coordinate_field" in build_list:
            self._build_coordinate_field(
                use_self_particle=use_locals, keep=use_locals
            )
        if "imager" in build_list:
            self._build_imager(use_local_field=use_locals)

    def _gen_reference(
        self, write=False, keep=False, ref_acq_pars=None, modality_wise=False
    ):
        """
        Calculate a reference image of the virtual sample by using the ideal parameters for each of the parameters to sweep.

        Requires the dictionary of the params2sweep. If modality_wise is True (default), a reference is calculated for each modality with that modality's parameters. Otherwise, it will use the Reference modality configuration to generate an idealised image as reference.

        Parameters
        ----------
        :param write : bool, optional
            Whether to write the reference image to disk. Default is False.
        :param keep : bool, optional
            If True, store the reference in the experiment. Default is False.
        :param ref_acq_pars : dict or None, optional
            Acquisition parameters for the reference. Default is None.
        :param modality_wise : bool, optional
            If True, calculate reference for each modality. Default is False.

        Returns
        -------
        tuple
            (_reference, _reference_parameters)
        """
        reference_pars = dict()
        output_name = "REFERENCE_"
        # get ideal parameters for virtual sample
        for param_name, param_settings in self.sweep_pars.items():
            reference_pars[param_name] = param_settings["ideal"]
        print(reference_pars)
        # generate particle
        tmp_particle = self._build_particle(
            lab_eff=reference_pars["labelling_efficiency"], keep=True
        )
        # use particle to create new field
        tmp_exported_field = self._build_coordinate_field()
        # create ideal image modality
        if modality_wise:
            self.imager.import_field(**tmp_exported_field)
            _reference = generate_multi_imaging_modalities(
                image_generator=self.imager,
                experiment_name=output_name,
                savingdir=self.output_directory,
                write=write,
            )
        else:
            reference_imager, ref_modality_parameters = create_imaging_system(
                modalities_id_list=["Reference"],
                config_dir=self.configuration_path,
                fluorophore_parameters=self.fluorophore_parameters
            )
            reference_imager.import_field(**tmp_exported_field)
            # make a copy
            _reference = dict()
            _reference_parameters = dict()
            for mod_name in list(self.imager.modalities.keys()):
                _reference_iteration = generate_multi_imaging_modalities(
                    image_generator=reference_imager,
                    experiment_name=output_name,
                    acquisition_param=ref_acq_pars,
                    savingdir=self.output_directory,
                    write=write,
                )
                _reference[mod_name] = _reference_iteration["Reference"]
                imager_scale = reference_imager.roi_params["scale"]
                scalefactor = np.ceil(
                    imager_scale / 1e-9
                )  # resulting pixel size in nanometers
                _reference_parameters[mod_name] = dict(
                    ref_pixelsize=reference_imager.modalities["Reference"][
                        "detector"
                    ]["pixelsize"]
                    * scalefactor
                )
        if keep:
            self.experiment_reference = _reference
            self.objects_created["output_reference"] = True
        return _reference, _reference_parameters

    def clear_labelled_structure(self):
        self.remove_probes()

    def clear_virtual_sample(self):
        if self.generators_status("coordinate_field"):
            self.exported_coordinate_field = None
            self.objects_created["exported_coordinate_field"] = False
        print("Virtual sample cleared")

    def create_example_experiment(self):
        with io.capture_output() as captured:
            self.select_structure("1XI5", build=False)
            self.add_probe()
            self.set_virtualsample_params()
            for modality_name in self.example_modalities:
                self.add_modality(modality_name)
            self.build()
        print("Experiment created with default parameters")

    def clear_experiment(self):
        """
        Clear the current experiment by resetting all parameters and objects.
        """
        self.clear_structure()
        self.clear_labelled_structure()
        self.clear_virtual_sample()
        self.clear_modalities()
        self.results = dict()

    def run_simulation(
        self,
        name="vlab4mic_experiment",
        acq_params=None,
        save=False,
        modality="All",
        **kwargs,
    ):
        """
        Run a simulation for the specified imaging modality or all modalities.

        Parameters
        ----------
        :param name : str, optional
            Name of the experiment or simulation. Default is "vlab4mic_experiment".
        :param acq_params : dict or None, optional
            Acquisition parameters for the simulation. If None and modality is "All", uses self.selected_mods.
        :param save : bool, optional
            Whether to save the simulation output to disk. Default is False.
        :param modality : str, optional
            The imaging modality to simulate. If "All", simulates all available modalities. Default is "All".
        :param **kwargs
            Additional keyword arguments passed to the imaging generator.

        Returns
        -------
        dict or None
            Simulation output as a dictionary mapping modality names to results. Returns None if the imager could not be created.
        """
        if not self.generators_status("imager"):
            self.build(modules="imager")
            if self.imager is None:
                print("Imager not created. Cannot run simulation.")
                return None
        if modality == "All":
            print("Simulating all modalities")
            if acq_params is None:
                acq_params = self.selected_mods
            if self.experiment_id:
                name = self.experiment_id
            simulation_output, simulation_output_noiseless = generate_multi_imaging_modalities(
                image_generator=self.imager,
                experiment_name=name,
                savingdir=self.output_directory,
                write=save,
                acquisition_param=acq_params,
            )
            self.results = simulation_output
            return simulation_output, simulation_output_noiseless
        else:
            print(f"Simulating: {modality}")
            acq_p = self.selected_mods[modality]
            timeseries_noise, beads_noise, timeseries_noiseless, beads_noiseless = self.imager.generate_imaging(
                modality=modality, **acq_p
            )
            simulation_output = {}
            simulation_output_noiseless = {}
            simulation_output[modality] = timeseries_noise
            simulation_output_noiseless[modality] = timeseries_noiseless
            return simulation_output, simulation_output_noiseless

    def remove_probes(self):
        """
        Remove all probe-related parameters and update the internal state accordingly.

        This method clears the probe parameters, updates the probes, and resets the
        particle and structure objects if their respective generators are active.
        It also clears any labels associated with the structure.

        Returns
        -------
        None
        """
        self.probe_parameters = dict()
        self._update_probes()
        if self.generators_status("particle"):
            self.particle = None
            self.objects_created["particle"] = False
            print("Labelled structure cleared")
        if self.generators_status("structure"):
            self.structure._clear_labels()
        print("Probes removed from experiment")

    def add_probe(
        self,
        probe_template: str = "NHS_ester",
        probe_name: str = None,
        probe_target_type: str = None,
        probe_target_value: str = None,
        probe_target_option: str = None,
        probe_distance_to_epitope: float = None,
        probe_model: str = None,
        probe_fluorophore: str = "AF647",
        probe_steric_hindrance=None,
        probe_paratope: str = None,
        probe_conjugation_target_info=None,
        probe_DoL: float = None,
        probe_seconday_epitope=None,
        probe_wobble_theta: float = None,
        labelling_efficiency: float = 1.0,
        as_primary=False,
        peptide_motif: dict = None,
        fluorophore_parameters = None,
        **kwargs,
    ):
        """
        Add a probe configuration to the experiment.

        This method loads a probe configuration from a YAML file, updates its parameters based on the provided arguments,
        and stores it in the experiment's probe parameters. It supports customization of probe properties such as target type,
        fluorophore, model, conjugation details, and more.

        Parameters
        ----------
        :param probe_name : str, optional
            Name of the probe. Default is "NHS_ester".
        :param probe_target_type : str, optional
            Type of the probe target (e.g., "Primary", "Secondary").
        :param probe_target_value : str, optional
            Value or identifier of the probe target.
        :param probe_target_option : str, optional
            Additional option for the probe target, used for secondary epitopes.
        :param probe_distance_to_epitope : float, optional
            Distance from the probe to the epitope.
        :param probe_model : str, optional
            List of model identifiers for the probe.
        :param probe_fluorophore : str, optional
            Identifier for the fluorophore. Default is "AF647".
        :param probe_steric_hindrance : Any, optional
            Steric hindrance value or configuration.
        :param probe_paratope : str, optional
            Paratope identifier or information.
        :param probe_conjugation_target_info : Any, optional
            Information about the conjugation target.
        :param probe_DoL : float, optional
            Efficiency of the probe conjugation.
        :param probe_seconday_epitope : Any, optional
            Information about a secondary epitope target.
        :param probe_wobbling : bool, optional
            Whether to enable probe wobbling. Default is False.
        :param labelling_efficiency : float, optional
            Efficiency of probe labelling. Default is 1.0.
        :param as_primary : bool, optional
            Whether to treat the probe as a primary linker. Default is False.
        :param **kwargs
            Additional keyword arguments for future extensions.

        Raises
        ------
        FileNotFoundError
            If the probe configuration YAML file does not exist.
        KeyError
            If required keys are missing in the configuration or probe parameters.

        Notes
        -----
        Updates the ``probe_parameters`` attribute with the new or modified probe configuration and calls the :meth:`_update_probes` method to refresh internal probe state.
        """
        probe_configuration = copy.deepcopy(
            self.config_probe_params[probe_template]
        )
        if probe_name is None:
            probe_name = probe_template
        else:
            probe_configuration["label_name"] = probe_name
        if peptide_motif is not None:
            protein_name, _1, site, sequence = (
                self.structure.get_peptide_motif(**peptide_motif)
            )
            if len(sequence) > 0:
                probe_target_type = "Sequence"
                probe_target_value = sequence
        if probe_target_type and probe_target_value:
            probe_configuration["target"] = dict(
                type=probe_target_type, value=probe_target_value
            )
            if probe_target_type == "Primary" and probe_target_option:
                if probe_target_value in self.probe_parameters.keys():
                    print(
                        f"Using {probe_target_option} as epitope on {probe_target_value}"
                    )
                    self.probe_parameters[probe_target_value][
                        "probe_seconday_epitope"
                    ] = probe_target_option
        elif (
            probe_configuration["target"]["type"] is None
            or probe_configuration["target"]["value"] is None
        ):
            print(
                "No target info provided for the probe. Retrieving random sequence."
            )
            # probe has no target info
            # a random target will be used
            protein_name, _1, site, sequence = (
                self.structure.get_peptide_motif(position="cterminal")
            )
            probe_configuration["target"] = dict(
                type="Sequence", value=sequence
            )
            # probe_configuration["target"]["type"] = "Sequence"
            # probe_configuration["target"]["value"] = sequence
        if probe_distance_to_epitope is not None:
            probe_configuration["distance_to_epitope"] = (
                probe_distance_to_epitope
            )
        if probe_fluorophore is not None:
            probe_configuration["fluorophore_id"] = probe_fluorophore
        if fluorophore_parameters is not None:
                self.add_fluorophore_parameters(fluorophoe_id=probe_fluorophore, **fluorophore_parameters)
        if labelling_efficiency is not None:
            probe_configuration["labelling_efficiency"] = labelling_efficiency
        if probe_model is not None:
            probe_configuration["model_ID"] = probe_model
        if probe_paratope is not None:
            probe_configuration["paratope"] = probe_paratope
        if probe_conjugation_target_info is not None:
            probe_configuration["conjugation_target_info"] = (
                probe_conjugation_target_info
            )
            probe_configuration["conjugation_sites"]["target"]["type"] = (
                probe_conjugation_target_info["type"]
            )
            probe_configuration["conjugation_sites"]["target"]["value"] = (
                probe_conjugation_target_info["value"]
            )
        if probe_DoL is not None:
            print(f"Setting probe DoL to: {probe_DoL}")
            probe_configuration["conjugation_sites"]["DoL"] = (
                probe_DoL
            )
        else:
            print("Add_probe: Using default probe DoL from template")
        if probe_seconday_epitope is not None:
            probe_configuration["epitope_target_info"] = probe_seconday_epitope
        if probe_wobble_theta is not None:
            probe_configuration["enable_wobble"] = True
            probe_configuration["wobble_theta"] = probe_wobble_theta
        if as_primary:
            print("Adding probe as primary linker")
            probe_configuration["as_linker"] = True
        else:
            probe_configuration["as_linker"] = False
        if probe_steric_hindrance is not None:
            probe_configuration["distance_between_epitope"] = (
                probe_steric_hindrance
            )
            probe_configuration["binding"]["distance"][
                "between_targets"
            ] = probe_steric_hindrance
        self.probe_parameters[probe_name] = probe_configuration
        self._update_probes()

    def _update_probes(self):
        """
        Update the structure_label attribute based on the current probe_parameters.

        Returns
        -------
        None
        """
        if len(self.probe_parameters.keys()) == 0:
            self.structure_label = None
        else:
            self.structure_label = list(self.probe_parameters.keys())

    def add_fluorophore_parameters(self,
                                fluorophoe_id: str = None,
                                emission = None,
                                blinking_rates = None,
                                **kwargs):
        if fluorophoe_id not in self.fluorophore_parameters.keys():
            self.fluorophore_parameters[fluorophoe_id] = dict()
            self.fluorophore_parameters[fluorophoe_id]["emission"] = dict()
            self.fluorophore_parameters[fluorophoe_id]["emission"]["type"] = "constant"
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"] = dict()
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"]["kon"] = 0.99
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"]["koff"] = 0.01
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"]["kbleach"] = 0.0
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"]["initial_state"] = 1
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"]["photons_per_blink"] = 1000
        if emission is not None:
            self.fluorophore_parameters[fluorophoe_id]["emission"] = copy.deepcopy(emission)
        if blinking_rates is not None:
            self.fluorophore_parameters[fluorophoe_id]["blinking_rates"] = blinking_rates

    def set_virtualsample_params(
        self,
        virtualsample_template="square1x1um_randomised",
        sample_dimensions: list[int, int, int] = None,
        number_of_particles: int = None,
        particle_positions: list = None,
        random_orientations: bool = None,
        random_placing: bool = None,
        minimal_distance: float = None,
        update_mode: bool = True,
        random_rotations=False,
        rotation_angles=None,
        xy_orientations = None,
        xz_orientations = None,
        yz_orientations = None,
        axial_offset=None,
        **kwargs,
    ):
        """
        Set parameters for the virtual sample by loading a template configuration and optionally overriding specific values.

        Parameters
        ----------
        :param virtualsample_template : str, optional
            Name of the virtual sample template YAML file (without extension) to load. Default is "square1x1um_randomised".
        :param sample_dimensions : list of int, optional
            Dimensions of the sample to override the template value.
        :param number_of_particles : int, optional
            Number of particles to override the template value.
        :param particle_positions : list, optional
            List of particle positions to override the template value.
        :param random_orientations : bool, optional
            Whether to randomize particle orientations, overrides the template value.
        :param random_placing : bool, optional
            Whether to randomize particle placement, overrides the template value.
        :param **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        None
        """
        # load default configuration for virtual sample
        try:
            particle_minimal_distance = self.coordinate_field.molecules_params[
                "minimal_distance"
            ]
        except:
            particle_minimal_distance = None
        if update_mode:
            pass
        else:
            virtual_sample_template = os.path.join(
                self.configuration_path,
                "virtualsample",
                virtualsample_template + ".yaml",
            )
            self.virtualsample_params = load_yaml(virtual_sample_template)
        if sample_dimensions is not None:
            self.virtualsample_params["sample_dimensions"] = sample_dimensions
        if number_of_particles is not None:
            self.virtualsample_params["number_of_particles"] = number_of_particles
        if particle_positions is not None:
            self.virtualsample_params["relative_positions"] = particle_positions
        if random_orientations is not None:
            self.virtualsample_params["random_orientations"] = random_orientations
        if random_placing is not None:
            self.virtualsample_params["random_placing"] = random_placing
        if minimal_distance is not None:
            self.virtualsample_params["minimal_distance"] = minimal_distance
        else:
            self.virtualsample_params["minimal_distance"] = (
                particle_minimal_distance
            )
        self.virtualsample_params["random_rotations"] = random_rotations
        if rotation_angles is not None:
            self.virtualsample_params["rotation_angles"] = rotation_angles
        if xy_orientations is not None:
            self.virtualsample_params["xy_orientations"] = xy_orientations
        if xz_orientations is not None:
            self.virtualsample_params["xz_orientations"] = xz_orientations
        if yz_orientations is not None:
            self.virtualsample_params["yz_orientations"] = yz_orientations
        if axial_offset is not None:
            self.virtualsample_params["axial_offset"] = axial_offset

    def use_image_for_positioning(
        self,
        img,
        mode="localmaxima",
        sigma=None,
        background=None,
        threshold=None,
        pixelsize=None,
        min_distance=None,
        **kwargs,
    ):
        """
        Generate and set relative positions for a virtual sample based on features detected in an input image.

        This method processes the provided image to extract feature coordinates using the specified detection mode
        and parameters. The resulting positions are stored in `self.virtualsample_params["relative_positions"]`,
        and the sample's physical dimensions are updated accordingly. The method then triggers the build process
        for the 'coordinate_field' and 'imager' modules.

        Parameters
        ----------
        :param img : array-like
            The input image (as a NumPy array or compatible format) to be analyzed for feature detection.
        :param mode : str, optional
            The feature detection mode to use (e.g., "localmaxima"). Default is "localmaxima".
        :param sigma : float or None, optional
            Standard deviation for Gaussian smoothing. If None, no smoothing is applied.
        :param background : float or None, optional
            Background intensity to subtract from the image. If None, no subtraction.
        :param threshold : float or None, optional
            Minimum intensity threshold for feature detection. If None, uses default.
        :param pixelsize : float or None, optional
            Pixel size of the input image. If None, uses default.
        :param min_distance : float or None, optional
            Minimum distance between detected features. If None, uses default.
        :param **kwargs
            Additional keyword arguments passed to the feature detection function.

        Returns
        -------
        None

        Notes
        -----
        Updates `self.virtualsample_params` with new relative positions and sample dimensions. Calls `self.build()` for the specified modules.
        """
        xyz_relative, image_physical_size = (
            coordinates_field.gen_positions_from_image(
                img=img,
                mode=mode,
                sigma=sigma,
                background=background,
                threshold=threshold,
                pixelsize=pixelsize,
                min_distance=min_distance,
                **kwargs,
            )
        )
        self.set_virtualsample_params(
            sample_dimensions=[
                image_physical_size[0],
                image_physical_size[1],
                100,
            ],
            particle_positions=xyz_relative,
            number_of_particles=len(xyz_relative),
            random_orientations=False,
            random_placing=False,
            minimal_distance=min_distance,
        )
        # self.virtualsample_params["relative_positions"] = xyz_relative
        # self.virtualsample_params["sample_dimensions"] = [
        #    image_physical_size[0],
        #    image_physical_size[1],
        #    100,
        # ]
        self.build(modules=["coordinate_field", "imager"])

    def current_settings(
        self, as_string=True, newline="<br>", modalities_acq_params=False
    ):
        """
        Print the current settings of the experiment, including structure, particle, coordinate field, and imaging modalities.

    Returns
    -------
    None
        """
        string = "Current settings of the experiment:" + newline
        string += f"Structure ID: {self.structure_id}" + newline
        string += f"Probes: {list(self.probe_parameters.keys())}" + newline
        string += f"Virtual sample: {self.coordinate_field_id}" + newline
        string += "Imaging Modalities: "
        for modality_name, acqparams in self.selected_mods.items():
            string += f"  {modality_name}"
            if modalities_acq_params:
                string += f": {acqparams}" + newline
            else:
                string += "; "
        string += newline
        if not as_string:
            print(string)
        else:
            return string

    def expand_virtual_sample(self, factor=1):
        if factor > 1:
            self.coordinate_field.expand_isotropically(factor=factor)
            self.build(modules=["imager",])

def generate_virtual_sample(
    structure: str = "1XI5",
    structure_is_path = False,
    probe_template: str = "NHS_ester",
    probe_name: str = None,
    probe_target_type: str = None,
    probe_target_value: str = None,
    probe_distance_to_epitope: float = None,
    probe_model: list[str] = None,
    probe_fluorophore: str = "AF647",
    probe_paratope: str = None,
    probe_conjugation_target_info=None,
    probe_DoL: float = None,
    probe_seconday_epitope=None,
    probe_wobble_theta=None,
    labelling_efficiency: float = 1.0,
    incomplete_labelling_small_cluster: float = None,
    incomplete_labelling_large_cluster: float = None,
    incomplete_labelling: float = None,
    virtual_sample_template: str = "square1x1um_randomised",
    sample_dimensions: list[float] = None,
    number_of_particles: int = None,
    particle_positions: list[np.array] = None,
    random_orientations=False,
    xy_orientations = None,
    xz_orientations = None,
    yz_orientations = None,
    axial_offset = None,
    random_placing=False,
    random_rotations=False,
    rotation_angles = None,
    expansion_factor=1,
    clear_probes=False,
    clear_experiment=False,
    # secondary
    primary_probe=None,
    secondary_probe=None,
    probe_list=None,
    **kwargs,
):
    """
    Create a virtual sample from a structure model.

    Parameters
    ----------
    :param structure : str, optional
        4-letter ID of PDB/CIF model. Default is "1XI5".
    :param structure_is_path : logical
        Use structure value as absolute path for the PDB/CIF file.
    :param probe_name : str, optional
        Name ID of probe configuration file (filename).
    :param probe_target_type : str, optional
        Options: "Sequence", "Atom_residue", or "Primary".
    :param probe_target_value : str or dict, optional
        For target type "Sequence" or "Primary", a string. For "Atom_residue", a dictionary with keys "Atom" and "Residue".
    :param probe_distance_to_epitope : float, optional
        Minimal distance set from epitope and probe paratope.
    :param probe_model : list of str, optional
        4-letter ID(s) of PDB/CIF model(s).
    :param probe_fluorophore : str, optional
        Fluorophore name (e.g., "AF647"). Default is "AF647".
    :param probe_paratope : str, optional
        Sequence of the paratope site for when probe includes a model.
    :param probe_conjugation_target_info : Any, optional
        Information about the probe conjugation target.
    :param probe_DoL : float, optional
        Efficiency of conjugation of emitters.
    :param probe_seconday_epitope : str, optional
        Sequence within probe model to be used as epitope for a secondary.
    :param probe_wobbling : bool, optional
        Enable probe wobbling. Default is False.
    :param labelling_efficiency : float, optional
        Labelling efficiency of probe. Default is 1.0.
    :param incomplete_labelling_small_cluster : float, optional
        In Å, distance used to group epitopes into multimers.
    :param incomplete_labelling_large_cluster : float, optional
        In Å, distance within multimers to consider neighbors.
    :param incomplete_labelling : float, optional
        Fraction of incomplete_labelling to model.
    :param virtual_sample_template : str, optional
        Name of the configuration file for template. Default is "square1x1um_randomised".
    :param sample_dimensions : list of float, optional
        In nanometers, define the X, Y, and Z sizes of the field.
    :param number_of_particles : int, optional
        Number of independent copies of a particle to create and distribute.
    :param particle_positions : list of np.array, optional
        Relative positions of particles in the field.
    :param random_orientations : bool, optional
        If True, each particle will be randomly assigned a new orientation. Default is False.
    :param random_placing : bool, optional
        Define if position in field is random or the center of field. Default is False.
    :param clear_probes : bool, optional
        If True, default parameters will be cleared. Default is False.
    :param **kwargs
        Additional keyword arguments.

    Returns
    -------
    tuple
        - dict: The virtual sample as exported format. This can be used as input for the `image_vsample` method.
        - ExperimentParametrisation: The experiment containing all modules that were generated to build the virtual sample, and the virtual sample module itself. This experiment can be further used and tweaked for subsequent analysis or branching workflows.
    """
    myexperiment = ExperimentParametrisation()
    if clear_experiment:
        myexperiment.clear_experiment()
    # select structure
    if structure_is_path:
        print(f"Selecting structure from path: {structure}")
        myexperiment.select_structure(
            structure_id=structure.split(".")[0], 
            structure_path=structure,
            build=True
        )
    else:
        print("Selecting structure from ID:", structure)
        myexperiment.select_structure(
            structure_id=structure, 
            structure_path=None,
            build=True
        )
    # load default configuration for probe
    if (primary_probe is not None) and (secondary_probe is not None):
        print("Adding primary and secondary probes")
        myexperiment.add_probe(as_primary=True, **primary_probe)
        myexperiment.add_probe(as_primary=False, **secondary_probe)
    elif probe_list is not None:
        for probe in probe_list:
            print(f"Adding probe: {probe['probe_name']}")
            myexperiment.add_probe(**probe)
    else:
        if not clear_probes:
            print(probe_template)
            probe_configuration_file = os.path.join(
                myexperiment.configuration_path, "probes", probe_template + ".yaml"
            )
            probe_configuration = load_yaml(probe_configuration_file)
            probe_configuration["probe_template"] = probe_template
            if probe_name is None:
                probe_name = probe_template
            else:
                probe_configuration["label_name"] = probe_name
            if probe_target_type and probe_target_value:
                print(probe_target_type, probe_target_value)
                probe_configuration["probe_target_type"] = probe_target_type
                probe_configuration["probe_target_value"] = probe_target_value
            if probe_distance_to_epitope is not None:
                probe_configuration["distance_to_epitope"] = (
                    probe_distance_to_epitope
                )
            if probe_fluorophore is not None:
                probe_configuration["fluorophore_id"] = probe_fluorophore
            if labelling_efficiency is not None:
                probe_configuration["labelling_efficiency"] = labelling_efficiency
            if probe_model is not None:
                probe_configuration["model_ID"] = probe_model
            if probe_paratope is not None:
                probe_configuration["paratope"] = probe_paratope
            if probe_conjugation_target_info is not None:
                probe_configuration["conjugation_target_info"] = (
                    probe_conjugation_target_info
                )
            if probe_DoL is not None:
                probe_configuration["probe_DoL"] = (
                    probe_DoL
                )
            if probe_seconday_epitope is not None:
                probe_configuration["epitope_target_info"] = probe_seconday_epitope
            if probe_wobble_theta is not None:
                probe_configuration["probe_wobble_theta"] = probe_wobble_theta
            myexperiment.add_probe(**probe_configuration)
    # load default configuration for virtual sample
    virtual_sample_template = os.path.join(
        myexperiment.configuration_path,
        "virtualsample",
        virtual_sample_template + ".yaml",
    )
    vsample_configuration = load_yaml(virtual_sample_template)
    #myexperiment.configuration_path
    if incomplete_labelling and incomplete_labelling_large_cluster and incomplete_labelling_small_cluster:
        myexperiment.incomplete_labelling_eps["eps1"] = incomplete_labelling_small_cluster
        myexperiment.incomplete_labelling_eps["eps2"] = incomplete_labelling_large_cluster
        myexperiment.incomplete_labelling_eps["incomplete_labelling"] = incomplete_labelling
        myexperiment.incomplete_labelling_eps["use_incomplete_labelling"] = True

    if sample_dimensions is not None:
        vsample_configuration["sample_dimensions"] = sample_dimensions
    if number_of_particles is not None:
        vsample_configuration["number_of_particles"] = number_of_particles
    if particle_positions is not None:
        vsample_configuration["relative_positions"] = particle_positions
    if random_orientations is not None:
        vsample_configuration["random_orientations"] = random_orientations
    if random_placing is not None:
        vsample_configuration["random_placing"] = random_placing
    if random_rotations:
        vsample_configuration["random_rotations"] = random_rotations
    vsample_configuration["rotation_angles"] = rotation_angles
    vsample_configuration["xy_orientations"] = xy_orientations
    vsample_configuration["xz_orientations"] = xz_orientations
    vsample_configuration["yz_orientations"] = yz_orientations
    vsample_configuration["axial_offset"] = axial_offset
    myexperiment.virtualsample_params = vsample_configuration
    myexperiment.build(modules=[
                "particle",
                "coordinate_field",
                "imager"], use_locals=True)
    if expansion_factor > 1:
        myexperiment.expand_virtual_sample(factor=expansion_factor)
    # myexperiment.coordinate_field_id = virtual_sample
    return myexperiment.exported_coordinate_field, myexperiment


def build_virtual_microscope(
    modality="STED", multimodal: list[str] = None, experiment=None, **kwargs
):
    """
    Initialise a virtual microscope for single or multimodal imaging.

    Parameters
    ----------
    :param modality : str, optional
        Modality name. Default is "STED".
    :param multimodal : list of str, optional
        List of modality names. Overrides the `modality` parameter if provided.
    :param experiment : ExperimentParametrisation, optional
        An Experiment object. If None, a new one is created.
    :param **kwargs
        Additional arguments passed to `add_modality`.

    Returns
    -------
    tuple
        - Imager: The virtual microscope with the specified modality models. Contains a default sample.
        - ExperimentParametrisation: The experiment containing the virtual microscope. All other modules are not initialised. This experiment can be further used and tweaked for subsequent analysis or branching workflows.
    """
    if experiment is None:
        experiment = ExperimentParametrisation()
    if multimodal is not None:
        for mod in multimodal:
            print(mod)
            experiment.add_modality(mod, **kwargs)
    else:
        experiment.add_modality(modality, **kwargs)
    if "use_local_field" in kwargs.keys():
        experiment._build_imager(use_local_field=kwargs["use_local_field"])
    else:
        experiment._build_imager()
    return experiment.imager, experiment


def image_vsample(
    vsample=None,
    modality: str = "STED",
    multimodal: list[str] = None,
    run_simulation: bool = True,
    structure: str = "1XI5",
    structure_is_path: bool = False,
    probe_template: str = "NHS_ester",
    probe_name: str = None,
    probe_target_type: str = None,
    probe_target_value: str = None,
    probe_distance_to_epitope: float = None,
    probe_model: list = None,
    probe_fluorophore: str = "AF647",
    probe_paratope: str = None,
    probe_conjugation_target_info = None,
    probe_DoL: float = None,
    probe_seconday_epitope = None,
    probe_wobble_theta = None,
    labelling_efficiency: float = 1.0,
    incomplete_labelling_small_cluster: float = None,
    incomplete_labelling_large_cluster: float = None,
    incomplete_labelling: float = None,
    virtual_sample_template: str = "square1x1um_randomised",
    sample_dimensions: list = None,
    number_of_particles: int = None,
    particle_positions: list = None,
    random_orientations = False,
    xy_orientations: list[int] = None,
    xz_orientations: list[int] = None,
    yz_orientations: list[int] = None,
    axial_offset: list[int] = None,
    random_placing = False,
    random_rotations = False,
    rotation_angles = None,
    expansion_factor = 1,
    clear_probes = False,
    clear_experiment = False,
    primary_probe = None,
    secondary_probe = None,
    probe_list = None,
    save: bool = False,
    **kwargs,
):
    """
    Generate imaging simulations of the specified virtual sample and imaging modalities.

    This function can either:
    1. Take an existing virtual sample (as a dictionary) and wrap it in a new experiment, adding the specified imaging modalities.
    2. Or, generate a new virtual sample and experiment from scratch using the provided parameters (matching those of `generate_virtual_sample`).

    Parameters
    ----------
    :param vsample : dict, optional
        Dictionary specifying sample parameters. Corresponds to Experiment attribute "exported_coordinate_field". If None, a new sample is generated.
    :param modality : str, optional
        Modality name to use for imaging. Default is "STED".
    :param multimodal : list of str, optional
        List of modality names. If provided, overrides the `modality` parameter and adds all listed modalities.
    :param run_simulation : bool, optional
        If True, generates image simulation for each modality set. Default is True.
    
    ### Parameters for virtual sample generation (see `generate_virtual_sample` for details):
    
    :param structure : str, optional
        4-letter ID of PDB/CIF model. Default is "1XI5".
    :param structure_is_path : bool, optional
        Use structure value as absolute path for the PDB/CIF file.
    :param probe_template : str, optional
        Name of probe configuration file (filename). Default is "NHS_ester".
    :param probe_name : str, optional
        Name for the probe configuration.
    :param probe_target_type : str, optional
        Options: "Sequence", "Atom_residue", or "Primary".
    :param probe_target_value : str or dict, optional
        For target type "Sequence" or "Primary", a string. For "Atom_residue", a dictionary with keys "Atom" and "Residue".
    :param probe_distance_to_epitope : float, optional
        Minimal distance set from epitope and probe paratope.
    :param probe_model : list, optional
        4-letter ID(s) of PDB/CIF model(s).
    :param probe_fluorophore : str, optional
        Fluorophore name (e.g., "AF647"). Default is "AF647".
    :param probe_paratope : str, optional
        Sequence of the paratope site for when probe includes a model.
    :param probe_conjugation_target_info : any, optional
        Information about the probe conjugation target.
    :param probe_DoL : float, optional
        Efficiency of conjugation of emitters.
    :param probe_seconday_epitope : str, optional
        Sequence within probe model to be used as epitope for a secondary.
    :param probe_wobble_theta : any, optional
        Enable probe wobbling.
    :param labelling_efficiency : float, optional
        Labelling efficiency of probe. Default is 1.0.
    :param incomplete_labelling_small_cluster : float, optional
        In Å, distance used to group epitopes into multimers.
    :param incomplete_labelling_large_cluster : float, optional
        In Å, distance within multimers to consider neighbors.
    :param incomplete_labelling : float, optional
        Fraction of incomplete_labelling to model.
    :param virtual_sample_template : str, optional
        Name of the configuration file for template. Default is "square1x1um_randomised".
    :param sample_dimensions : list, optional
        In nanometers, define the X, Y, and Z sizes of the field.
    :param number_of_particles : int, optional
        Number of independent copies of a particle to create and distribute.
    :param particle_positions : list, optional
        Relative positions of particles in the field.
    :param random_orientations : bool, optional
        If True, each particle will be randomly assigned a new orientation. Default is False.
    :param xy_orientations, xz_orientations, yz_orientations : any, optional
        Orientation parameters for the sample.
    :param axial_offset : any, optional
        Axial offset for the sample.
    :param random_placing : bool, optional
        Define if position in field is random or the center of field. Default is False.
    :param random_rotations : bool, optional
        If True, apply random rotations to particles. Default is False.
    :param rotation_angles : any, optional
        Rotation angles for the sample.
    :param clear_probes : bool, optional
        If True, default probe parameters will be cleared. Default is False.
    :param clear_experiment : bool, optional
        If True, clear the experiment before generating a new sample. Default is False.
    :param primary_probe, secondary_probe : any, optional
        Dictionaries for primary and secondary probe configuration.

    Returns
    -------
    tuple
        - dict: Image simulations (per modality).
        - dict: Noiseless image simulations (per modality).
        - ExperimentParametrisation: The experiment object containing all modules and configuration.
    """
    if vsample is None:
        vsample, sample_experiment = generate_virtual_sample(
            structure=structure,
            structure_is_path=structure_is_path,
            probe_template=probe_template,
            probe_name=probe_name,
            probe_target_type=probe_target_type,
            probe_target_value=probe_target_value,
            probe_distance_to_epitope=probe_distance_to_epitope,
            probe_model=probe_model,
            probe_fluorophore=probe_fluorophore,
            probe_paratope=probe_paratope,
            probe_conjugation_target_info=probe_conjugation_target_info,
            probe_DoL=probe_DoL,
            probe_seconday_epitope=probe_seconday_epitope,
            probe_wobble_theta=probe_wobble_theta,
            labelling_efficiency=labelling_efficiency,
            incomplete_labelling_small_cluster=incomplete_labelling_small_cluster,
            incomplete_labelling_large_cluster=incomplete_labelling_large_cluster,
            incomplete_labelling=incomplete_labelling,
            virtual_sample_template=virtual_sample_template,
            sample_dimensions=sample_dimensions,
            number_of_particles=number_of_particles,
            particle_positions=particle_positions,
            random_orientations=random_orientations,
            xy_orientations=xy_orientations,
            xz_orientations=xz_orientations,
            yz_orientations=yz_orientations,
            axial_offset=axial_offset,
            random_placing=random_placing,
            random_rotations=random_rotations,
            rotation_angles=rotation_angles,
            clear_probes=clear_probes,
            clear_experiment=clear_experiment,
            primary_probe=primary_probe,
            secondary_probe=secondary_probe,
            probe_list=probe_list,
            expansion_factor=expansion_factor,
        )
        sample_experiment.clear_modalities()
        if multimodal is not None:
            for mod in multimodal:
                try:
                    mod_template = kwargs[mod]["modality_template"]
                    print(f"Adding modality: {mod}, {mod_template}")
                    sample_experiment.add_modality(mod, modality_template=mod_template)
                except:
                    sample_experiment.add_modality(mod)
        else:
            sample_experiment.add_modality(modality)
        sample_experiment.build(
            modules=[
                "imager",
            ]
        )
        for modality_name in sample_experiment.imaging_modalities.keys():
            if modality_name in kwargs.keys():
                sample_experiment.update_modality(
                    modality_name, **kwargs[modality_name]
                )
                sample_experiment.set_modality_acq(
                    modality_name, **kwargs[modality_name]
                )
    else:
        # need to create an experiment for it
        if multimodal is not None:
            vmicroscope, sample_experiment = build_virtual_microscope(
                multimodal=multimodal, use_local_field=True
            )
        else:
            vmicroscope, sample_experiment = build_virtual_microscope(
                modality=modality, use_local_field=True
            )
        for modality_name in sample_experiment.imaging_modalities.keys():
            if modality_name in kwargs.keys():
                sample_experiment.update_modality(
                    modality_name, **kwargs[modality_name]
                )
                sample_experiment.set_modality_acq(
                    modality_name, **kwargs[modality_name]
                )
        sample_experiment.imager.import_field(**vsample)
    imaging_output = dict()
    imaging_output_noiseless = dict()
    if run_simulation:
        imaging_output, imaging_output_noiseless = sample_experiment.run_simulation(save=save)
    return imaging_output, imaging_output_noiseless, sample_experiment
