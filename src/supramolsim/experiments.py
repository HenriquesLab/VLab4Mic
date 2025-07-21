from dataclasses import dataclass, field, fields
from typing import List, Dict
import numpy as np
import supramolsim
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
from supramolsim.utils.io.yaml_functions import load_yaml
import numpy as np
import os
import copy
from supramolsim.utils.io import yaml_functions
from IPython.utils import io
from pathlib import Path

output_path = Path.home() / "vlab4mic_outputs"

if not os.path.exists(output_path):
    os.makedirs(output_path)


@dataclass
class ExperimentParametrisation:
    experiment_id: str = "vLab4mic_experiment"
    structure_id: str = "1XI5"
    configuration_path: str = ""
    structure_label: str = "NHS_ester"
    fluorophore_id: str = "AF647"
    coordinate_field_id: str = "square1x1um_randomised"
    selected_mods: Dict[str, int] = field(default_factory=dict)
    imaging_modalities: Dict[str, int] = field(default_factory=dict)
    probe_parameters: Dict[str, int] = field(default_factory=dict)
    defect_eps: Dict[str, int] = field(default_factory=dict)
    sweep_pars: Dict[str, int] = field(default_factory=dict)
    objects_created: Dict[str, int] = field(default_factory=dict)
    output_directory: str = str(output_path)
    example_structures = ["3J3Y", "7R5K", "1XI5", "8GMO"]
    example_modalities = ["Widefield", "Confocal", "STED", "SMLM"]

    def __post_init__(self):
        pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
        local_dir = os.path.join(pck_dir, "configs")
        self.configuration_path = local_dir
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
        self.defect_eps["defect"] = None
        self.defect_eps["eps1"] = None
        self.defect_eps["eps2"] = None
        self.defect_eps["use_defects"] = False
        # read information of local modalities configuration
        modalities_dir = os.path.join(local_dir, "modalities")
        modalities_names_list = []
        modality_parameters = {}
        for mods in os.listdir(modalities_dir):
            if os.path.splitext(mods)[-1] == ".yaml" and "_template" not in mods:
                modalities_names_list.append(os.path.splitext(mods)[0])
        for mod in modalities_names_list:
            mod_info = configuration_format.compile_modality_parameters(mod, local_dir)
            modality_parameters[mod] = mod_info
        self.local_modalities_names = modalities_names_list
        self.local_modalities_parameters = modality_parameters
        probes_dir = os.path.join(local_dir, "probes")
        structure_dir = os.path.join(local_dir, "structures")
        self.config_probe_params = {}
        self.config_probe_models_names = []
        self.config_global_probes_names = []
        self.config_probe_per_structure_names = {}
        for p_file in os.listdir(probes_dir):
            if os.path.splitext(p_file)[-1] == ".yaml" and "_template" not in p_file:
                label_config_path = os.path.join(probes_dir, p_file)
                label_parmeters = supramolsim.load_yaml(label_config_path)
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
                        if struct in self.config_probe_per_structure_names.keys():
                            self.config_probe_per_structure_names[struct].append(
                                lablname
                            )
                        else:
                            self.config_probe_per_structure_names[struct] = [
                                lablname,
                            ]
        self.demo_structures = []
        # get available structure IDs
        self.structures_info_list = dict()
        structure_dir = os.path.join(self.configuration_path, "structures")
        fluorophores_dir = os.path.join(self.configuration_path, "fluorophores")
        probes_dir = os.path.join(self.configuration_path, "probes")
        modalities_dir = os.path.join(self.configuration_path, "modalities")
        for file in os.listdir(structure_dir):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
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

    def select_structure(self, structure_id="1XI5", build=True):
        """
        Select a molecular structure by its identifier and optionally build the structure module.

        Parameters
        ----------
        structure_id : str, optional
            The identifier of the structure to select. Default is "1XI5".
        build : bool, optional
            If True, triggers the build process for the structure module. Default is True.

        Returns
        -------
        None
        """
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
        if self.generators_status("structure"):
            self.structure = None
        self.objects_created["structure"] = False
        print("Structure cleared")

    def add_modality(self, modality_name, save=False, **kwargs):
        """
        Add a new imaging modality to the experiment.

        If the specified modality name exists in the local modalities, it copies the parameters from the local template.
        Otherwise, it uses the 'Widefield' modality parameters as a template to create a new modality.
        Additional parameters for the modality can be specified via keyword arguments and will override the defaults.
        Optionally, the modality output can be saved by setting `save` to True.

        Parameters
        ----------
        modality_name : str
            The name of the modality to add.
        save : bool, optional
            If True, saves the modality output. Default is False.
        **kwargs
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
                f"Modality {modality_name} not found in demo modalities. Using Widefield params as template to create new."
            )
            self.imaging_modalities[modality_name] = copy.deepcopy(
                self.local_modalities_parameters["Widefield"]
            )
            for param, value in kwargs.items():
                self.imaging_modalities[modality_name][param] = value
            self.set_modality_acq(modality_name, save=save)

    def update_modality(
        self,
        modality_name,
        pixelsize_nm: int = None,
        lateral_resolution_nm: int = None,
        axial_resolution_nm: int = None,
        psf_voxel_nm: int = None,
        remove=False,
    ):
        """
        Update or remove an imaging modality's parameters.

        This method allows updating specific parameters of an imaging modality, such as pixel size, lateral and axial resolution, and PSF voxel size.
        If `remove` is True, the modality is removed from the internal dictionaries.

        Parameters
        ----------
        modality_name : str
            The name of the imaging modality to update or remove.
        pixelsize_nm : int, optional
            The new pixel size in nanometers. If provided, updates the detector pixel size.
        lateral_resolution_nm : int, optional
            The new lateral resolution in nanometers. If provided, updates the lateral standard deviations of the PSF.
        axial_resolution_nm : int, optional
            The new axial resolution in nanometers. If provided, updates the axial standard deviation of the PSF.
        psf_voxel_nm : int, optional
            The new PSF voxel size in nanometers. If provided, updates the PSF voxel size for all axes.
        remove : bool, optional
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
                self.imaging_modalities[modality_name]["detector"]["pixelsize"] = (
                    pixelsize_nm / 1000
                )
                changes = True
            if lateral_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name]["psf_params"][
                    "voxelsize"
                ][0]
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][0] = (
                    lateral_resolution_nm / voxel_size
                )
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][1] = (
                    lateral_resolution_nm / voxel_size
                )
                changes = True
            if axial_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name]["psf_params"][
                    "voxelsize"
                ][0]
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][2] = (
                    axial_resolution_nm / voxel_size
                )
                changes = True
            if psf_voxel_nm is not None:
                self.imaging_modalities[modality_name]["psf_params"]["voxelsize"] = [
                    psf_voxel_nm,
                    psf_voxel_nm,
                    psf_voxel_nm,
                ]
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
        modality_name : str
            The name of the imaging modality to configure.
        exp_time : float, optional
            Exposure time for the acquisition in seconds. Default is 0.001.
        noise : bool, optional
            Whether to include noise in the acquisition. Default is True.
        save : bool, optional
            Whether to save the acquired data. Default is False.
        nframes : int, optional
            Number of frames to acquire. Default is 1.
        channels : list of str, optional
            List of channel names to use. Default is ["ch0"].
        **kwargs
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
            self.update_modality(
                modality_name, remove=True
            )
        if self.generators_status("imager"):
            self.imager = None
            self.objects_created["imager"] = False
        print("All modalities cleared")

    def reset_to_defaults(self, module="acquisitions", **kwargs):
        """
        Reset acquisition parameters or other module settings to their default values.

        Parameters
        ----------
        module : str, optional
            The module to reset. Default is "acquisitions".
        **kwargs
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
        keep : bool, optional
            If True, store the structure in the experiment. Default is True.

        Returns
        -------
        None
        """
        if self.structure_id:
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
        lab_eff : float, optional
            Labelling efficiency. Default is 1.
        keep : bool, optional
            If True, keep the label. Default is False.
        **kwargs
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

    def _check_if_defects(self):
        """
        Check if defect parameters are set and update the use_defects flag.

        Returns
        -------
        None
        """
        if (
            self.defect_eps["defect"]
            and self.defect_eps["eps1"]
            and self.defect_eps["eps2"]
        ):
            self.defect_eps["use_defects"] = True

    def _build_particle(self, lab_eff=1.0, defect_build=None, keep=False):
        """
        Build the particle object for the experiment, optionally adding defects.

        Parameters
        ----------
        lab_eff : float, optional
            Labelling efficiency. Default is 1.0.
        defect_build : float or None, optional
            Defect parameter to use. Default is None.
        keep : bool, optional
            If True, store the particle in the experiment. Default is False.

        Returns
        -------
        particle or None
            The particle object if created, else None.
        """
        if self.generators_status("structure"):
            self._check_if_defects()
            labels_list = self._build_label(lab_eff=lab_eff)
            if len(labels_list) > 0:
                particle, label_params_list = particle_from_structure(
                    self.structure, labels_list, self.configuration_path
                )
                if self.defect_eps["use_defects"]:
                    print("adding defects")
                    if defect_build is not None:
                        defect = defect_build
                    else:
                        defect = self.defect_eps["defect"]
                    particle.add_defects(
                        eps1=self.defect_eps["eps1"],
                        xmer_neigh_distance=self.defect_eps["eps2"],
                        deg_dissasembly=defect,
                    )
                else:
                    print("Particle without defects")
                if keep:
                    self.particle = particle
                    self.objects_created["particle"] = True
                return particle

    def _build_coordinate_field(
        self, use_self_particle=True, keep=False, coordinate_field_path=None, **kwargs
    ):
        """
        Build the coordinate field for the experiment, either from the current particle or as a minimal field.

        Parameters
        ----------
        use_self_particle : bool, optional
            If True, use the current particle to build the field. Default is True.
        keep : bool, optional
            If True, store the field in the experiment. Default is False.
        coordinate_field_path : str or None, optional
            Path to a coordinate field file. Default is None.
        **kwargs
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
            self.virtualsample_params["minimal_distance"] = fieldobject.molecules_params["minimal_distance"]
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
        use_local_field : bool, optional
            If True, use the local exported coordinate field. Default is False.
        prints : bool, optional
            If True, print status messages. Default is True.

        Returns
        -------
        None
        """
        if self.imaging_modalities:
            # print(f"Using selected mods: {self.imaging_modalities.keys()}")
            mods_list = list(self.imaging_modalities.keys())
            if use_local_field and self.generators_status("exported_coordinate_field"):
                self.imager, modality_parameters = create_imaging_system(
                    exported_field=self.exported_coordinate_field,
                    modalities_id_list=mods_list,
                    mod_params=self.imaging_modalities,
                    config_dir=self.configuration_path,
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
                )
            self.objects_created["imager"] = True
        else:
            self.imager = None
            print("No modalities")

    def generators_status(self, generator_name):
        """
        Return the status of the specified generator.

        Parameters
        ----------
        generator_name : str
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
        use_locals : bool, optional
            Determines whether to use local variables or attributes when building objects. Default is True.
        modules : list, optional
            List of module names to build. If "all" is included, builds all available modules: "structure", "particle", "coordinate_field", and "imager". Default is ["all"].
        **kwargs
            Additional keyword arguments for building modules.

        Notes
        -----
        The order of building is: structure, particle, coordinate_field, imager.
        """
        print("Building objects")
        if "all" in modules:
            build_list = ["structure", "particle", "coordinate_field", "imager"]
        else:
            build_list = modules
        if "structure" in build_list:
            self._build_structure()
        if "particle" in build_list:
            self._build_particle(keep=use_locals)
        if "coordinate_field" in build_list:
            self._build_coordinate_field(use_self_particle=use_locals, keep=use_locals)
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
        write : bool, optional
            Whether to write the reference image to disk. Default is False.
        keep : bool, optional
            If True, store the reference in the experiment. Default is False.
        ref_acq_pars : dict or None, optional
            Acquisition parameters for the reference. Default is None.
        modality_wise : bool, optional
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
                modalities_id_list=["Reference"], config_dir=self.configuration_path
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
                    ref_pixelsize=reference_imager.modalities["Reference"]["detector"][
                        "pixelsize"
                    ]
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
        name : str, optional
            Name of the experiment or simulation. Default is "vlab4mic_experiment".
        acq_params : dict or None, optional
            Acquisition parameters for the simulation. If None and modality is "All", uses self.selected_mods.
        save : bool, optional
            Whether to save the simulation output to disk. Default is False.
        modality : str, optional
            The imaging modality to simulate. If "All", simulates all available modalities. Default is "All".
        **kwargs
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
            simulation_output = generate_multi_imaging_modalities(
                image_generator=self.imager,
                experiment_name=name,
                savingdir=self.output_directory,
                write=save,
                acquisition_param=acq_params,
            )
            self.results = simulation_output
            return simulation_output
        else:
            print(f"Simulating: {modality}")
            acq_p = self.selected_mods[modality]
            timeseries, _ = self.imager.generate_imaging(modality=modality, **acq_p)
            simulation_output = {}
            simulation_output[modality] = timeseries
            return simulation_output

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
        probe_conjugation_efficiency: float = None,
        probe_seconday_epitope=None,
        probe_wobbling=False,
        labelling_efficiency: float = 1.0,
        as_primary=False,
        peptide_motif: dict = None,
        **kwargs,
    ):
        """
        Add a probe configuration to the experiment.

        This method loads a probe configuration from a YAML file, updates its parameters based on the provided arguments,
        and stores it in the experiment's probe parameters. It supports customization of probe properties such as target type,
        fluorophore, model, conjugation details, and more.

        Parameters
        ----------
        probe_name : str, optional
            Name of the probe. Default is "NHS_ester".
        probe_target_type : str, optional
            Type of the probe target (e.g., "Primary", "Secondary").
        probe_target_value : str, optional
            Value or identifier of the probe target.
        probe_target_option : str, optional
            Additional option for the probe target, used for secondary epitopes.
        probe_distance_to_epitope : float, optional
            Distance from the probe to the epitope.
        probe_model : str, optional
            List of model identifiers for the probe.
        probe_fluorophore : str, optional
            Identifier for the fluorophore. Default is "AF647".
        probe_steric_hindrance : Any, optional
            Steric hindrance value or configuration.
        probe_paratope : str, optional
            Paratope identifier or information.
        probe_conjugation_target_info : Any, optional
            Information about the conjugation target.
        probe_conjugation_efficiency : float, optional
            Efficiency of the probe conjugation.
        probe_seconday_epitope : Any, optional
            Information about a secondary epitope target.
        probe_wobbling : bool, optional
            Whether to enable probe wobbling. Default is False.
        labelling_efficiency : float, optional
            Efficiency of probe labelling. Default is 1.0.
        as_primary : bool, optional
            Whether to treat the probe as a primary linker. Default is False.
        **kwargs
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
        probe_configuration = copy.deepcopy(self.config_probe_params[probe_template])
        if probe_name is None:
            probe_name = probe_template
        else:
            probe_configuration["label_name"] = probe_name
        if peptide_motif is not None:
            protein_name, _1, site, sequence = self.structure.get_peptide_motif(**peptide_motif)
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
        elif probe_configuration["target"]["type"] is None or probe_configuration["target"]["value"] is None:
            print("No target info provided for the probe. Retrieving random sequence.")
            # probe has no target info
            # a random target will be used
            protein_name, _1, site, sequence = self.structure.get_peptide_motif(position="cterminal") 
            probe_configuration["target"] = dict(
                type="Sequence", value=sequence
            )
            #probe_configuration["target"]["type"] = "Sequence"
            #probe_configuration["target"]["value"] = sequence  
        if probe_distance_to_epitope is not None:
            probe_configuration["distance_to_epitope"] = probe_distance_to_epitope
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
        if probe_conjugation_efficiency is not None:
            probe_configuration["conjugation_efficiency"] = probe_conjugation_efficiency
        if probe_seconday_epitope is not None:
            probe_configuration["epitope_target_info"] = probe_seconday_epitope
        if probe_wobbling:
            probe_configuration["enable_wobble"] = probe_wobbling
        if as_primary:
            print("Adding probe as primary linker")
            probe_configuration["as_linker"] = True
        else:
            probe_configuration["as_linker"] = False
        if probe_steric_hindrance is not None:
            probe_configuration["distance_between_epitope"] = probe_steric_hindrance
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

    def set_virtualsample_params(
        self,
        virtualsample_template="square1x1um_randomised",
        sample_dimensions: list[int, int, int] = None,
        number_of_particles: int = None,
        particle_positions: list = None,
        random_orientations: bool = None,
        random_placing: bool = None,
        minimal_distance: float = None,
        **kwargs,
    ):
        """
        Set parameters for the virtual sample by loading a template configuration and optionally overriding specific values.

        Parameters
        ----------
        virtualsample_template : str, optional
            Name of the virtual sample template YAML file (without extension) to load. Default is "square1x1um_randomised".
        sample_dimensions : list of int, optional
            Dimensions of the sample to override the template value.
        number_of_particles : int, optional
            Number of particles to override the template value.
        particle_positions : list, optional
            List of particle positions to override the template value.
        random_orientations : bool, optional
            Whether to randomize particle orientations, overrides the template value.
        random_placing : bool, optional
            Whether to randomize particle placement, overrides the template value.
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        None
        """
        # load default configuration for virtual sample
        virtual_sample_template = os.path.join(
            self.configuration_path,
            "virtualsample",
            virtualsample_template + ".yaml",
        )
        vsample_configuration = load_yaml(virtual_sample_template)
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
        vsample_configuration["minimal_distance"] = minimal_distance
        self.virtualsample_params = vsample_configuration

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
        img : array-like
            The input image (as a NumPy array or compatible format) to be analyzed for feature detection.
        mode : str, optional
            The feature detection mode to use (e.g., "localmaxima"). Default is "localmaxima".
        sigma : float or None, optional
            Standard deviation for Gaussian smoothing. If None, no smoothing is applied.
        background : float or None, optional
            Background intensity to subtract from the image. If None, no subtraction.
        threshold : float or None, optional
            Minimum intensity threshold for feature detection. If None, uses default.
        pixelsize : float or None, optional
            Pixel size of the input image. If None, uses default.
        min_distance : float or None, optional
            Minimum distance between detected features. If None, uses default.
        **kwargs
            Additional keyword arguments passed to the feature detection function.

        Returns
        -------
        None

        Notes
        -----
        Updates `self.virtualsample_params` with new relative positions and sample dimensions. Calls `self.build()` for the specified modules.
        """
        xyz_relative, image_physical_size = coordinates_field.gen_positions_from_image(
            img=img,
            mode=mode,
            sigma=sigma,
            background=background,
            threshold=threshold,
            pixelsize=pixelsize,
            min_distance=min_distance,
            **kwargs,
        )
        self.set_virtualsample_params(
            sample_dimensions=[
                image_physical_size[0],
                image_physical_size[1],
                100,],
            particle_positions=xyz_relative,
            number_of_particles=len(xyz_relative),
            random_orientations=False,
            random_placing=False,
            minimal_distance=min_distance,
        )
        #self.virtualsample_params["relative_positions"] = xyz_relative
        #self.virtualsample_params["sample_dimensions"] = [
        #    image_physical_size[0],
        #    image_physical_size[1],
        #    100,
        #]
        self.build(modules=["coordinate_field", "imager"])

    def current_settings(self, as_string=True, newline="<br>", modalities_acq_params=False):
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
        


def generate_virtual_sample(
    structure: str = "1XI5",
    probe_name: str = None,
    probe_target_type: str = None,
    probe_target_value: str = None,
    probe_distance_to_epitope: float = None,
    probe_model: list[str] = None,
    probe_fluorophore: str = "AF647",
    probe_paratope: str = None,
    probe_conjugation_target_info=None,
    probe_conjugation_efficiency: float = None,
    probe_seconday_epitope=None,
    probe_wobbling=False,
    labelling_efficiency: float = 1.0,
    defect_small_cluster: float = None,
    defect_large_cluster: float = None,
    defect: float = None,
    virtual_sample_template: str = "square1x1um_randomised",
    sample_dimensions: list[float] = None,
    number_of_particles: int = None,
    particle_positions: list[np.array] = None,
    random_orientations=False,
    random_placing=False,
    clear_probes=False,
    **kwargs,
):
    """
    Create a virtual sample from a structure model.

    Parameters
    ----------
    structure : str, optional
        4-letter ID of PDB/CIF model. Default is "1XI5".
    probe_name : str, optional
        Name ID of probe configuration file (filename).
    probe_target_type : str, optional
        Options: "Sequence", "Atom_residue", or "Primary".
    probe_target_value : str or dict, optional
        For target type "Sequence" or "Primary", a string. For "Atom_residue", a dictionary with keys "Atom" and "Residue".
    probe_distance_to_epitope : float, optional
        Minimal distance set from epitope and probe paratope.
    probe_model : list of str, optional
        4-letter ID(s) of PDB/CIF model(s).
    probe_fluorophore : str, optional
        Fluorophore name (e.g., "AF647"). Default is "AF647".
    probe_paratope : str, optional
        Sequence of the paratope site for when probe includes a model.
    probe_conjugation_target_info : Any, optional
        Information about the probe conjugation target.
    probe_conjugation_efficiency : float, optional
        Efficiency of conjugation of emitters.
    probe_seconday_epitope : str, optional
        Sequence within probe model to be used as epitope for a secondary.
    probe_wobbling : bool, optional
        Enable probe wobbling. Default is False.
    labelling_efficiency : float, optional
        Labelling efficiency of probe. Default is 1.0.
    defect_small_cluster : float, optional
        In Å, distance used to group epitopes into multimers.
    defect_large_cluster : float, optional
        In Å, distance within multimers to consider neighbors.
    defect : float, optional
        Fraction of defect to model.
    virtual_sample_template : str, optional
        Name of the configuration file for template. Default is "square1x1um_randomised".
    sample_dimensions : list of float, optional
        In nanometers, define the X, Y, and Z sizes of the field.
    number_of_particles : int, optional
        Number of independent copies of a particle to create and distribute.
    particle_positions : list of np.array, optional
        Relative positions of particles in the field.
    random_orientations : bool, optional
        If True, each particle will be randomly assigned a new orientation. Default is False.
    random_placing : bool, optional
        Define if position in field is random or the center of field. Default is False.
    clear_probes : bool, optional
        If True, default parameters will be cleared. Default is False.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    tuple
        - dict: The virtual sample as exported format. This can be used as input for the `image_vsample` method.
        - ExperimentParametrisation: The experiment containing all modules that were generated to build the virtual sample, and the virtual sample module itself. This experiment can be further used and tweaked for subsequent analysis or branching workflows.
    """
    myexperiment = ExperimentParametrisation()
    # load default configuration for probe
    if not clear_probes:
        if probe_name is None:
            probe_name = "NHS_ester"
        probe_configuration_file = os.path.join(
            myexperiment.configuration_path, "probes", probe_name + ".yaml"
        )
        probe_configuration = load_yaml(probe_configuration_file)
        probe_configuration["probe_name"] = probe_name
        if probe_target_type and probe_target_value:
            probe_configuration["target_info"] = dict(
                type=probe_target_type, value=probe_target_value
            )
        if probe_distance_to_epitope is not None:
            probe_configuration["distance_to_epitope"] = probe_distance_to_epitope
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
        if probe_conjugation_efficiency is not None:
            probe_configuration["conjugation_efficiency"] = probe_conjugation_efficiency
        if probe_seconday_epitope is not None:
            probe_configuration["epitope_target_info"] = probe_seconday_epitope
        if probe_wobbling:
            probe_configuration["enable_wobble"] = probe_wobbling
        myexperiment.add_probe(**probe_configuration)
    # load default configuration for virtual sample
    virtual_sample_template = os.path.join(
        myexperiment.configuration_path,
        "virtualsample",
        virtual_sample_template + ".yaml",
    )
    vsample_configuration = load_yaml(virtual_sample_template)
    myexperiment.configuration_path
    myexperiment.structure_id = structure
    # myexperiment.structure_label = probe_name
    # myexperiment.probe_parameters[probe_name] = probe_configuration

    if defect and defect_large_cluster and defect_small_cluster:
        myexperiment.defect_eps["eps1"] = defect_small_cluster
        myexperiment.defect_eps["eps2"] = defect_large_cluster
        myexperiment.defect_eps["defect"] = defect
        myexperiment.defect_eps["use_defects"] = True

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
    myexperiment.virtualsample_params = vsample_configuration
    myexperiment.build(use_locals=True)

    # myexperiment.coordinate_field_id = virtual_sample
    return myexperiment.exported_coordinate_field, myexperiment


def build_virtual_microscope(
    modality="STED", multimodal: list[str] = None, experiment=None, **kwargs
):
    """
    Initialise a virtual microscope for single or multimodal imaging.

    Parameters
    ----------
    modality : str, optional
        Modality name. Default is "STED".
    multimodal : list of str, optional
        List of modality names. Overrides the `modality` parameter if provided.
    experiment : ExperimentParametrisation, optional
        An Experiment object. If None, a new one is created.
    **kwargs
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
    modality="STED",
    multimodal: list[str] = None,
    run_simulation=True,
    **kwargs,
):
    """
    Generate imaging simulations of the specified virtual sample and imaging modalities.

    If a virtual sample is provided, a virtual microscope is created around it. In this case, the resulting experiment object will only contain the imager and the virtual sample loaded into it.

    If a virtual sample is not provided, a default sample will be created along with any keyword provided that specifies structure, probes, or virtual sample parameters. In this case, the resulting experiment will also contain initialised modules for structure, probes, particle, and coordinates field (generator of the virtual sample).

    Parameters
    ----------
    vsample : dict, optional
        Dictionary specifying sample parameters. Corresponds to Experiment attribute "exported_coordinated_field".
    modality : str, optional
        Modality name. Default is "STED".
    multimodal : list of str, optional
        List of modality names. Overrides the `modality` parameter if provided.
    run_simulation : bool, optional
        If True, generates image simulation for each modality set. Default is True.
    **kwargs
        Additional arguments passed to `add_modality`.

    Returns
    -------
    tuple
        - dict: Image simulations.
        - ExperimentParametrisation: The experiment object as described above.
    """
    if vsample is None:
        vsample, sample_experiment = generate_virtual_sample(**kwargs)
        if multimodal is not None:
            for mod in multimodal:
                sample_experiment.add_modality(mod, **kwargs)
        else:
            sample_experiment.add_modality(modality, **kwargs)
        sample_experiment.build(
            modules=[
                "imager",
            ]
        )
    else:
        # need to create an experiment for it
        if multimodal is not None:
            vmicroscope, sample_experiment = build_virtual_microscope(
                multimodal=multimodal, use_local_field=True, **kwargs
            )
        else:
            vmicroscope, sample_experiment = build_virtual_microscope(
                modality=modality, use_local_field=True, **kwargs
            )
        sample_experiment.imager.import_field(**vsample)
    imaging_output = dict()
    if run_simulation:
        imaging_output = sample_experiment.run_simulation()
    return imaging_output, sample_experiment
