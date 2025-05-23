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

from pathlib import Path

output_path = Path.home() / "vlab4mic_outputs"

if not os.path.exists(output_path):
    os.makedirs(output_path)

@dataclass
class ExperimentParametrisation:
    experiment_id: str = ""
    structure_id: str = ""
    configuration_path: str = ""
    structure_label: str = ""
    fluorophore_id: str = ""
    coordinate_field_id: str = None
    selected_mods: Dict[str, int] = field(default_factory=dict)
    imaging_modalities: Dict[str, int] = field(default_factory=dict)
    probe_parameters: Dict[str, int] = field(default_factory=dict)
    defect_eps: Dict[str, int] = field(default_factory=dict)
    sweep_pars: Dict[str, int] = field(default_factory=dict)
    objects_created: Dict[str, int] = field(default_factory=dict)
    output_directory: str = str(output_path)


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
            mod_info = configuration_format.compile_modality_parameters(
                mod, local_dir
            )
            modality_parameters[mod] = mod_info
        self.local_modalities_names = modalities_names_list
        self.local_modalities_parameters = modality_parameters
        #self.imaging_modalities = dict()

    def add_modality(self, modality_name, save=False, **kwargs):
        if modality_name in self.local_modalities_names:
            self.imaging_modalities[modality_name] = copy.deepcopy(self.local_modalities_parameters[modality_name])
            for param, value in kwargs.items():
                self.imaging_modalities[modality_name][param] = value
            self.set_modality_acq(modality_name, save=save)
        else:
            print(f"Modality {modality_name} not found in demo modalities. Using Widefield params as template to create new.")
            self.imaging_modalities[modality_name] = copy.deepcopy(self.local_modalities_parameters["Widefield"])
            for param, value in kwargs.items():
                self.imaging_modalities[modality_name][param] = value
            self.set_modality_acq(modality_name, save=save)

    def update_modality(self,
                        modality_name,
                        pixelsize_nm: int = None,
                        lateral_resolution_nm: int = None,
                        axial_resolution_nm: int = None,
                        psf_voxel_nm: int = None,
                        remove = False,
                        ):
        if remove:
            self.imaging_modalities.pop(modality_name, None)
        else:
            changes = False
            if pixelsize_nm is not None:
                self.imaging_modalities[modality_name]["detector"]["pixelsize"] = pixelsize_nm / 1000
                changes = True
            if lateral_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name]["psf_params"]["voxelsize"][0]
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][0] = lateral_resolution_nm/voxel_size
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][1] = lateral_resolution_nm/voxel_size
                changes = True
            if axial_resolution_nm is not None:
                voxel_size = self.imaging_modalities[modality_name]["psf_params"]["voxelsize"][0]
                self.imaging_modalities[modality_name]["psf_params"]["std_devs"][2] = axial_resolution_nm/voxel_size
                changes = True
            if psf_voxel_nm is not None:
                self.imaging_modalities[modality_name]["psf_params"]["voxelsize"] = [psf_voxel_nm, psf_voxel_nm, psf_voxel_nm]
                changes = True
            if changes:
                #mod_pars = self.imaging_modalities[modality_name]
                self.imager.set_imaging_modality(**self.imaging_modalities[modality_name])

    def set_modality_acq(
            self,
            modality_name, 
            exp_time=0.001,
            noise=True,
            save=False,
            nframes=2,
            channels=["ch0",],
            **kwargs
        ):
        if modality_name in self.imaging_modalities.keys():
            self.selected_mods[modality_name] = configuration_format.format_modality_acquisition_params(
               exp_time,
               noise,
               save,
               nframes,
               channels
            )
        else:
            print("Modality not selected")
    
    def reset_to_defaults(self, module = "acquisitions", **kwargs):
        if module == "acquisitions":
            for mod_name in self.selected_mods.keys():
                self.set_modality_acq(mod_name, **kwargs)

    def _build_structure(self, keep=True):
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
        Create a list with one dictionary with
        the set up info for creating a label.
        It assumes there is only one label store in
        the attrubute "structure_label" and one in
        "fluorophore_id".
        """
        labels_list = []
        for probe_name, probe_params in self.probe_parameters.items():
            labels_list.append(
                label_builder_format(
                    label_id=probe_name,
                    **probe_params
                )
            )
        return labels_list
        #probe_name = self.structure_label
        #if probe_name in self.probe_parameters.keys():
        #    probeparams = self.probe_parameters[probe_name]
        #    labels_list.append(
        #        label_builder_format(
        #            label_id=probe_name,
        #            **probeparams
        #        )
        #    )
        #else:
        #    labels_list.append(
        #        label_builder_format(
        #            label_id=self.structure_label,
        #            fluorophore_id=self.fluorophore_id,
        #            labelling_efficiency=lab_eff,
        #        )
        #    )
        #if keep:
        #    pass
        #else:
        #    return labels_list

    def _check_if_defects(self):
        if self.defect_eps["defect"] and self.defect_eps["eps1"] and self.defect_eps["eps2"]:
            self.defect_eps["use_defects"] = True

    def _build_particle(self, lab_eff=1.0, defect_build=None, keep=False):
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
        self.exported_coordinate_field = None
        if use_self_particle and self.generators_status("particle"):
            print("creating field from existing particle")
            exported_field, fieldobject = field_from_particle(
                self.particle, **self.virtualsample_params, **kwargs
            )
            if keep:
                self.exported_coordinate_field = exported_field
                self.objects_created["exported_coordinate_field"] = True
                self.coordinate_field = fieldobject
                self.objects_created["coordinate_field"] = True
            return exported_field
        else:
            # create minimal field
            fieldobject = coordinates_field.create_min_field(**kwargs)
            exported_field = fieldobject.export_field()
            if keep:
                self.exported_coordinate_field = exported_field
                self.objects_created["exported_coordinate_field"] = True
                self.coordinate_field = fieldobject
                self.objects_created["coordinate_field"] = True

    def _build_imager(self, use_local_field=False, prints=True):
        if self.imaging_modalities:
            #print(f"Using selected mods: {self.imaging_modalities.keys()}")
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
                    print("Local field missing or unused. Creating imager without particles")
                self.imager, modality_parameters = create_imaging_system(
                    modalities_id_list=mods_list,
                    mod_params=self.imaging_modalities,
                    config_dir=self.configuration_path
                )
            self.objects_created["imager"] = True
        else:
            print("No modalities")

    def generators_status(self, generator_name):
        return self.objects_created[generator_name]

    def build(self, use_locals=True, modules: list = ["all",], **kwargs):
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
        # self._param_linspaces()

    def _gen_reference(self, write=False, keep=False, ref_acq_pars=None, modality_wise=False):
        """
        Calculate a reference image of the virtual sample by using the ideal
        parameters for each of the parameters to sweep. Requires the 
        dictionary of the params2sweep


        If modality_wise is true (default), a reference is calculated for each modality
        with that modality's parameters. Otherwhise, it will use a the Reference modality 
        configuration to generate an idealised image as reference.
        """
        reference_pars = dict()
        output_name = "REFERENCE_"
        # get ideal parameters for virtual sample
        for param_name, param_settings in self.sweep_pars.items():
            # print(param_name, param_settings)
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
                config_dir=self.configuration_path
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
                scalefactor = np.ceil(imager_scale / 1e-9)  # resulting pixel size in nanometers
                _reference_parameters[mod_name] = dict(
                    ref_pixelsize=reference_imager.modalities["Reference"]["detector"]["pixelsize"]*scalefactor 
                )
        #
        if keep:
            self.experiment_reference = _reference
            self.objects_created["output_reference"] = True
        return _reference, _reference_parameters

    def run_simulation(self, name="NONAME", acq_params=None, save=False, modality="All", **kwargs):
        # imager will run regardless, since by default
        # has a minimal coordinate field
        if modality == "All":
            print("Simulating all modalities")
            if acq_params is None:
                acq_params=self.selected_mods
            if self.experiment_id:
                name = self.experiment_id
            if self.generators_status("imager"):
                simulation_output = generate_multi_imaging_modalities(
                    image_generator=self.imager,
                    experiment_name=name,
                    savingdir=self.output_directory,
                    write=save,
                    # acq_params is a value in selected mods
                    acquisition_param=acq_params,
                )
                return simulation_output
            else:
                print("Missing attributes")
        else:
            print(f"Simulating: {modality}")
            acq_p = self.selected_mods[modality]
            timeseries, _ = self.imager.generate_imaging(
                modality=modality, **acq_p
            )
            simulation_output = {}
            simulation_output[modality] = timeseries
            return simulation_output

    def remove_probes(self):
        self.probe_parameters = dict()
        self._update_probes()
        if self.generators_status("particle"):
            self.particle = None
        if self.generators_status("structure"):
            self.structure._clear_labels()
            #self.structure.label_targets = dict()

    def add_probe(
            self,
            probe_name: str = "NHS_ester",
            probe_target_type: str = None,
            probe_target_value: str = None,
            probe_target_option: str = None,
            probe_distance_to_epitope: float = None,
            probe_model: list[str] = None,
            probe_fluorophore: str = "AF647",
            probe_steric_hindrance = None,
            probe_paratope: str = None,
            probe_conjugation_target_info = None,
            probe_conjugation_efficiency: float = None,
            probe_seconday_epitope = None,
            probe_wobbling = False,
            labelling_efficiency: float = 1.0,
            as_primary=False,
            **kwargs):
        probe_configuration_file = os.path.join(self.configuration_path, "probes", probe_name + ".yaml")
        probe_configuration = load_yaml(probe_configuration_file)
        if probe_target_type and probe_target_value:
            probe_configuration["target_info"] = dict(type=probe_target_type, value=probe_target_value)
            if probe_target_type == "Primary" and probe_target_option:
                # check if there is a primary probe with the name of value
                if probe_target_value in self.probe_parameters.keys():
                    print(f"Using {probe_target_option} as epitope on {probe_target_value}")
                    self.probe_parameters[probe_target_value]["probe_seconday_epitope"] = probe_target_option
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
            probe_configuration["conjugation_target_info"] = probe_conjugation_target_info
        if probe_conjugation_efficiency is not None:
            probe_configuration["conjugation_efficiency"] = probe_conjugation_efficiency
        if probe_seconday_epitope is not None:
            probe_configuration["epitope_target_info"] = probe_seconday_epitope
        if probe_wobbling:
            probe_configuration["enable_wobble"] = probe_wobbling
        if as_primary:
            probe_configuration["as_linker"] = as_primary
        if probe_steric_hindrance is not None:
            probe_configuration["distance_between_epitope"] = probe_steric_hindrance
        self.probe_parameters[probe_name] = probe_configuration
        self._update_probes()

    def _update_probes(self):
        if len(self.probe_parameters.keys()) == 0:
            self.structure_label = None
        else: 
            self.structure_label = list(self.probe_parameters.keys())
    
    def use_image_for_positioning(self,
                                  img,
                                  mode="localmaxima",
                                  sigma=None,
                                  background=None,
                                  threshold=None,
                                  pixelsize=None,
                                  min_distance=None,
                                  **kwargs):
        xyz_relative, image_physical_size = coordinates_field.gen_positions_from_image(
            img=img,
            mode=mode,
            sigma=sigma,
            background=background,
            threshold=threshold,
            pixelsize=pixelsize,
            min_distance=min_distance
            )
        self.virtualsample_params["relative_positions"] = xyz_relative
        self.virtualsample_params["sample_dimensions"] = [image_physical_size[0],image_physical_size[1], 100]
        self.build(modules=["coordinate_field", "imager"])

def generate_virtual_sample(
        structure: str = "1XI5",
        probe_name: str = None,
        probe_target_type: str = None,
        probe_target_value: str = None,
        probe_distance_to_epitope: float = None,
        probe_model: list[str] = None,
        probe_fluorophore: str = "AF647",
        probe_paratope: str = None,
        probe_conjugation_target_info = None,
        probe_conjugation_efficiency: float = None,
        probe_seconday_epitope = None,
        probe_wobbling = False,
        labelling_efficiency: float = 1.0,
        defect_small_cluster: float = None,
        defect_large_cluster: float = None,
        defect: float = None,
        virtual_sample_template: str = "square1x1um_randomised",
        sample_dimensions: list[float] = None,
        number_of_particles: int = None,
        particle_positions: list[np.array] = None,
        random_orientations = False,
        random_placing = False,
        clear_probes=False,
        **kwargs
):
    myexperiment = ExperimentParametrisation()
    # load default configuration for probe
    if not clear_probes:
        if probe_name is None:
            probe_name = "NHS_ester"
        probe_configuration_file = os.path.join(myexperiment.configuration_path, "probes", probe_name + ".yaml")
        probe_configuration = load_yaml(probe_configuration_file)
        probe_configuration["probe_name"] = probe_name
        if probe_target_type and probe_target_value:
            probe_configuration["target_info"] = dict(type=probe_target_type, value=probe_target_value)
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
            probe_configuration["conjugation_target_info"] = probe_conjugation_target_info
        if probe_conjugation_efficiency is not None:
            probe_configuration["conjugation_efficiency"] = probe_conjugation_efficiency
        if probe_seconday_epitope is not None:
            probe_configuration["epitope_target_info"] = probe_seconday_epitope
        if probe_wobbling:
            probe_configuration["enable_wobble"] = probe_wobbling
        myexperiment.add_probe(**probe_configuration)
    # load default configuration for virtual sample
    virtual_sample_template = os.path.join(myexperiment.configuration_path, "virtualsample", virtual_sample_template + ".yaml")
    vsample_configuration = load_yaml(virtual_sample_template)
    myexperiment.configuration_path
    myexperiment.structure_id = structure
    #myexperiment.structure_label = probe_name
    #myexperiment.probe_parameters[probe_name] = probe_configuration
    

    if defect and defect_large_cluster and defect_small_cluster:
        myexperiment.defect_eps["eps1"] = defect_small_cluster
        myexperiment.defect_eps["eps2"] = defect_large_cluster
        myexperiment.defect_eps["defect"] = defect
        myexperiment.defect_eps["use_defects"] = True

    if sample_dimensions is not None:
        vsample_configuration["sample_dimensions"] = sample_dimensions
    if number_of_particles is not None:
        vsample_configuration["nparticles"] = number_of_particles
    if particle_positions is not None:
        vsample_configuration["relative_positions"] = particle_positions
    if random_orientations is not None:
        vsample_configuration["random_orientations"] = random_orientations
    if random_placing is not None:
        vsample_configuration["random_placing"] = random_placing
    myexperiment.virtualsample_params = vsample_configuration
    myexperiment.build(use_locals=True)

    #myexperiment.coordinate_field_id = virtual_sample
    return myexperiment.exported_coordinate_field, myexperiment


def build_virtual_microscope(
        modality = "STED",
        multimodal: list[str] = None,
        experiment = None,
        **kwargs
):
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
        vsample = None, 
        modality = "STED",
        multimodal: list[str] = None,
        run_simulation = True,
        **kwargs
        ):
    if vsample is None:
        vsample, sample_experiment = generate_virtual_sample(**kwargs)
        if multimodal is not None:
            for mod in multimodal:
                sample_experiment.add_modality(mod, **kwargs)
        else:
            sample_experiment.add_modality(modality, **kwargs)
        sample_experiment.build(modules=["imager",])
    else:
        # need to create an experiment for it
        if multimodal is not None:
            vmicroscope, sample_experiment = build_virtual_microscope(
                multimodal=multimodal,
                use_local_field=True,
                **kwargs)
        else:
            vmicroscope, sample_experiment = build_virtual_microscope(
                modality=modality,
                use_local_field=True,
                **kwargs)
        sample_experiment.imager.import_field(**vsample)
    imaging_output = dict()
    if run_simulation:
        imaging_output = sample_experiment.run_simulation()
    return imaging_output, sample_experiment