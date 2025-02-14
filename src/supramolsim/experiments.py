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
from .utils.data_format.structural_format import label_builder_format
from .utils.data_format import configuration_format
import os


@dataclass
class ExperimentParametrisation:
    experiment_id: str = ""
    structure_id: str = ""
    configuration_path: str = ""
    structure_label: str = ""
    fluorophore_id: str = ""
    coordinate_field_id: str = None
    selected_mods: Dict[str, int] = field(default_factory=dict)
    defect_eps: Dict[str, int] = field(default_factory=dict)
    sweep_pars: Dict[str, int] = field(default_factory=dict)
    objects_created: Dict[str, int] = field(default_factory=dict)
    output_directory: str = None

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
        labels_list.append(
            label_builder_format(
                label_id=self.structure_label,
                fluorophore_id=self.fluorophore_id,
                labelling_efficiency=lab_eff,
            )
        )
        if keep:
            pass
        else:
            return labels_list

    def _build_particle(self, lab_eff=1.0, defect=0.0, keep=False):
        if self.generators_status("structure"):
            labels_list = self._build_label(lab_eff=lab_eff)
            particle, label_params_list = particle_from_structure(
                self.structure, labels_list, self.configuration_path
            )
            if self.defect_eps:
                particle.add_defects(
                    eps1=self.defect_eps["eps1"],
                    xmer_neigh_distance=self.defect_eps["eps2"],
                    deg_dissasembly=defect,
                )
            if keep:
                self.particle = particle
                self.objects_created["particle"] = True
            return particle

    def _build_coordinate_field(
        self, use_self_particle=True, keep=False, coordinate_field_path=None, **kwargs
    ):
        if use_self_particle and self.generators_status("particle"):
            print("creating field from existing particle")
            exported_field, fieldobject = field_from_particle(
                self.particle, field_config=coordinate_field_path, **kwargs
            )
            if keep:
                self.exported_coordinate_field = exported_field
                self.objects_created["exported_coordinate_field"] = True
                self.coordinate_field = fieldobject
                self.objects_created["coordinate_field"] = True
            return exported_field
        else:
            # create minimal field
            print("else")
            pass

    def _build_imager(self, use_local_field=False):
        if self.selected_mods:
            print(f"Using selected mods: {self.selected_mods}")
            mods_list = list(self.selected_mods.keys())
            if use_local_field and self.generators_status("exported_coordinate_field"):
                self.imager, modality_parameters = create_imaging_system(
                    exported_field=self.exported_coordinate_field,
                    modalities_id_list=mods_list,
                    config_dir=self.configuration_path,
                )
            else:
                print("Local field missing or unused. Creating imager without particles")
                self.imager, modality_parameters = create_imaging_system(
                    modalities_id_list=mods_list, config_dir=self.configuration_path
                )
            self.objects_created["imager"] = True
        else:
            print("No modalities")

    def _param_linspaces(self):
        if self.sweep_pars:
            self.sweep_linspaces = dict()
            for param_name, pars in self.sweep_pars.items():
                self.sweep_linspaces[param_name] = np.linspace(
                    pars["start"], pars["end"], pars["nintervals"]
                )
        else:
            print("No parameters set")

    def generators_status(self, generator_name):
        return self.objects_created[generator_name]

    def build(self, use_locals=False):
        print("Building objects")
        self._build_structure()
        self._build_particle(keep=use_locals)
        self._build_coordinate_field(use_self_particle=use_locals, keep=use_locals)
        self._build_imager(use_local_field=use_locals)
        self._param_linspaces()

    def gen_reference(self, write=False, keep=False, ref_acq_pars=None, modality_wise=False):
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

    def run_simulation(self, name="NONAME", acq_params=None, save=False):
        # imager will run regardless, since by default
        # has a minimal coordinate field
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


def create_experiment_parametrisation(
    structure_and_labels: dict = None,
    modalities_acquisition: dict = None,
    field_params: dict = None,
    savging: dict = None,
    defects_params: dict = None,
    params2sweep: dict = None,
    use_locals=False,
):
    generator = ExperimentParametrisation()
    # Structural parameters
    if structure_and_labels:
        generator.structure_id = structure_and_labels["structure_id"]
        generator.structure_label = structure_and_labels["structure_label"]
        generator.fluorophore_id = structure_and_labels["fluorophore_id"]
    if modalities_acquisition:
        for mod, acquisition in modalities_acquisition.items():
            if acquisition is None:
                generator.selected_mods[mod] = configuration_format.format_modality_acquisition_params()
            else:
                generator.selected_mods[mod] = acquisition
    if field_params:
        #for field_param, value in field_params.items():
        generator.coordinate_field_id=field_params
    if defects_params:
        for key, val in defects_params.items():
            generator.defect_eps[key] = val
    if params2sweep:
        for parametername, pars in params2sweep.items():
            generator.sweep_pars[parametername] = pars
    # writing
    if savging:
        generator.experiment_id = savging["experiment_id"]
        generator.output_directory = savging["output_directory"]
    generator.build(use_locals=use_locals)
    return generator
