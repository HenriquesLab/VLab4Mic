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
    output_directory: str = ""

    def __post_init__(self):
        pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
        local_dir = os.path.join(pck_dir, "configuration")
        self.configuration_path = local_dir
        # keep track of objects created
        self.objects_created = dict(
            structure=False,
            particle=False,
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
        labels_list = self._build_label(lab_eff=lab_eff)
        particle = particle_from_structure(
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

    def _build_coordinate_field(self, use_self_particle=True, keep=False):
        if use_self_particle and self.particle:
            exported_field, fieldobject = field_from_particle(
                self.particle, field_config=self.coordinate_field_id
            )
        else:
            pass
        if keep:
            self.exported_coordinate_field = exported_field
        return exported_field

    def _build_imager(self):
        if self.selected_mods:
            mods_list = list(self.selected_mods.keys())
            self.imager = create_imaging_system(
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

    def build(self):
        self._build_imager()
        self._build_structure()
        self._param_linspaces()

    def gen_reference(self, write=False, keep=False):
        reference_pars = dict()
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
        self.imager.import_field(**tmp_exported_field)
        #
        output_name = "REFERENCE_"
        _reference = generate_multi_imaging_modalities(
            image_generator=self.imager,
            experiment_name=output_name,
            savingdir=self.output_directory,
            write=write,
        )
        if keep:
            self.experiment_reference = _reference
            self.objects_created["output_reference"] = True
        return _reference


def create_experiment_parametrisation(
    structure_and_labels: dict,
    modalities_acquisition: dict,
    defects_params: dict,
    params2sweep: dict,
    savging: dict,
):
    generator = ExperimentParametrisation()
    # Structural parameters
    generator.structure_id = structure_and_labels["structure_id"]
    generator.structure_label = structure_and_labels["structure_label"]
    generator.fluorophore_id = structure_and_labels["fluorophore_id"]
    for mod, acquisition in modalities_acquisition.items():
        generator.selected_mods[mod] = acquisition
    for key, val in defects_params.items():
        generator.defect_eps[key] = val
    for parametername, pars in params2sweep.items():
        generator.sweep_pars[parametername] = pars
    # writing
    generator.experiment_id = savging["experiment_id"]
    generator.output_directory = savging["output_directory"]
    generator.build()
    return generator
