from dataclasses import dataclass, field, fields
from typing import List, Dict
import numpy as np
import supramolsim
from supramolsim.workflows import *
import os


@dataclass
class ExperimentParametrisation:
    experiment_id: str = ""
    structure_id: str = ""
    configuration_path: str = ""
    structure_label: str = ""
    fluorophore_id: str = ""
    selected_mods: Dict[str, int] = field(default_factory=dict)
    defect_eps: Dict[str, int] = field(default_factory=dict)
    sweep_pars: Dict[str, int] = field(default_factory=dict)
    output_directory: str = ""

    def __post_init__(self):
        pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
        local_dir = os.path.join(pck_dir, "configuration")
        self.configuration_path = local_dir

    def _build_structure(self):
        if self.structure_id:
            struct, struct_param = load_structure(
                self.structure_id, self.configuration_path
            )
            self.structure = struct
        else:
            print("No structure ID")

    def _build_imager(self):
        if self.selected_mods:
            mods_list = list(selected_mods.keys())
            self.imager = create_imaging_system(
                modalities_id_list=mods_list, config_dir=self.configuration_path
            )
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
    return generator
