from .experiments import ExperimentParametrisation
from .analysis import sweep
import matplotlib.pyplot as plt
from IPython.utils import io
import os
from pathlib import Path

output_dir = Path.home() / "vlab4mic_outputs"

class sweep_generator:
    #sweep parameters
    experiment = ExperimentParametrisation()
    structures = [
        "1XI5",
    ]
    probes = [
        "NHS_ester",
    ]
    modalities = ["Widefield", "Confocal", "STED", "SMLM"]
    sweep_repetitions = 3
    probe_parameters = None
    defect_parameters = None
    vsample_parameters = None
    modality_parameters = None
    acquisition_parameters = None
    # parameter dictionary
    params_by_group = {}
    # reference parameters
    reference_structure = "1XI5"
    reference_probe = "NHS_ester"
    reference_probe_parameters = {"labelling_efficiency": 1.0}
    ####  outputs
    ouput_directory = output_dir
    reference_virtual_sample = None
    reference_virtual_sample_params = None
    reference_image = None
    virtual_samples = None
    virtual_samples_parameters = None
    acquisition_outputs = None
    acquisition_outputs_parameters = None
    # analysis
    analysis = {}
    analysis["measurements"] = None
    analysis["inputs"] = None
    analysis["extended_dataframe"] = None

    def __init__(self):
        self.params_by_group["probe"] = {}
        self.params_by_group["virtual_sample"] = {}
        self.params_by_group["particle_defect"] = {}
        self.params_by_group["modality"] = {}
        self.params_by_group["acquisition"] = {}

    # generators
    def generate_virtual_samples(self):
        self.create_parameters_iterables()
        self.experiment, self.virtual_samples, self.virtual_samples_parameters = (
            sweep.sweep_vasmples(
                structures=self.structures,
                probes=self.probes,
                probe_parameters=self.probe_parameters,
                virtual_samples=self.vsample_parameters,
                repetitions=self.sweep_repetitions,
            )
        )

    def generate_acquisitions(self):
        # generate virtual samples
        with io.capture_output() as captured:
            if self.virtual_samples is None:
                self.generate_virtual_samples()
        # acquisition of virtual samples
        (
            self.experiment,
            self.acquisition_outputs,
            self.acquisition_outputs_parameters,
            mod_pixelsizes,
        ) = sweep.sweep_modalities_updatemod(
            experiment=self.experiment,
            vsample_outputs=self.virtual_samples,
            vsampl_pars=self.virtual_samples_parameters,
            modalities=self.modalities,
            modality_acq_prams=self.acquisition_parameters,
        )

    def set_reference_parameters(
            self,
            reference_structure: str = None,
            reference_probe: str = None,
            reference_probe_parameters: dict = None,
            **kwargs):
        if reference_structure is not None:
            self.reference_structure = reference_structure
        if reference_probe is not None:
            self.reference_probe = reference_probe
        if reference_probe_parameters is not None:
            self.reference_probe_parameters = reference_probe_parameters
        

    def generate_reference_sample(self):
        self.reference_virtual_sample, self.reference_virtual_sample_params = (
            sweep.generate_global_reference_sample(
                structure=self.reference_structure,
                probe=self.reference_probe,
                probe_parameters=self.reference_probe_parameters,
            )
        )

    def generate_reference_image(self):
        if self.reference_image is None:
            self.generate_reference_sample()
        self.reference_image, self.reference_image_parameters = (
            sweep.generate_global_reference_modality(
                reference_vsample=self.reference_virtual_sample,
                reference_vsample_params=self.reference_virtual_sample_params,
            )
        )

    # previews
    def preview_acquisition_output(self, return_image=False):
        param_id = list(self.acquisition_outputs_parameters.keys())[0]
        repetition = 0
        frame = 0
        if return_image:
            return self.acquisition_outputs[param_id][repetition][frame]
        else:
            plt.imshow(self.acquisition_outputs[param_id][repetition][frame])
            print(self.acquisition_outputs_parameters[param_id])

    def preview_reference_image(self, return_image=False):
        if return_image:
            return self.reference_image
        else:
            plt.imshow(self.reference_image[0])
            print(self.reference_image_parameters)

    # set and change parameters
    def set_param_range(self, param_group, param_name, param_type, first=None, last=None, option=None):
        if param_type == "numeric":
            self.params_by_group[param_group][param_name] = [first, last, option]
        if param_type == "logical":
            logic_list = []
            if option == "True":
                logic_list = [True,]
            elif option == "True":
                logic_list = [False,]
            elif option == "Both":
                logic_list = [True, False,]
            self.params_by_group[param_group][param_name] = logic_list

    def create_parameters_iterables(self):
        if self.params_by_group["probe"]:
            self.probe_parameters = sweep.probe_parameters_sweep(
                **self.params_by_group["probe"]
            )
        if self.params_by_group["particle_defect"]:
            pass
        if self.params_by_group["virtual_sample"]:
            self.vsample_parameters = sweep.virtual_sample_parameters_sweep(
                **self.params_by_group["virtual_sample"]
            )
        if self.params_by_group["modality"]:
            self.modality_parameters = sweep.acquisition_parameters_sweep(
                **self.params_by_group["modality"]
            )
        if self.params_by_group["acquisition"]:
            self.acquisition_parameters = sweep.acquisition_parameters_sweep(
                **self.params_by_group["acquisition"]
            )

    def run_analysis(self):
        self.analysis["measurements"], self.analysis["inputs"] = sweep.analyse_sweep_single_reference(
            self.acquisition_outputs, 
            self.acquisition_outputs_parameters, 
            self.reference_image[0], 
            self.reference_image_parameters)
    

    def gen_analysis_dataframe(self):
        self.analysis["dframe"], self.analysis["extended_dataframe"] = sweep.measurements_dataframe(
            measurement_vectors=self.analysis["measurements"],
            probe_parameters=self.probe_parameters,
            p_defects=self.defect_parameters,
            sample_params=self.vsample_parameters,
            mod_acq=self.acquisition_parameters,
            mod_names=self.modalities,
            mod_params=self.modality_parameters)

    # methods to retrieve attributes    
    def get_analysis_output(self, keyname="extended_dataframe"):
        return self.analysis[keyname]
            
    def save_analysis(self, keyname="extended_dataframe", output_name="results", output_directory=None):
        if output_directory is None:
            output_directory = self.ouput_directory
        if keyname == "extended_dataframe":
            results = self.get_analysis_output(keyname)
            title = output_name
            file_name = title + "_dataframe.csv"
            results.to_csv(os.path.join(output_directory, file_name), index=False)