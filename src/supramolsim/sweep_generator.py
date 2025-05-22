from .experiments import ExperimentParametrisation
from .analysis import sweep
import matplotlib.pyplot as plt
from IPython.utils import io
import os
from pathlib import Path
from .utils.io.yaml_functions import load_yaml
from .analysis import _plots
import numpy as np

output_dir = Path.home() / "vlab4mic_outputs"


class sweep_generator:
    # sweep parameters
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
    # analysis["measurements"] = None
    # analysis["inputs"] = None
    # analysis["extended_dataframe"] = None
    analysis["unsorted"] = {}
    analysis["dataframes"] = None
    analysis["plots"] = {}
    # saving
    output_directory: str = None

    def __init__(self):
        # parameter dictionary
        self.params_by_group = {}
        self.params_by_group["probe"] = {}
        self.params_by_group["virtual_sample"] = {}
        self.params_by_group["particle_defect"] = {}
        self.params_by_group["modality"] = {}
        self.params_by_group["acquisition"] = {}
        self.output_directory = self.experiment.output_directory
        self.configuration_directory = self.experiment.configuration_path
        param_settings_file = os.path.join(
            self.configuration_directory, "parameter_settings.yaml"
        )
        self.parameter_settings = load_yaml(param_settings_file)
        self.analysis_parameters = {}
        self.analysis_parameters["zoom_in"] = 0
        self.analysis_parameters["metrics_list"] = ["ssim", ]
        self.plot_parameters = {}
        self.plot_parameters["heatmaps"] = {}
        self.plot_parameters["heatmaps"]["category"] = "modality_name"
        self.plot_parameters["heatmaps"]["param1"] = None
        self.plot_parameters["heatmaps"]["param2"] = None

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
        **kwargs,
    ):
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
    def set_parameter_values(self, param_group, param_name, values=None):
        if param_group in list(self.params_by_group.keys()):
            if values is None:
                # create defaults or ignore this parameter for it to not be used
                pass
            elif type(values) == list:
                # set list directly as the params to use
                self.params_by_group[param_group][param_name] = values
            elif type(values) == tuple:
                # 3 values are expected: min, max, steps
                # generate a linspace
                param_iterables = np.linspace(values[0], values[1], values[2])
                self.params_by_group[param_group][param_name] = param_iterables
        else:
            print(f"{param_group} is not a valid parameter group")

    def create_parameters_iterables(self):
        param_groups = list(self.params_by_group.keys())
        no_params_set = True
        for group_name in param_groups:
            if len(self.params_by_group[group_name]) > 0:
                no_params_set = False
        if no_params_set:
            # set default with minimal options to iterate
            self.set_parameter_values(
                "probe", "labelling_efficiency", values=(0.5, 1, 2)
            )

        self.probe_parameters = sweep.create_param_combinations(
            **self.params_by_group["probe"]
        )
        self.defect_parameters = sweep.create_param_combinations(
            **self.params_by_group["particle_defect"]
        )
        self.vsample_parameters = sweep.create_param_combinations(
            **self.params_by_group["virtual_sample"]
        )
        self.acquisition_parameters = sweep.create_param_combinations(
            **self.params_by_group["acquisition"]
        )
        self.modality_parameters = sweep.create_param_combinations(
            **self.params_by_group["modality"]
        )

    def set_analysis_parameters(self, param_name, value):
        if param_name in list(self.analysis_parameters.keys()):
            self.analysis_parameters[param_name] = value
        else:
            print(f"input name {param_name} is not a valid analysis parameter")

    def set_plot_parameters(self, plot_type, **kwargs):
        if plot_type in self.plot_parameters.keys():
            for key, val in kwargs.items():
                self.plot_parameters[plot_type][key] = val


    def run_analysis(
        self, save=False, output_name=None, output_directory=None, plots=False
    ):
        if self.acquisition_outputs is None:
            self.generate_acquisitions()
        if self.reference_image is None:
            self.generate_reference_image()
        measurement_vectors, inputs, metric = sweep.analyse_sweep_single_reference(
            img_outputs=self.acquisition_outputs,
            img_params=self.acquisition_outputs_parameters,
            reference_image=self.reference_image[0],
            reference_params=self.reference_image_parameters,
            **self.analysis_parameters
        )
        self.analysis["unsorted"]["measurement_vectors"] = measurement_vectors
        self.analysis["unsorted"]["inputs"] = inputs
        self.gen_analysis_dataframe()
        if plots:
            for metric_name in self.analysis_parameters["metrics_list"]:
                self.generate_analysis_plots(
                    plot_type="heatmaps",
                    return_figure=True,
                    metric_name=metric_name)
        if save:
            self.save_analysis(
                output_name=output_name, output_directory=output_directory
            )

    def gen_analysis_dataframe(self):
        unformatted_df, self.analysis["dataframes"] = sweep.measurements_dataframe(
            measurement_vectors=self.analysis["unsorted"]["measurement_vectors"],
            probe_parameters=self.probe_parameters,
            p_defects=self.defect_parameters,
            sample_params=self.vsample_parameters,
            mod_acq=self.acquisition_parameters,
            mod_names=self.modalities,
            mod_params=self.modality_parameters,
            metric_names=self.analysis_parameters["metrics_list"]
        )

    # methods to retrieve attributes
    def get_analysis_output(self, keyname="dataframes"):
        return self.analysis[keyname]

    def generate_analysis_plots(
        self,
        plot_type=None,
        metric_name = None,
        **kwargs,
    ):
        if plot_type is None:
            plot_type = "heatmaps"
        plot_params = self.plot_parameters[plot_type]
        if plot_type == "heatmaps":
            metric_plot = self._gen_heatmaps(
                metric_name=metric_name,
                return_figure=True,
                **plot_params)
            self.analysis["plots"][metric_name] = metric_plot
        else:
            pass
        
    def _gen_heatmaps(
            self,
            metric_name=None,
            category: str = None,
            param1: str = None,
            param2: str = None,
            return_figure = False,
            **kwargs):
        print("into heatmaps")
        if metric_name is None:
            metric_name = self.analysis_parameters["metrics_list"][0]
        if category is None:
            category = "modality_name"
        if param1 is None:
            param1 = "labelling_efficiency"
        if param2 is None:
            param2 = "probe_n"
        analysis_resut_df = self.get_analysis_output(keyname="dataframes")
        df_categories, titles = sweep.pivot_dataframes_byCategory(
            dataframe=analysis_resut_df,
            category_name=category,
            param1=param1,
            param2=param2,
            metric_name=metric_name
        )
        plot = _plots.sns_heatmap_pivots(
                df_categories, titles, annotations=True, return_figure=return_figure
            )
        return plot

    def save_analysis(
        self, output_name=None, output_directory=None, analysis_type=None
    ):
        if analysis_type is None:
            analysis_type = ["dataframes", "plots"]
        if output_name is None:
            output_name = "vLab4mic_results"
        if output_directory is None:
            output_directory = self.ouput_directory
        for keyname in analysis_type:
            if keyname == "dataframes":
                df = self.get_analysis_output(keyname)
                df_name = output_name + "_dataframe.csv"
                df.to_csv(os.path.join(output_directory, df_name), index=False)
            elif keyname == "plots":
                plots_dictionary = self.get_analysis_output(keyname)
                for metric, plot in plots_dictionary.items():
                    figure_name = output_name + "_" + metric + "_figure.png"
                    plot.savefig(os.path.join(output_directory, figure_name))
