from .experiments import ExperimentParametrisation
from .analysis import sweep
import matplotlib.pyplot as plt
from IPython.utils import io
import os
from pathlib import Path
from .utils.io.yaml_functions import load_yaml
from .analysis import _plots
import numpy as np
from datetime import datetime
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import copy
from pandas.api.types import is_numeric_dtype

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
        self.parameters_with_set_values = []
        self.plot_parameters["lineplots"] = {}
        self.plot_parameters["lineplots"]["x_param"] = None
        self.plot_parameters["lineplots"]["hue"] = 'modality_name'
        self.plot_parameters["lineplots"]["style"] = None
        self.plot_parameters["lineplots"]["estimator"] = "mean"
        self.plot_parameters["lineplots"]["errorbar"] = "ci"


    def set_number_of_repetitions(self, repeats: int = 3 ):
        self.sweep_repetitions = repeats

    # generators
    def generate_virtual_samples(self):
        """
        Generate virtual samples from the specified parameter combinations

        """
        self.create_parameters_iterables()
        self.experiment, self.virtual_samples, self.virtual_samples_parameters = (
            sweep.sweep_vasmples(
                structures=self.structures,
                probes=self.probes,
                probe_parameters=self.probe_parameters,
                particle_defects=self.defect_parameters,
                virtual_samples=self.vsample_parameters,
                repetitions=self.sweep_repetitions,
            )
        )

    def generate_acquisitions(self):
        """
        Generate image simulation acquisition for all virtual samples
        generated with generate_virtual_samples
        """
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
        """
        Specify parameters to use for both reference sample and 
        reference image used at analysis.

        :param reference_structure: 4-letter ID of PDB/CIF model.
        :type reference_structure: str
        :param reference_probe: Name ID of probe configuration file (filename).
        :type reference_probe: str
        :param reference_probe_parameters: probe parameters to use
        :type reference_probe: dict
        
        """
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
        """
        Specify list or range of values to use per parameter

        :param param_group: group of parameter e.g probe, virtual_sample
        :type save: str
        :param param_name: name of parameter to use for sweep
        :type save: str
        :param values: parameter values to use. 
            If values is a tuple,
            a sequence of values will be generated asssuming the tuple contains
            the parameter for numpy linspace method (start, stop, num)
            If values is a list, the list is set as the final values to use
            If values is a None, a default list will be generated 
            according to the parameter set
        :type values: [None, tuple, list]

        """
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
            self.parameters_with_set_values.append(param_name)
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

    def set_analysis_parameters(self, metrics_list: list[str] = None, zoom_in: float = None, **kwargs):
        if metrics_list is not None and type(metrics_list) == list:
            self.analysis_parameters["metrics_list"] = metrics_list
        if zoom_in is not None:
            self.analysis_parameters["zoom_in"] = zoom_in

    def set_plot_parameters(self, plot_type, **kwargs):
        if plot_type in self.plot_parameters.keys():
            for key, val in kwargs.items():
                self.plot_parameters[plot_type][key] = val


    def run_analysis(
        self, save=True, output_name=None, output_directory=None, plots=False
    ):
        """
        Analyse image simulations against the specified image reference.
        This method generates a dataframe containing the calculated metrics
        per parameter combination

        :param save: specify if output analysis will be saved in output_directory
        :type save: bool
        :param output_name: name to identify analysis output files
        :type save: str
        :param output_directory: path for writing outputs
        :type output_directory: str
        :param plots: Generate heatmas and lineplots from the results dataframe 
        :type plots: bool
        """
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
                for plot_type in self.plot_parameters.keys():
                    self.generate_analysis_plots(
                        plot_type=plot_type,
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
        decimals:int = None,
        **kwargs,
    ):
        if decimals is None:
            decimals = "%.3f"
        else:
            decimals = "%." + str(decimals) + "f"
        data = self.get_analysis_output(keyname="dataframes")
        if plot_type is None:
            plot_type = "heatmaps"
        if plot_type not in self.analysis["plots"].keys():
            self.analysis["plots"][plot_type] = {}
        plot_params = self.plot_parameters[plot_type]
        print(plot_type)
        print(plot_params)
        if plot_type == "heatmaps":
            metric_plot = self._gen_heatmaps(
                metric_name=metric_name,
                return_figure=True,
                decimals=decimals,
                **plot_params)
            self.analysis["plots"][plot_type][metric_name] = metric_plot
        elif plot_type == "lineplots":
            metric_plot = self._gen_lineplots(
                data=data,
                metric_name=metric_name,
                decimals=decimals,
                **plot_params 
                )
            self.analysis["plots"][plot_type][metric_name] = metric_plot
        
    def _gen_heatmaps(
            self,
            metric_name=None,
            category: str = None,
            param1: str = None,
            param2: str = None,
            return_figure = False,
            decimals='%.4f',
            **kwargs):
        if metric_name is None:
            metric_name = self.analysis_parameters["metrics_list"][0]
        if category is None:
            category = "modality_name"
        if param1 is None:
            param1 = "labelling_efficiency"
        if param2 is None:
            param2 = "probe_n"
        analysis_resut_df = self.get_analysis_output(keyname="dataframes")
        df = copy.deepcopy(self.analysis["dataframes"])
        if is_numeric_dtype(df[param1]):
            df[param1] = df[param1].round(3)
        if is_numeric_dtype(df[param2]):
            df[param2] = df[param2].round(3)
        df_categories, titles = sweep.pivot_dataframes_byCategory(
            dataframe=df,
            category_name=category,
            param1=param1,
            param2=param2,
            metric_name=metric_name
        )
        plot = _plots.sns_heatmap_pivots(
                df_categories, titles, 
                annotations=True, 
                return_figure=return_figure, 
                metric_name=metric_name, 
                decimals=decimals
            )
        return plot
    
    def _gen_lineplots(self,
                       data=None,
                       x_param=None,
                       metric_name=None,
                       hue='modality_name',
                       style=None,
                       markers=None, 
                       estimator = "mean", 
                       errorbar="ci",
                       figsize=[10,10],
                       decimals='%.4f',
                       **kwargs):
        if data is None:
            data = self.analysis["dataframes"]
        if x_param is None and len(self.parameters_with_set_values) > 0:
            x_param = self.parameters_with_set_values[0]
        else: 
            x_param = "labelling_efficiency"
        if metric_name is None:
            metric_name = self.analysis_parameters["metrics_list"][0]
        if style is None and len(self.parameters_with_set_values) > 1:
            style = self.parameters_with_set_values[1]
        fig, axes = plt.subplots(figsize=figsize)
        sns.lineplot(data=data,
                     x=x_param,
                     y=metric_name,
                     hue=hue,
                     style=style,
                     markers=markers, 
                     estimator=estimator,
                     errorbar=errorbar,
                     ax=axes)
        axes.yaxis.set_major_formatter(FormatStrFormatter(decimals))
        axes.xaxis.set_major_formatter(FormatStrFormatter(decimals))
        title = estimator + " " + metric_name + " for " + x_param
        if style is not None:
            title = title + "per " + style
        plt.title(title)
        plt.close()  
        return fig

    def save_analysis(
        self, output_name=None, output_directory=None, analysis_type=None
    ):
        now = datetime.now()  # dd/mm/YY H:M:S
        dt_string = now.strftime("%Y%m%d")
        if analysis_type is None:
            analysis_type = ["dataframes", "plots"]
        if output_name is None:
            output_name = "vLab4mic_results_"
        output_name = output_name + dt_string
        if output_directory is None:
            output_directory = self.ouput_directory
        for keyname in analysis_type:
            if keyname == "dataframes":
                df = self.get_analysis_output(keyname)
                df_name = output_name + "_dataframe.csv"
                df.to_csv(os.path.join(output_directory, df_name), index=False)
            elif keyname == "plots":
                plots_dictionary = self.get_analysis_output(keyname)
                for plot_type, plots_by_metric in plots_dictionary.items():
                    for metric, plot in plots_by_metric.items():
                        figure_name = output_name + "_" + metric + "_" + plot_type + ".png"
                        plot.savefig(os.path.join(output_directory, figure_name))
