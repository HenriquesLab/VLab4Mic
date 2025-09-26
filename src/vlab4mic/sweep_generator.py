from .experiments import ExperimentParametrisation
from .analysis import sweep
import matplotlib.pyplot as plt
from IPython.utils import io
import os
from pathlib import Path
from .utils.io.yaml_functions import load_yaml, save_yaml
from .analysis import _plots
import numpy as np
from datetime import datetime
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import copy
import tifffile as tiff
from pandas.api.types import is_numeric_dtype

#output_dir = Path.home() / "vlab4mic_outputs"


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
    output_directory = None
    reference_virtual_sample = None
    reference_virtual_sample_params = None
    reference_image = None
    reference_image_mask = None
    virtual_samples = None
    virtual_samples_parameters = None
    acquisition_outputs = None
    acquisition_outputs_masks = None
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
    #output_directory: str = None

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
        self.analysis_parameters["metrics_list"] = [
            "ssim",
        ]
        self.plot_parameters = {}
        self.plot_parameters["heatmaps"] = {}
        self.plot_parameters["heatmaps"]["category"] = "modality_name"
        self.plot_parameters["heatmaps"]["param1"] = None
        self.plot_parameters["heatmaps"]["param2"] = None
        self.parameters_with_set_values = []
        self.plot_parameters["lineplots"] = {}
        self.plot_parameters["lineplots"]["x_param"] = None
        self.plot_parameters["lineplots"]["hue"] = "modality_name"
        self.plot_parameters["lineplots"]["style"] = None
        self.plot_parameters["lineplots"]["estimator"] = "mean"
        self.plot_parameters["lineplots"]["errorbar"] = "ci"
        self.structures_info_list = self.experiment.structures_info_list
        # Use the directly loaded parameter_settings instead of experiment.param_settings
        # to ensure all parameter groups (including particle_defect) are available
        self.param_settings = self.parameter_settings
        self.use_experiment_structure = False
        self.reference_parameters_unsorted = dict()
        print("vLab4mic sweep generator initialised")

    def set_number_of_repetitions(self, repeats: int = 3):
        """
        Sets the number of repetitions for each parameter configuration in a sweep.

        Parameters
        ----------
        repeats : int, optional
            Number of times to repeat each parameter combination. Defaults to 3.

        Returns
        -------
        None
        """
        self.sweep_repetitions = repeats

    def select_structures(self, structures: list = None, **kwargs):
        """
        Select structures to use in the sweep.

        Parameters
        ----------
        structures : list of str, optional
            List of 4-letter PDB/CIF IDs or paths to structure files. If None, uses all available structures.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        None
        """
        if structures is not None and type(structures) == list:
            self.structures = structures
            self.use_experiment_structure = False

    def select_probe_templates(self, probe_templates: list = None, **kwargs):
        """
        Select probe templates to use in the sweep.

        Parameters
        ----------
        probe_templates : list of str, optional
            List of probe configuration filenames. If None, uses all available probes.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        None
        """
        if probe_templates is not None and type(probe_templates) == list:
            self.probes = probe_templates

    def select_modalities(self, modalities: list = None, **kwargs):
        """
        Select imaging modalities to use in the sweep.

        Parameters
        ----------
        modalities : list of str, optional
            List of modality names. If None, uses all available modalities.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        None
        """
        if modalities is not None and type(modalities) == list:
            self.modalities = modalities

    def set_output_directory(self, output_directory: str = None):
        """
        Set the output directory for saving results.

        Parameters
        ----------
        output_directory : str, optional
            Path to the desired output directory. If None, uses the default output directory.

        Returns
        -------
        None
        """
        if output_directory is not None:
            self.output_directory = output_directory

    # generators
    def generate_virtual_samples(self):
        """
        Generate virtual samples from the specified parameter combinations.

        Returns
        -------
        None
        """
        self.create_parameters_iterables()
        if self.use_experiment_structure:
                       (
                self.experiment,
                self.virtual_samples,
                self.virtual_samples_parameters,
            ) = sweep.sweep_vasmples(
                experiment=self.experiment,
                structures=self.structures,
                probes=self.probes,
                probe_parameters=self.probe_parameters,
                particle_defects=self.defect_parameters,
                virtual_samples=self.vsample_parameters,
                repetitions=self.sweep_repetitions,
                use_experiment_structure=True
            )
        else: 
            (
                self.experiment,
                self.virtual_samples,
                self.virtual_samples_parameters,
            ) = sweep.sweep_vasmples(
                experiment=self.experiment,
                structures=self.structures,
                probes=self.probes,
                probe_parameters=self.probe_parameters,
                particle_defects=self.defect_parameters,
                virtual_samples=self.vsample_parameters,
                repetitions=self.sweep_repetitions,
            )

    def generate_acquisitions(self):
        """
        Generate image simulation acquisition for all virtual samples generated with generate_virtual_samples.

        Returns
        -------
        None
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
            self.acquisition_outputs_masks
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
        Specify parameters to use for both reference sample and reference image used at analysis.

        Parameters
        ----------
        reference_structure : str, optional
            4-letter ID of PDB/CIF model.
        reference_probe : str, optional
            Name ID of probe configuration file (filename).
        reference_probe_parameters : dict, optional
            Probe parameters to use.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        None
        """
        if reference_structure is None:
            self.reference_structure = self.structures[0]
        else:
            self.reference_structure = reference_structure
        if reference_probe is not None:
            self.reference_probe = reference_probe
        if reference_probe_parameters is not None:
            self.reference_probe_parameters = reference_probe_parameters
        if "relative_positions" in kwargs.keys():
            self.reference_parameters_unsorted["particle_positions"] = kwargs.pop("relative_positions")
        for key, value in kwargs.items():
            self.reference_parameters_unsorted[key] = value

    def generate_reference_sample(self, use_experiment_probe=True, use_experiment_vsample=True):
        """
        Generate the reference virtual sample and its parameters.

        Returns
        -------
        None
        """
        if self.use_experiment_structure:
            structure_is_path = True
        else: 
            structure_is_path = False
        if use_experiment_probe:
            reference_probe = list(self.probes)[0]
        else:
            reference_probe = None
        if use_experiment_vsample:
            reference_sample_params = copy.copy(self.experiment.virtualsample_params)
        self.set_reference_parameters(
            reference_probe=reference_probe,
            **reference_sample_params
            )
        self.reference_virtual_sample, self.reference_virtual_sample_params = (
            sweep.generate_global_reference_sample(
                structure=self.reference_structure,
                probe=self.reference_probe,
                probe_parameters=self.reference_probe_parameters,
                structure_is_path=structure_is_path,
                **self.reference_parameters_unsorted
            )
        )

    def generate_reference_image(self, override=False):
        """
        Generate the reference image using the reference virtual sample.

        Parameters
        ----------
        override : bool, optional
            If True, regenerate the reference image even if it exists. Defaults to False.

        Returns
        -------
        None
        """
        if self.reference_image is None or override:
            if self.use_experiment_structure:
                self.reference_structure = copy.copy(self.structure_path)
                self.generate_reference_sample()
            else:
                self.generate_reference_sample()
        self.reference_image, self.reference_image_parameters, self.reference_image_mask = (
            sweep.generate_global_reference_modality(
                reference_vsample=self.reference_virtual_sample,
                reference_vsample_params=self.reference_virtual_sample_params,
            )
        )

    def load_reference_image(
        self, ref_image_path=None, ref_pixelsize=None, override=True
    ):
        """
        Load a reference image from a specified path.
        """
        reference_parameters = {}
        reference_parameters["Vector"] = None
        ref_pixelsize
        reference_parameters["ref_pixelsize"] = ref_pixelsize
        ref_image = tiff.imread(ref_image_path)
        if override:
            self.reference_image = ref_image
            self.reference_image_parameters = reference_parameters

    # previews
    def preview_image_output_by_ID(
        self,
        probe_template=0,
        probe_parameters=0,
        defect_parameters=0,
        virtual_sample_parameters=0,
        modality_template=0,
        modality_parameters=0,
        acquisition_parameters=0,
        replica_number=0,
        frame=0,
        return_image=False,
        cmap="Grays_r"
    ):
        """
        Preview or return the first acquisition output image.

        Parameters
        ----------
        return_image : bool, optional
            If True, return the image array instead of displaying it. Defaults to False.

        Returns
        -------
        numpy.ndarray or None
            The image array if `return_image` is True, otherwise None.
        """
        parameter_id = (
            str(probe_template)
            + "_"
            + str(probe_parameters)
            + "_"
            + str(defect_parameters)
            + "_"
            + str(virtual_sample_parameters)
            + "_"
            + str(modality_template)
            + "_"
            + str(modality_parameters)
            + "_"
            + str(acquisition_parameters)
        )
        if (
            parameter_id is None
            or parameter_id not in self.acquisition_outputs_parameters.keys()
        ):
            parameter_id = list(self.acquisition_outputs_parameters.keys())[0]
        if (
            len(self.acquisition_outputs[parameter_id][replica_number].shape)
            == 3
        ):
            image = self.acquisition_outputs[parameter_id][replica_number][
                frame
            ]
        else:
            image = self.acquisition_outputs[parameter_id][replica_number]
        if return_image:
            return image, self.acquisition_outputs_parameters[parameter_id]
        else:
            plt.imshow(image, cmap=cmap)
            print(self.acquisition_outputs_parameters[parameter_id])

    def preview_image_output_by_parameter_values(self, **kwargs):
        pass

    def preview_reference_image(self, return_image=False, cmap="Grays_r"):
        """
        Preview or return the reference image.

        Parameters
        ----------
        return_image : bool, optional
            If True, return the image array instead of displaying it. Defaults to False.

        Returns
        -------
        numpy.ndarray or None
            The image array if `return_image` is True, otherwise None.
        """
        if return_image:
            return self.reference_image
        else:
            plt.imshow(self.reference_image[0], cmap=cmap)
            print(self.reference_image_parameters)

    # set and change parameters
    def set_parameter_values(self, param_group, param_name, values=None):
        """
        Specify list or range of values to use per parameter.

        Parameters
        ----------
        param_group : str
            Group of parameter, e.g., 'probe', 'virtual_sample'.
        param_name : str
            Name of parameter to use for sweep.
        values : None, tuple, or list, optional
            Parameter values to use. If tuple, interpreted as (start, stop, num) for linspace.
            If list, used directly. If None, defaults are generated.

        Returns
        -------
        None
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
                if param_group in self.parameter_settings.keys():
                    if (
                        param_name
                        in self.parameter_settings[param_group].keys()
                    ):
                        if (
                            self.parameter_settings[param_group][param_name][
                                "wtype"
                            ]
                            == "int_slider"
                        ):
                            step = np.ceil((values[1] - values[0]) / values[2])
                            self.params_by_group[param_group][param_name] = (
                                np.arange(
                                    start=values[0],
                                    stop=values[1],
                                    step=step,
                                    dtype=int,
                                )
                            )
                        else:
                            param_iterables = np.linspace(
                                values[0], values[1], values[2]
                            )
                            self.params_by_group[param_group][
                                param_name
                            ] = param_iterables
            self.parameters_with_set_values.append(param_name)
        else:
            print(f"{param_group} is not a valid parameter group")

    def clear_sweep_parameters(self):
        """
        Clear all parameters set for the sweep, resetting to default values.

        Returns
        -------
        None
        """
        for group_name in self.params_by_group.keys():
            self.params_by_group[group_name] = {}
        self.parameters_with_set_values = []

    def create_parameters_iterables(self):
        """
        Create iterables for all parameter groups based on set values.

        Returns
        -------
        None
        """
        param_groups = list(self.params_by_group.keys())
        no_params_set = True
        for group_name in param_groups:
            if len(self.params_by_group[group_name]) > 0:
                no_params_set = False
        if no_params_set:
            # set default with minimal options to iterate
            self.set_parameter_values(
                "probe", "labelling_efficiency", values=[0.5, 1]
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

    def set_analysis_parameters(
        self, metrics_list: list[str] = None, zoom_in: float = None, **kwargs
    ):
        """
        Set analysis parameters for the simulation.

        Parameters
        ----------
        metrics_list : list of str, optional
            A list of metric names to be used in the analysis.
        zoom_in : float, optional
            A zoom factor to be applied during analysis.
        **kwargs
            Additional keyword arguments for future extensibility.

        Returns
        -------
        None
        """
        if metrics_list is not None and type(metrics_list) == list:
            self.analysis_parameters["metrics_list"] = metrics_list
        if zoom_in is not None:
            self.analysis_parameters["zoom_in"] = zoom_in

    def set_plot_parameters(self, plot_type, **kwargs):
        """
        Set plotting parameters for a given plot type.

        Parameters
        ----------
        plot_type : str
            The type of plot to configure.
        **kwargs
            Plot-specific keyword arguments.

        Returns
        -------
        None
        """
        if plot_type in self.plot_parameters.keys():
            for key, val in kwargs.items():
                self.plot_parameters[plot_type][key] = val

    def run_analysis(
        self,
        save=True,
        output_name=None,
        output_directory=None,
        plots=False,
        save_images=True,
        **kwargs,
    ):
        """
        Analyse image simulations against the specified image reference.

        This method generates a dataframe containing the calculated metrics per parameter combination.

        Parameters
        ----------
        save : bool, optional
            Specify if output analysis will be saved in output_directory. Defaults to True.
        output_name : str, optional
            Name to identify analysis output files.
        output_directory : str, optional
            Path for writing outputs.
        plots : bool, optional
            Generate heatmaps and lineplots from the results dataframe. Defaults to False.
        save_images : bool, optional
            Whether to save images. Defaults to True.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        None
        """
        if self.acquisition_outputs is None:
            self.generate_acquisitions()
        if self.reference_image is None:
            self.generate_reference_image()
        if len(self.reference_image.shape) == 3:
            # if reference image is 3D, take the first frame
            reference_image = self.reference_image[0]
            reference_image_mask = self.reference_image_mask
        else:
            reference_image = self.reference_image
            reference_image_mask = self.reference_image_mask
        measurement_vectors, inputs, metric = (
            sweep.analyse_sweep_single_reference(
                img_outputs=self.acquisition_outputs,
                img_outputs_masks=self.acquisition_outputs_masks,
                img_params=self.acquisition_outputs_parameters,
                reference_image=reference_image,
                reference_image_mask=reference_image_mask,
                reference_params=self.reference_image_parameters,
                **self.analysis_parameters,
            )
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
                        metric_name=metric_name,
                        filter_dictionary=None,
                    )
        if save:
            self.save_analysis(
                output_name=output_name, output_directory=output_directory
            )

    def gen_analysis_dataframe(self):
        """
        Generate the analysis dataframe from measurement vectors and parameters.

        Returns
        -------
        None
        """
        unformatted_df, self.analysis["dataframes"] = (
            sweep.measurements_dataframe(
                measurement_vectors=self.analysis["unsorted"][
                    "measurement_vectors"
                ],
                probe_parameters=self.probe_parameters,
                p_defects=self.defect_parameters,
                sample_params=self.vsample_parameters,
                mod_acq=self.acquisition_parameters,
                mod_names=self.modalities,
                mod_params=self.modality_parameters,
                metric_names=self.analysis_parameters["metrics_list"],
            )
        )

    # methods to retrieve attributes
    def get_analysis_output(self, keyname="dataframes"):
        """
        Retrieve the analysis output associated with the specified key.

        Parameters
        ----------
        keyname : str, optional
            The key corresponding to the desired analysis output. Defaults to "dataframes".

        Returns
        -------
        Any
            The analysis output associated with the given key.

        Raises
        ------
        KeyError
            If the specified key does not exist in the analysis dictionary.
        """
        return self.analysis[keyname]

    def generate_analysis_plots(
        self,
        plot_type=None,
        metric_name=None,
        decimals: int = None,
        return_figure=True,
        filter_dictionary=None,
        **kwargs,
    ):
        """
        Generates and stores analysis plots based on the specified plot type and metric.

        Parameters
        ----------
        plot_type : str, optional
            The type of plot to generate. Supported types are "heatmaps" and "lineplots".
            If not specified, defaults to "heatmaps".
        metric_name : str, optional
            The name of the metric to plot. Required for generating the plot.
        decimals : int, optional
            Number of decimal places to format metric values. If not specified, defaults to 3.
        **kwargs
            Additional keyword arguments passed to the underlying plotting functions.

        Notes
        -----
        The generated plot is stored in the `self.analysis["plots"]` dictionary under the corresponding plot type and metric name.
        Plot parameters are retrieved from `self.plot_parameters` based on the plot type.
        Requires that analysis output dataframes are available via `self.get_analysis_output(keyname="dataframes")`.

        Returns
        -------
        None
        """
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
        if plot_type == "heatmaps":
            metric_plot = self._gen_heatmaps(
                metric_name=metric_name,
                return_figure=return_figure,
                decimals=decimals,
                filter_dictionary=filter_dictionary,
                **plot_params,
                **kwargs
            )
            self.analysis["plots"][plot_type][metric_name] = metric_plot
        elif plot_type == "lineplots":
            metric_plot = self._gen_lineplots(
                data=data,
                metric_name=metric_name,
                decimals=decimals,
                return_figure=return_figure,
                filter_dictionary=filter_dictionary,
                **plot_params,
                **kwargs
            )
            self.analysis["plots"][plot_type][metric_name] = metric_plot

    def _gen_heatmaps(
        self,
        metric_name=None,
        category: str = None,
        param1: str = None,
        param2: str = None,
        return_figure=False,
        decimals="%.4f",
        filter_dictionary=None,
        **kwargs,
    ):
        """
        Generate heatmap plots for the specified metric and parameters.

        Parameters
        ----------
        metric_name : str, optional
            Name of the metric to plot.
        category : str, optional
            Category for pivoting data (e.g., modality_name).
        param1 : str, optional
            First parameter for the heatmap axes.
        param2 : str, optional
            Second parameter for the heatmap axes.
        return_figure : bool, optional
            If True, return the matplotlib figure. Defaults to False.
        decimals : str, optional
            Format string for decimal places. Defaults to '%.4f'.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        matplotlib.figure.Figure
            The generated heatmap figure.
        """
        if metric_name is None:
            metric_name = self.analysis_parameters["metrics_list"][0]
        if category is None:
            category = "modality_name"
        if param1 is None:
            param1 = self.parameters_with_set_values[0]
        if param2 is None:
            param2 = "probe_n"
        analysis_resut_df = self.get_analysis_output(keyname="dataframes")
        df = copy.deepcopy(self.analysis["dataframes"])
        if is_numeric_dtype(df[param1]):
            df[param1] = df[param1].round(3)
        if is_numeric_dtype(df[param2]):
            df[param2] = df[param2].round(3)
        pre_filter_dt = True
        if filter_dictionary is not None:
             pre_filter_dt = True
        df_categories, titles = sweep.pivot_dataframes_byCategory(
            dataframe=df,
            category_name=category,
            param1=param1,
            param2=param2,
            metric_name=metric_name,
            pre_filter_dt=pre_filter_dt,
            filter_dictionary=filter_dictionary,
            **kwargs
        )
        conditions = list(df_categories.keys())
        nconditions = len(conditions)
        my_palette = sns.diverging_palette(20, 200, s=80, l=50, as_cmap=True)
        plot = _plots.sns_heatmap_pivots(
            df_categories,
            titles,
            annotations=True,
            return_figure=return_figure,
            metric_name=metric_name,
            conditions_cmaps=[my_palette]*nconditions,
            decimals=decimals,
            **kwargs
        )
        return plot

    def _gen_lineplots(
        self,
        data=None,
        x_param=None,
        metric_name=None,
        hue="modality_name",
        style=None,
        markers=None,
        estimator="mean",
        errorbar="ci",
        figsize=[10, 10],
        decimals="%.4f",
        return_figure=True,
        filter_dictionary=None,
        **kwargs,
    ):
        """
        Generate line plots for the specified metric and parameters.

        Parameters
        ----------
        data : pandas.DataFrame, optional
            Data to plot. If None, uses self.analysis["dataframes"].
        x_param : str, optional
            Parameter for the x-axis.
        metric_name : str, optional
            Metric to plot on the y-axis.
        hue : str, optional
            Variable that defines subsets of the data, which will be drawn on separate lines.
        style : str, optional
            Variable that determines the line style.
        markers : bool or list, optional
            Markers for the lines.
        estimator : str or callable, optional
            Method for aggregating data.
        errorbar : str, optional
            Error bar type.
        figsize : list, optional
            Figure size.
        decimals : str, optional
            Format string for decimal places.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        matplotlib.figure.Figure
            The generated line plot figure.
        """
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
        sns.lineplot(
            data=data,
            x=x_param,
            y=metric_name,
            hue=hue,
            style=style,
            markers=markers,
            estimator=estimator,
            errorbar=errorbar,
            ax=axes,
            **kwargs
        )
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
        """
        Saves analysis results, including dataframes and plots, to the specified output directory.

        Parameters
        ----------
        output_name : str, optional
            Base name for the output files. If None, defaults to "vLab4mic_results_".
            The current date (YYYYMMDD) will be appended to the name.
        output_directory : str, optional
            Directory where the output files will be saved. If None, uses the instance's `output_directory` attribute.
        analysis_type : list of str, optional
            Types of analysis outputs to save. Can include "dataframes" and/or "plots".
            If None, both "dataframes" and "plots" are saved.

        Notes
        -----
        Dataframes are saved as CSV files.
        Plots are saved as PNG files, named according to their metric and plot type.
        The method relies on `get_analysis_output` to retrieve the relevant dataframes and plots.

        Returns
        -------
        None
        """
        now = datetime.now()  # dd/mm/YY H:M:S
        dt_string = now.strftime("%Y%m%d")
        if analysis_type is None:
            analysis_type = ["dataframes", "plots"]
        if output_name is None:
            output_name = "vLab4mic_results_"
        output_name = output_name + dt_string
        if output_directory is None:
            output_directory = self.output_directory
        for keyname in analysis_type:
            if keyname == "dataframes":
                df = self.get_analysis_output(keyname)
                df_name = output_name + "_dataframe.csv"
                df.to_csv(os.path.join(output_directory, df_name), index=False)
            elif keyname == "plots":
                plots_dictionary = self.get_analysis_output(keyname)
                for plot_type, plots_by_metric in plots_dictionary.items():
                    for metric, plot in plots_by_metric.items():
                        figure_name = (
                            output_name
                            + "_"
                            + metric
                            + "_"
                            + plot_type
                            + ".png"
                        )
                        plot.savefig(
                            os.path.join(output_directory, figure_name)
                        )

    def save_images(self, output_name=None, output_directory=None, floats_as=float):
        """
        Saves simulated images and an optional reference image to disk in TIFF format.

        Parameters
        ----------
        output_name : str, optional
            Base name for the output image files. Defaults to "vLab4mic_images_".
        output_directory : str, optional
            Directory where images will be saved. If not provided, uses a default path
            based on `self.output_directory` under a "simulated_images" subdirectory. The directory is created if it does not exist.

        Returns
        -------
        None

        Behavior
        --------
        - Iterates over all parameter combinations in `self.acquisition_outputs`.
        - For each combination, concatenates all replicate images and saves the result as a TIFF file named after the parameter combination ID.
        - If `self.reference_image` is present, saves it as "reference.tiff" in the output directory.
        """
        if output_name is None:
            output_name = "vLab4mic_images_"
        if output_directory is None:
            output_directory = os.path.join(
                self.output_directory, "simulated_images", ""
            )
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
        for (
            param_combination_id,
            replicates,
        ) in self.acquisition_outputs.items():
            nreps = len(replicates)
            image = replicates[0]
            for i in range(1, nreps):
                image = np.concatenate((image, replicates[i]))
            name = output_directory + param_combination_id + ".tiff"
            tiff.imwrite(name, image)
        if floats_as is not None and callable(floats_as):
            copy_of_params = copy.deepcopy(self.acquisition_outputs_parameters)
            for combination_id, list_of_parameters in copy_of_params.items():
                for parameter in list_of_parameters:
                    if type(parameter) is dict:
                        for parameter_name in parameter.keys():
                            if type(parameter[parameter_name]) is np.float64:
                                parameter[parameter_name] = floats_as(parameter[parameter_name])
            save_yaml(
                data=copy_of_params,
                name="acquisition_parameters",
                output_directory=output_directory,
            )            
        else:
            save_yaml(
                data=self.acquisition_outputs_parameters,
                name="acquisition_parameters",
                output_directory=output_directory,
            )
        if self.reference_image is not None:
            # save reference image
            name_ref = output_directory + "reference.tiff"
            tiff.imwrite(name_ref, self.reference_image)


def run_parameter_sweep(
        structures: list[str] = None,
        probe_templates: list[str] = None,
        modalities: list[str] = None,
        output_directory: str = None,
        output_name: str = None,
        save_analysis_results: bool = True,
        analysis_plots: bool = True,
        save_sweep_images: bool = True,
        sweep_repetitions: int = 3,
        **kwargs
):
    sweep_gen = sweep_generator()
    sweep_gen.select_structures(structures=structures)
    sweep_gen.select_probe_templates(probe_templates=probe_templates)
    sweep_gen.select_modalities(modalities=modalities)
    sweep_gen.set_output_directory(output_directory=output_directory)
    sweep_gen.set_number_of_repetitions(sweep_repetitions)
    sweep_gen.run_analysis(
        save=save_analysis_results, 
        plots=analysis_plots,
        output_name=output_name
        )
    if save_sweep_images:
        sweep_gen.save_images(
            output_name=output_name, 
            output_directory=output_directory
        )