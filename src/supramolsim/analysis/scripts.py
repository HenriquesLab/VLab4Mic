from tqdm import tqdm
from supramolsim.experiments import ExperimentParametrisation
import itertools
import numpy as np
from IPython.utils import io
from ..utils.transform.datatype import truncate
import matplotlib.pyplot as plt
import os
import pandas as pd
from .metrics import img_compare
from sklearn import preprocessing as pre


def parameter_sweep_reps(
    Experiment: ExperimentParametrisation,
    sweep_parameters,
    write=False,
    repetitions=1,
    reference_parameters = None,
    **kwargs,
):
    """
    Runs a parameter sweep with multiple repetitions
    for each parameter combination.

    Args:
        Experiment (ExperimentParametrisation):
            The experiment to run parameter sweeps on
        sweep_parameters (dict):
            Dictionary containing parameter names and their sweep ranges
        write (bool, optional):
            Whether to save outputs to disk. Defaults to False.
        repetitions (int, optional):
            Number of repetitions per parameter combination. Defaults to 1.
        **kwargs:
            Additional keyword arguments

    Returns:
        tuple: Contains:
            - sweep_outputs (list):
                List of simulation outputs for each parameter combination
                and repetition
            - iteration_params (list):
                List of parameter dictionaries used for each iteration
            - reference (dict):
                Reference images generated for the experiment

    """
    sweep_outputs = []
    # setup sweep parameters
    for parametername, pars in sweep_parameters.items():
        Experiment.sweep_pars[parametername] = pars
    Experiment._param_linspaces()
    if reference_parameters:
        reference = Experiment.gen_reference(**reference_parameters)
    else:
        reference = Experiment.gen_reference()
    out_dir = Experiment.output_directory
    # prepare combination of parameters
    linspaces_dict = Experiment.sweep_linspaces
    linspace_list = [linspaces_dict[p_name] for p_name in linspaces_dict.keys()]
    total_lengths = [len(i) for i in linspace_list]
    total_len = np.prod(np.array(total_lengths))
    # iterate over parameter combination
    print("Parameter sweep")
    iteration_params = []
    for combination in tqdm(itertools.product(*linspace_list), total=total_len):
        with io.capture_output() as captured:
            labeff = combination[0]
            defect = combination[1]
            iteration_output_reps = []
            for rep in range(repetitions):
                iteration_name = (
                    Experiment.experiment_id
                    + "LEff_"
                    + str(truncate(labeff, 3))
                    + "_Defect_"
                    + str(truncate(defect, 3))
                )
                # start of replicate
                _particle = Experiment._build_particle(
                    lab_eff=labeff, defect=defect, keep=True
                )
                if write:
                    _particle.show_instance(with_sources=True)
                    fig_name = iteration_name + "rep" + str(rep) + ".png"
                    name_path = os.path.join(out_dir, fig_name)
                    plt.savefig(name_path)
                    plt.close()
                exported_field = Experiment._build_coordinate_field(
                    use_self_particle=True
                )
                Experiment.imager.import_field(**exported_field)
                iteration_output = Experiment.run_simulation(
                    name=iteration_name, save=write
                )
                iteration_output_reps.append(iteration_output)
            sweep_outputs.append(iteration_output_reps)
            # end of replicate
            # dictionary keys to use sweep_parameters keys
            iteration_params.append(dict(labelling_efficiency=labeff, defect=defect))
    return sweep_outputs, iteration_params, reference


def _reformat_img_stack(img, use_first=True, **kwargs):
    if use_first:
        if len(img.shape) == 3:
            # print("input is stack")
            single_img = img[0]
        else:
            # print("input is image")
            single_img = img
        formated_im = _crop_image(single_img, **kwargs)
        return formated_im
    else:
        pass


def _crop_image(img, subregion=False, normalise=True, **kwargs):
    if len(img.shape) == 2:
        single_img = img
        if normalise:
            # single_img = pre.MinMaxScaler().fit_transform(single_img)
            single_img = (single_img - single_img.min()) / (
                single_img.max() - single_img.min()
            )
        if subregion:
            single_img = single_img[
                subregion[0] : subregion[1], subregion[0] : subregion[1]
            ]
        return single_img


def analyse_sweep(img_outputs, img_params, reference, analysis_case_params, **kwargs):
    """
    Analyzes parameter sweep results by comparing simulated images to references.

    This function processes the outputs of a parameter sweep simulation,
    comparing each simulated image to corresponding reference images
    using specified metrics. It organizes the results into a structured
    format for analysis.

    Args:
        img_outputs (list):
            List of simulation image outputs
            for each parameter combination and repetition
        img_params (list):
            List of parameter dictionaries used for each iteration
        reference (dict):
            Dictionary of reference images for each experimental condition
        analysis_case_params (dict):
            Parameters for analyzing each experimental condition/case
        **kwargs:
            Additional keyword arguments

    Returns:
        tuple: Contains:
            - data_frame (pd.DataFrame):
                DataFrame containing analysis results with columns for parameters,
              conditions,
                metrics and replica numbers
            - queries (dict):
                Dictionary organizing the analyzed images by parameter combination
                and replica
            - references (dict):
                Dictionary of processed reference images for each condition
    """
    conditions = list(reference.keys())
    references = dict()
    measurement_vectors = list()
    param_names = list(img_params[0].keys())
    n_param_names = len(param_names)
    queries = dict()
    for case in conditions:
        reference_img = _reformat_img_stack(
            reference[case], **analysis_case_params[case]
        )
        references[case] = reference_img
    for i in range(len(img_params)):
        param_values = [truncate(v, 6) for k, v in img_params[i].items()]
        combination_name = str(param_values)
        queries[combination_name] = dict()
        rep = 0
        for img_rep in img_outputs[i]:
            queries[combination_name][rep] = dict()
            for case in conditions:
                item_vector = []
                case_iteration = _reformat_img_stack(
                    img_rep[case], **analysis_case_params[case]
                )
                rep_measurement = img_compare(
                    references[case], case_iteration, **analysis_case_params[case]
                )
                #
                item_vector.append(case)
                item_vector.append(img_params[i]["labelling_efficiency"])
                item_vector.append(img_params[i]["defect"])
                item_vector.append(rep_measurement)
                item_vector.append(rep)
                #
                measurement_vectors.append(item_vector)
                queries[combination_name][rep][case] = case_iteration
            rep = rep + 1
    measurement_array = np.array(measurement_vectors)
    data_frame = pd.DataFrame(
        data={
            "Labelling_efficiency": np.array(measurement_array[:, 1], dtype=np.float32),
            "Fractional_defect": np.array(measurement_array[:, 2], dtype=np.float32),
            "Condition": measurement_array[:, 0],
            "Metric": np.array(measurement_array[:, 3], dtype=np.float32),
            "replica": measurement_array[:, 4],
        }
    )
    return data_frame, queries, references
