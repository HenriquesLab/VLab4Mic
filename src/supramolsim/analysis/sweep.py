from ..experiments import ExperimentParametrisation
from ..utils.io.yaml_functions import load_yaml
import os
import supramolsim
from . import metrics
import numpy as np
import pandas as pd
import itertools
import copy
from tqdm.notebook import tqdm
from IPython.utils import io


def sweep_vasmples(
    experiment: ExperimentParametrisation = None,
    structures=None,
    probes=None,
    probe_parameters=None,
    particle_defects=None,
    virtual_samples=None,
    repetitions=1,
    **kwargs,
):
    """
    Generate virtual samples for all combinations of structures, probes, probe parameters, defects, and virtual sample parameters.

    Parameters
    ----------
    experiment : ExperimentParametrisation, optional
        Experiment object to use. If None, a new one is created.
    structures : list, optional
        List of structure IDs.
    probes : list, optional
        List of probe names.
    probe_parameters : dict, optional
        Dictionary of probe parameter sets.
    particle_defects : dict, optional
        Dictionary of defect parameter sets.
    virtual_samples : dict, optional
        Dictionary of virtual sample parameter sets.
    repetitions : int, optional
        Number of repetitions for each combination.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    experiment : ExperimentParametrisation
        The experiment object used.
    vsample_outputs : dict
        Dictionary of generated virtual samples.
    vsample_params : dict
        Dictionary of parameter combinations for each sample.
    """
    # empty lists, but fill up with default
    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    local_dir = os.path.join(pck_dir, "configs")
    default_probe = "NHS_ester"
    default_fluorophore = "AF647"
    if experiment is None:
        experiment = ExperimentParametrisation()
    if structures is None:
        structures = [
            "1XI5",
        ]
    if probes is None:
        print("Probe list is None. Using default probe")
        probes = [
            default_probe,
        ]
    if probe_parameters is None:
        probe_filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        default_params = load_yaml(probe_filepath)
        probe_parameters = dict()
        probe_parameters[0] = None
    if particle_defects is None:
        particle_defects = dict()
        particle_defects[0] = {"use_defects": False}
    if virtual_samples is None:
        # vprobe_filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        # default_vsample =  load_yaml(vprobe_filepath)
        virtual_samples = dict()
        virtual_samples[0] = {"random_orientations": False}
        #virtual_samples = [
        #    None,
        #]
    vsample_params = dict()
    vsample_outputs = dict()

    for struct in structures:
        experiment.structure_id = struct
        experiment._build_structure()
        #for rep in range(repetitions):
        probe_n = 0
        for probe in probes:
            print(f"probe: {probe}")
            # probe_param_n = 0
            for probe_param_n, p_param in probe_parameters.items():
                # copy p_param?
                experiment.remove_probes()
                if p_param is None:
                    experiment.add_probe(
                        probe_template=probe,
                    )
                    #experiment.probe_parameters[probe] = default_params
                else:
                    p_param_copy = copy.deepcopy(p_param)
                    p_param_copy["fluorophore_id"] = default_fluorophore
                    experiment.add_probe(
                        probe_template=probe,
                        **p_param_copy
                    )
                for defect_n, defects_pars in particle_defects.items():
                    particle_defects_copy = copy.deepcopy(defects_pars)
                    for d_key, d_val in particle_defects_copy.items():
                        if d_key=="defect_small_cluster":
                            experiment.defect_eps["eps1"] = d_val
                        if d_key=="defect_large_cluster":
                            experiment.defect_eps["eps2"] = d_val
                        if d_key=="defect":
                            experiment.defect_eps["defect"] = d_val
                    # print(experiment.defect_eps)
                    print(experiment.probe_parameters)
                    experiment._build_particle(keep=True)
                    if experiment.generators_status("particle"):
                        if len(experiment.particle.emitters) == 0:
                            print(f"Skipping {probe}. No emitters were generated")
                            break
                    else:
                        break
                    #vsample_n = 0
                    for vsample_n, vsample_pars in virtual_samples.items():
                        _exported_field = None
                        # combination += str(vsample_n)
                        experiment.set_virtualsample_params(
                            **vsample_pars
                        )
                        #_exported_field = experiment._build_coordinate_field(
                        #    keep=False, use_self_particle=True, **vsample_pars
                        #)
                        for rep in range(repetitions):
                            experiment.clear_virtual_sample()
                            experiment.build(modules=["coordinate_field"], use_self_particle=True)
                            _exported_field = experiment.coordinate_field.export_field()
                            combination_n = (
                                str(probe_n)
                                + "_"
                                + str(probe_param_n)
                                + "_"
                                + str(defect_n)
                                + "_"
                                + str(vsample_n)
                            )
                            _parameters = [
                                struct,
                                probe,
                                p_param,
                                defects_pars,
                                vsample_pars,
                            ]
                            if combination_n not in vsample_params.keys():
                                vsample_params[combination_n] = _parameters
                                vsample_outputs[combination_n] = []
                            vsample_outputs[combination_n].append(_exported_field)
                        #
                        #
                        #vsample_n += 1
                    # defect_n += 1
            probe_n += 1
    return experiment, vsample_outputs, vsample_params


def sweep_modalities_updatemod(
    experiment: ExperimentParametrisation = None,
    vsample_outputs=None,
    vsampl_pars=None,
    modalities: list = None,
    modality_params: dict = None,
    modality_acq_prams: dict =None,
):
    """
    Simulate image acquisitions for all virtual samples and modalities.

    Parameters
    ----------
    experiment : ExperimentParametrisation, optional
        Experiment object to use. If None, a new one is created.
    vsample_outputs : dict, optional
        Dictionary of virtual sample outputs.
    vsampl_pars : dict, optional
        Dictionary of virtual sample parameters.
    modalities : list, optional
        List of modality names.
    modality_params : dict, optional
        Dictionary of modality parameter sets.
    modality_acq_prams : dict, optional
        Dictionary of modality acquisition parameter sets.

    Returns
    -------
    experiment : ExperimentParametrisation
        The experiment object used.
    mod_outputs : dict
        Dictionary of simulated image outputs.
    mod_params : dict
        Dictionary of parameter combinations for each modality.
    pixelsizes : dict
        Dictionary of pixel sizes for each modality.
    """
    default_mod = "Confocal"
    default_aqc = dict(
        nframes=2,
        exp_time=0.005
    )
    default_vsample = "None"
    mod_outputs = dict()
    mod_params = dict()
    if experiment is None:
        experiment = ExperimentParametrisation()
    if vsample_outputs is None:
        vsample_outputs = [
            default_vsample,
        ]
        # needs default sample. can be a minimal field with a single emitter
    if modalities == "all":
        list_of_locals = ["Widefield", "Confocal", "SMLM", "STED"]
        for local_mode in list_of_locals:
            experiment.add_modality(local_mode)
    elif modalities is None:
        experiment.add_modality("STED")
    else:
        for modality_name in modalities:
            experiment.add_modality(modality_name)
    if modality_params is None:
        modality_params = {}
        modality_params[0] = {}
    if modality_acq_prams is None:
        modality_acq_prams = {}
        modality_acq_prams[0] = None
    experiment._build_imager(use_local_field=False, prints=False)
    # print(experiment.objects_created["imager"])
    pixelsizes = dict()
    imager_scale = experiment.imager.roi_params["scale"]
    scalefactor = np.ceil(imager_scale / 1e-9)  # in nanometers
    for mod_name, parameters in experiment.imaging_modalities.items():
        pixelsizes[mod_name] = np.array(
            parameters["detector"]["pixelsize"] * scalefactor
        )
    for vsmpl_id in tqdm(
        list(vsampl_pars.keys()),
        position=0,
        leave=True,
        desc="Unique parameter combination",
    ):
        for virtualsample in tqdm(
            list(vsample_outputs[vsmpl_id]), position=1, leave=False, desc="Repeats"
        ):
            experiment.imager.import_field(**virtualsample)
            with io.capture_output() as captured:
                mod_n = 0
                for modality_name in experiment.selected_mods.keys():
                    for mod_pars_number, mod_pars in modality_params.items():
                        experiment.update_modality(modality_name, **mod_pars)
                        for mod_acq_number, acq_pars in modality_acq_prams.items():
                            if acq_pars:
                                experiment.set_modality_acq(
                                    modality_name=modality_name,
                                    **acq_pars)
                            # iteration_name = combination
                            modality_timeseries = experiment.run_simulation(name="", save=False, modality=modality_name)
                            mod_comb = vsmpl_id + "_" + str(mod_n) + "_" + str(mod_pars_number) + "_" + str(mod_acq_number)
                            mod_parameters = copy.copy(vsampl_pars[vsmpl_id])
                            mod_parameters.append(modality_name)
                            pxsize = experiment.imager.modalities[modality_name]["detector"]["pixelsize"]*1000
                            mod_params_copy = copy.deepcopy(mod_pars)
                            mod_params_copy["pixelsize"] = pxsize
                            mod_parameters.append(mod_params_copy)
                            mod_parameters.append(acq_pars)
                            if mod_comb not in mod_params.keys():
                                mod_params[mod_comb] = mod_parameters
                                mod_outputs[mod_comb] = []
                            mod_outputs[mod_comb].append(modality_timeseries[modality_name])
                            mod_parameters = None
                    mod_n += 1   
    return experiment, mod_outputs, mod_params, pixelsizes


def generate_global_reference_sample(
    experiment: ExperimentParametrisation = None,
    structure=None,
    probe=None,
    probe_parameters=None,
    virtual_sample=None,
):
    """
    Generate a global reference virtual sample.

    Parameters
    ----------
    experiment : ExperimentParametrisation, optional
        Experiment object to use. If None, a new one is created.
    structure : str, optional
        Structure ID.
    probe : str, optional
        Probe name.
    probe_parameters : dict, optional
        Probe parameters.
    virtual_sample : dict, optional
        Virtual sample parameters.

    Returns
    -------
    reference_vsample : dict
        The generated reference virtual sample.
    refernece_parameters : list
        List of parameters used for the reference.
    """
    # empty lists, but fill up with default
    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    local_dir = os.path.join(pck_dir, "configs")
    default_reference_probe = "NHS_ester"
    default_reference_fluorophore = "AF647"
    if experiment is None:
        experiment = ExperimentParametrisation()
    if probe is None:
        probe = default_reference_probe
    if probe_parameters is None:
        probe_filepath = os.path.join(
            local_dir, "probes", default_reference_probe + ".yaml"
        )
        probe_parameters = load_yaml(probe_filepath)
        probe_parameters["fluorophore_id"] = default_reference_fluorophore
    else:
        probe_parameters["fluorophore_id"] = default_reference_fluorophore
    experiment.structure_id = structure
    experiment._build_structure()
    experiment.add_probe(
        probe_template=probe,
        **probe_parameters
    )
    #experiment.structure_label = probe
    #experiment.probe_parameters[probe] = probe_parameters
    experiment._build_particle(keep=True)
    # combination += str(vsample_n)
    refernece_parameters = [structure, probe, probe_parameters, virtual_sample]
    reference_vsample = experiment._build_coordinate_field(
        keep=False, use_self_particle=True
    )
    return reference_vsample, refernece_parameters


def generate_global_reference_modality(
    experiment: ExperimentParametrisation = None,
    reference_vsample=None,
    reference_vsample_params=None,
    modality=None,
    modality_acquisition = None,
):
    """
    Generate a global reference image for a given virtual sample and modality.

    Parameters
    ----------
    experiment : ExperimentParametrisation, optional
        Experiment object to use. If None, a new one is created.
    reference_vsample : dict, optional
        Reference virtual sample.
    reference_vsample_params : dict, optional
        Parameters for the reference virtual sample.
    modality : str, optional
        Modality name.
    modality_acquisition : dict, optional
        Acquisition parameters for the modality.

    Returns
    -------
    reference_output : numpy.ndarray
        The simulated reference image.
    reference_parameters : dict
        Parameters used for the reference image.
    """
    if experiment is None:
        experiment = ExperimentParametrisation()
    if modality is None:
        modality = "Reference"
        #modality_acquisition = None
    if modality_acquisition is None:
        experiment.add_modality(modality)
        
        #experiment.selected_mods[modality] = (
        #    configuration_format.format_modality_acquisition_params()
        #)
    else:
        experiment.add_modality(modality)
        experiment.set_modality_acq(modality, **modality_acquisition)
        #experiment.selected_mods[modality] = modality_acquisition
    experiment._build_imager(use_local_field=False)
    print(experiment.generators_status)
    experiment.imager.import_field(**reference_vsample)
    reference_parameters = dict()
    reference_parameters["Vector"] = [
        reference_vsample_params,
        modality,
        modality_acquisition,
    ]
    reference_output = experiment.run_simulation(name="", save=False)
    imager_scale = experiment.imager.roi_params["scale"]
    scalefactor = np.ceil(imager_scale / 1e-9)  # resulting pixel size in nanometers
    reference_parameters["ref_pixelsize"] = (
        experiment.imager.modalities[modality]["detector"]["pixelsize"] * scalefactor
    )
    return reference_output[modality], reference_parameters


def analyse_image_sweep(img_outputs, img_params, reference, analysis_case_params=None):
    """
    Analyse a sweep of images against a reference image.

    Parameters
    ----------
    img_outputs : dict
        Dictionary of simulated image outputs.
    img_params : dict
        Dictionary of image parameters.
    reference : numpy.ndarray
        Reference image.
    analysis_case_params : dict, optional
        Additional parameters for analysis.

    Returns
    -------
    measurement_vectors : list
        List of measurement results for each image.
    inputs : dict
        Dictionary of input images and used references.
    """
    measurement_vectors = []
    # ref_pixelsize = analysis_case_params["ref_pixelsize"]
    inputs = dict()
    for params_id in img_params.keys():
        inputs[params_id] = dict()
        rep_number = 0
        mod_name = img_params[params_id][5]  # 5th item corresponds to Modality
        for img_r in img_outputs[params_id]:
            im1 = img_r[0]
            im_ref = reference[0]
            rep_measurement, ref_used, qry_used = metrics.img_compare(
                im_ref, im1, **analysis_case_params[mod_name]
            )
            measurement_vectors.append([params_id, rep_number, rep_measurement])
            inputs[params_id][rep_number] = [qry_used, im1]
            rep_number += 1
    return measurement_vectors, inputs

def analyse_sweep_single_reference(img_outputs, img_params, reference_image, reference_params, zoom_in=0, metrics_list:list =["ssim",], **kwargs):
    """
    Analyse a sweep of images against a single reference image using specified metrics.

    Parameters
    ----------
    img_outputs : dict
        Dictionary of simulated image outputs.
    img_params : dict
        Dictionary of image parameters.
    reference_image : numpy.ndarray
        Reference image.
    reference_params : dict
        Parameters for the reference image.
    zoom_in : float, optional
        Zoom factor for analysis. Defaults to 0.
    metrics_list : list of str, optional
        List of metric names to compute. Defaults to ["ssim"].
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    measurement_vectors : list
        List of measurement results for each image.
    inputs : dict
        Dictionary of input images and used references.
    metrics_list : list
        List of metric names used.
    """
    measurement_vectors = []
    # ref_pixelsize = analysis_case_params["ref_pixelsize"]
    inputs = dict()
    for params_id in img_params.keys():
        inputs[params_id] = dict()
        rep_number = 0
        mod_name = img_params[params_id][5]  # 5th item corresponds to Modality
        modality_pixelsize = img_params[params_id][6]["pixelsize"]
        for img_r in img_outputs[params_id]:
            im1 = img_r[0]
            im_ref = reference_image
            rep_measurement, ref_used, qry_used = metrics.img_compare(
                im_ref, im1,
                modality_pixelsize = modality_pixelsize,
                ref_pixelsize = reference_params["ref_pixelsize"],
                force_match=True,
                zoom_in=zoom_in,
                metric=metrics_list
            )
            r_vector = list([params_id, rep_number]) + list([*rep_measurement])
            measurement_vectors.append(r_vector)
            #measurement_vectors = measurement_vectors + rep_measurement[0]
            inputs[params_id][rep_number] = [qry_used, im1, ref_used]
            rep_number += 1
    return measurement_vectors, inputs, metrics_list




def measurements_dataframe(
    measurement_vectors, probe_parameters=None, p_defects=None, mod_names = None, sample_params=None, mod_params=None, mod_acq=None, metric_names=None
):
    """
    Create a pandas DataFrame from measurement vectors and parameter dictionaries.

    Parameters
    ----------
    measurement_vectors : list
        List of measurement results.
    probe_parameters : dict, optional
        Dictionary of probe parameter sets.
    p_defects : dict, optional
        Dictionary of defect parameter sets.
    mod_names : list, optional
        List of modality names.
    sample_params : dict, optional
        Dictionary of sample parameter sets.
    mod_params : dict, optional
        Dictionary of modality parameter sets.
    mod_acq : dict, optional
        Dictionary of modality acquisition parameter sets.
    metric_names : list, optional
        List of metric names.

    Returns
    -------
    data_frame : pandas.DataFrame
        DataFrame with basic combination and replica info.
    df_combined : pandas.DataFrame
        DataFrame with expanded parameter and metric columns.
    """
    measurement_array = np.array(measurement_vectors)
    nrows = len(measurement_array[:, 1])
    ids = [i.split("_") for i in measurement_array[:, 0]]
    ids_array = np.array(ids)
    data_frame = pd.DataFrame(
        data={
            "Combination_id": measurement_array[:, 0],
            "probe_n": ids_array[:, 0],
            "probe_param_n": ids_array[:, 1],
            "defects": ids_array[:, 2],
            "vsample": ids_array[:, 3],
            "modality": ids_array[:, 4],
            "modality_parameters": ids_array[:, 5],
            "modality_acqusition": ids_array[:, 6],
            "Replica": measurement_array[:, 1]
        }
    )
    df_combined = data_frame
    # 
    nmetrics = len(metric_names)
    metrics_dictionary = dict()
    for metric_number in range(nmetrics):
        metricvector = measurement_array[:, 2+metric_number]
        metrics_dictionary[metric_names[metric_number]] = np.array(metricvector, dtype=np.float64)
    metrics_df = pd.DataFrame(metrics_dictionary)
    df_combined = df_combined.join(metrics_df)
    if probe_parameters:
        probe_param_names = probe_parameters[0].keys()
        tmp_df1 = dict()
        for column_name in probe_param_names:
            tmp_df1[column_name] = []
        for i in range(nrows):
            probe_par_comb_id = int(data_frame.iloc[i]["probe_param_n"])
            for column_name in probe_param_names:
                tmp_df1[column_name].append(
                    probe_parameters[probe_par_comb_id][column_name]
                )
        tmp1 = pd.DataFrame(tmp_df1)
        df_combined = df_combined.join(tmp1)
    if p_defects:
        defect_param_names = p_defects[0].keys()
        tmp_df2 = dict()
        for column_name in defect_param_names:
            tmp_df2[column_name] = []
        for i in range(nrows):
            defect_par_comb_id = int(data_frame.iloc[i]["defects"])
            for column_name in defect_param_names:
                tmp_df2[column_name].append(p_defects[defect_par_comb_id][column_name])
        tmp2 = pd.DataFrame(tmp_df2)
        df_combined = df_combined.join(tmp2)
    if sample_params:
        sample_param_names = sample_params[0].keys()
        tmp_df3 = dict()
        for column_name in sample_param_names:
            tmp_df3[column_name] = []
        for i in range(nrows):
            sample_par_comb_id = int(data_frame.iloc[i]["vsample"])
            for column_name in sample_param_names:
                tmp_df3[column_name].append(sample_params[sample_par_comb_id][column_name])
        tmp3 = pd.DataFrame(tmp_df3)
        df_combined = df_combined.join(tmp3)
    if mod_params:
        param_names = mod_params[0].keys()
        tmp_df4 = dict()
        for column_name in param_names:
            tmp_df4[column_name] = []
        for i in range(nrows):
            acq_par_comb_id = int(data_frame.iloc[i]["modality_parameters"])
            for column_name in param_names:
                tmp_df4[column_name].append(mod_params[acq_par_comb_id][column_name])
        tmp4 = pd.DataFrame(tmp_df4)
        df_combined = df_combined.join(tmp4)
    if mod_acq:
        param_names = mod_acq[0].keys()
        tmp_df5 = dict()
        for column_name in param_names:
            tmp_df5[column_name] = []
        for i in range(nrows):
            acq_par_comb_id = int(data_frame.iloc[i]["modality_acqusition"])
            for column_name in param_names:
                tmp_df5[column_name].append(mod_acq[acq_par_comb_id][column_name])
        tmp5 = pd.DataFrame(tmp_df5)
        df_combined = df_combined.join(tmp5)
    if mod_names:
        tmp_df6 = {}
        tmp_df6["modality_name"] = []
        mod_names_id = {}
        for m in range(len(mod_names)):
            mod_names_id[m] = mod_names[m]
        for i in range(nrows):
            acq_par_comb_id = int(data_frame.iloc[i]["modality"])
            tmp_df6["modality_name"].append(mod_names_id[acq_par_comb_id])
        tmp6 = pd.DataFrame(tmp_df6)
        df_combined = df_combined.join(tmp6)




    return data_frame, df_combined


def create_param_combinations(**kwargs):
    """
    Generate all combinations of parameter values.

    Parameters
    ----------
    **kwargs
        Parameter names and their possible values (lists).

    Returns
    -------
    combinations_dict : dict
        Dictionary mapping combination index to parameter dictionary.
    """
    if kwargs:
        # Generate all combinations using itertools.product
        combinations = list(itertools.product(*kwargs.values()))
        # Create a new dictionary with unique integers as keys
        # Each value in the dictionary will be another dictionary where the keys are the parameter names
        combinations_dict = {
            i: {key: value for key, value in zip(kwargs.keys(), combination)}
            for i, combination in enumerate(combinations)
        }
        return combinations_dict


def pivot_dataframes_byCategory(dataframe, category_name, param1, param2, metric_name, **kwargs):
    """
    Pivot a DataFrame by category and two parameters, summarizing a metric.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        DataFrame to pivot.
    category_name : str
        Name of the category column.
    param1 : str
        First parameter for pivot index.
    param2 : str
        Second parameter for pivot columns.
    metric_name : str
        Name of the metric to summarize.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    df_categories : dict
        Dictionary of pivoted DataFrames for each category.
    titles : dict
        Dictionary with category, param1, and param2 names.
    """
    df_categories = dict()
    for category, group in dataframe.groupby(category_name):
        summarised_group = None
        summarised_group = (
            group.groupby([param1, param2])
            .agg(
                Mean_Value=(metric_name, "mean"),
                Std_Dev=(metric_name, "std"),
                Replicas_Count=("Replica", "count"),
            )
            .reset_index()
        )
        # get mean and std accross parameter combinations of axes_param_names
        condition_mean_pivot = summarised_group.pivot(
            index=param1, columns=param2, values="Mean_Value"
        ).round(4)
        condition_sd_pivot = summarised_group.pivot(
            index=param1, columns=param2, values="Std_Dev"
        ).round(4)
        df_categories[str(category)] = [condition_mean_pivot, condition_sd_pivot]
    titles = dict(category=category_name, param1=param1, param2=param2)
    return df_categories, titles


def probe_parameters_sweep(
        probe_target_type: list[str] = None,
        probe_target_value: list[str] = None,
        probe_distance_to_epitope: float = None,
        probe_model: list[str] = None,
        probe_fluorophore: str = None,
        probe_paratope: str = None,
        probe_conjugation_target_info = None,
        probe_conjugation_efficiency: list[float] = None,
        probe_seconday_epitope = None,
        probe_wobbling = None,
        labelling_efficiency: list[float] = None,
    ):
    """
    Generate combinations of probe parameters for a sweep.

    Parameters
    ----------
    probe_target_type : list of str, optional
    probe_target_value : list of str, optional
    probe_distance_to_epitope : float, optional
    probe_model : list of str, optional
    probe_fluorophore : str, optional
    probe_paratope : str, optional
    probe_conjugation_target_info : any, optional
    probe_conjugation_efficiency : list of float, optional
    probe_seconday_epitope : any, optional
    probe_wobbling : any, optional
    labelling_efficiency : list of float, optional

    Returns
    -------
    probe_parameters : dict or None
        Dictionary of parameter combinations, or None if no parameters.
    """
    local_params = locals()
    probe_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                probe_parameters_vectors[par] = value
            else:
                if isinstance(value[0], (str, bool)):
                    probe_parameters_vectors[par] = value
                else:
                    sequence = np.linspace(value[0],value[1],value[2])
                    probe_parameters_vectors[par] = sequence
    probe_parameters = None
    if bool(probe_parameters_vectors):
        probe_parameters = create_param_combinations(**probe_parameters_vectors)
    return probe_parameters

def virtual_sample_parameters_sweep(
        virtual_sample_template: str = None,
        sample_dimensions: list[float] = None,
        number_of_particles: int = None,
        particle_positions: list[np.array] = None,
        random_orientations = False,
        random_placing = False,
    ):
    """
    Generate combinations of virtual sample parameters for a sweep.

    Parameters
    ----------
    virtual_sample_template : str, optional
    sample_dimensions : list of float, optional
    number_of_particles : int, optional
    particle_positions : list of numpy.array, optional
    random_orientations : bool, optional
    random_placing : bool, optional

    Returns
    -------
    field_parameters : dict or None
        Dictionary of parameter combinations, or None if no parameters.
    """
    local_params = locals()
    field_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                field_parameters_vectors[par] = value
            else:
                if isinstance(value[0], (str, bool)):
                    field_parameters_vectors[par] = value
                else:
                    sequence = np.linspace(value[0],value[1],value[2])
                    field_parameters_vectors[par] = sequence
    field_parameters = None
    if bool(field_parameters_vectors):
        field_parameters = create_param_combinations(**field_parameters_vectors)
    return field_parameters

def defects_parameters_sweep(
        defect_small_cluster: float = None,
        defect_large_cluster: float = None,
        defect: float = None,
    ):
    """
    Generate combinations of defect parameters for a sweep.

    Parameters
    ----------
    defect_small_cluster : float, optional
    defect_large_cluster : float, optional
    defect : float, optional

    Returns
    -------
    defects_parameters : dict or None
        Dictionary of parameter combinations, or None if no parameters.
    """
    local_params = locals()
    defects_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                defects_parameters_vectors[par] = value
            else:
                if isinstance(value[0], (str, bool)):
                    defects_parameters_vectors[par] = value
                else:
                    sequence = np.linspace(value[0],value[1],value[2])
                    defects_parameters_vectors[par] = sequence
    defects_parameters = None
    if bool(defects_parameters_vectors):
        defects_parameters = create_param_combinations(**defects_parameters_vectors)
    return defects_parameters

def modality_parameters_sweep(
        #modality_name: str = None,
        pixelsize_nm: float = None,
        lateral_resolution_nm: float = None,
        axial_resolution_nm: float = None,
        psf_voxel_nm: int = None,
    ):
    """
    Generate combinations of modality parameters for a sweep.

    Parameters
    ----------
    pixelsize_nm : float, optional
    lateral_resolution_nm : float, optional
    axial_resolution_nm : float, optional
    psf_voxel_nm : int, optional

    Returns
    -------
    modality_parameters : dict or None
        Dictionary of parameter combinations, or None if no parameters.
    """
    local_params = locals()
    modality_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                modality_parameters_vectors[par] = value
            else:
                if isinstance(value[0], (str, bool)):
                    modality_parameters_vectors[par] = value
                else:
                    sequence = np.linspace(value[0],value[1],value[2])
                    modality_parameters_vectors[par] = sequence
    modality_parameters = None
    if bool(modality_parameters_vectors):
        modality_parameters = create_param_combinations(**modality_parameters_vectors)
    return modality_parameters


def acquisition_parameters_sweep(
        exp_time: str = None,
        noise: float = None,
        nframes: float = None,
        channels: float = None,
    ):
    """
    Generate combinations of acquisition parameters for a sweep.

    Parameters
    ----------
    exp_time : str, optional
    noise : float, optional
    nframes : float, optional
    channels : float, optional

    Returns
    -------
    acq_parameters : dict or None
        Dictionary of parameter combinations, or None if no parameters.
    """
    local_params = locals()
    acq_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                acq_parameters_vectors[par] = value
            else:
                if isinstance(value[0], (str, bool)):
                    acq_parameters_vectors[par] = value
                else:
                    sequence = np.linspace(value[0],value[1],value[2])
                    acq_parameters_vectors[par] = sequence
    acq_parameters = None
    if bool(acq_parameters_vectors):
        acq_parameters = create_param_combinations(**acq_parameters_vectors)
    return acq_parameters