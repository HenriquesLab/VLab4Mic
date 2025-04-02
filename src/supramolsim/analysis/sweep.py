from ..experiments import ExperimentParametrisation
from ..utils.data_format import configuration_format
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
        probes = [
            default_probe,
        ]
    if probe_parameters is None:
        probe_filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        default_params = load_yaml(probe_filepath)
        probe_parameters = dict()
        probe_parameters[0] = default_params
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
        for rep in range(repetitions):
            probe_n = 0
            for probe in probes:
                print(f"probe: {probe}")
                # probe_param_n = 0
                for probe_param_n, p_param in probe_parameters.items():
                    # copy p_param?
                    experiment.remove_probes()
                    if p_param is None:
                        experiment.probe_parameters[probe] = default_params
                    else:
                        p_param_copy = copy.deepcopy(p_param)
                        p_param_copy["fluorophore_id"] = default_fluorophore
                        experiment.structure_label = probe
                        experiment.probe_parameters[probe] = p_param_copy
                    for defect_n, defects_pars in particle_defects.items():
                        particle_defects_copy = copy.deepcopy(defects_pars)
                        for d_key, d_val in particle_defects_copy.items():
                            experiment.defect_eps[d_key] = d_val
                        # print(experiment.defect_eps)
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
                            _exported_field = experiment._build_coordinate_field(
                                keep=False, use_self_particle=True, **vsample_pars
                            )
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
                                particle_defects,
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


def sweep_modalities(
    experiment: ExperimentParametrisation = None,
    vsample_outputs=None,
    vsampl_pars=None,
    modalities: list = None,
    modality_params=None
):
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
    experiment._build_imager(use_local_field=False)
    print(experiment.objects_created["imager"])
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
            with io.capture_output() as captured:
                experiment.imager.import_field(**virtualsample)
                # iteration_name = combination
                iteration_output = experiment.run_simulation(name="", save=False)
                mod_outputs
                mod_n = 0
                for mod, acq in experiment.selected_mods.items():
                    print(f"modality and acq: {mod}, {acq}")
                    mod_comb = vsmpl_id + "_" + str(mod_n)
                    mod_parameters = copy.copy(vsampl_pars[vsmpl_id])
                    mod_parameters.append(mod)
                    mod_parameters.append(acq)
                    # mod_parameters = [vsampl_pars[vsmpl_id], mod, acq]
                    if mod_comb not in mod_params.keys():
                        mod_params[mod_comb] = mod_parameters
                        mod_outputs[mod_comb] = []
                    mod_outputs[mod_comb].append(iteration_output[mod])
                    mod_n += 1
                    mod_parameters = None
    return experiment, mod_outputs, mod_params, pixelsizes

def sweep_modalities_updatemod(
    experiment: ExperimentParametrisation = None,
    vsample_outputs=None,
    vsampl_pars=None,
    modalities: list = None,
    modality_params: dict = None,
    modality_acq_prams: dict =None,
):
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
        modality_acq_prams[0] = {"modality_name": "NA"}
    experiment._build_imager(use_local_field=False)
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
                            experiment.set_modality_acq(
                                modality_name=modality_name,
                                **acq_pars)
                            # iteration_name = combination
                            modality_timeseries = experiment.run_simulation(name="", save=False, modality=modality_name)
                            mod_comb = vsmpl_id + "_" + str(mod_n) + "_" + str(mod_pars_number) + "_" + str(mod_acq_number)
                            mod_parameters = copy.copy(vsampl_pars[vsmpl_id])
                            mod_parameters.append(modality_name)
                            mod_parameters.append(mod_pars)
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
    experiment.structure_label = probe
    experiment.probe_parameters[probe] = probe_parameters
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
):
    if experiment is None:
        experiment = ExperimentParametrisation()
    if modality is None:
        modality = "Reference"
        modality_acquisition = None
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
        experiment.imager.modalities["Reference"]["detector"]["pixelsize"] * scalefactor
    )
    return reference_output[modality], reference_parameters


def analyse_image_sweep(img_outputs, img_params, reference, analysis_case_params=None):
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


def measurements_dataframe(
    measurement_vectors, probe_parameters=None, p_defects=None, mod_names = None, sample_params=None, mod_params=None, mod_acq=None
):
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
            "Replica": measurement_array[:, 1],
            "Metric": np.array(measurement_array[:, 2], dtype=np.float32),
        }
    )
    df_combined = data_frame
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


def pivot_dataframe(dataframe, param1, param2):
    # extract individual dataframe per condition
    summarised_df = None
    summarised_df = (
        dataframe.groupby([param1, param2])
        .agg(
            Mean_Value=("Metric", "mean"),
            Std_Dev=("Metric", "std"),
            Replicas_Count=("Replica", "count"),
        )
        .reset_index()
    )
    # get mean and std accross parameter combinations of axes_param_names
    condition_mean_pivot = summarised_df.pivot(
        index=param1, columns=param2, values="Mean_Value"
    ).round(4)
    condition_sd_pivot = summarised_df.pivot(
        index=param1, columns=param2, values="Std_Dev"
    ).round(4)
    return condition_mean_pivot, condition_sd_pivot


def pivot_dataframes_byCategory(dataframe, category_name, param1, param2, **kwargs):
    df_categories = dict()
    for category, group in dataframe.groupby(category_name):
        summarised_group = None
        summarised_group = (
            group.groupby([param1, param2])
            .agg(
                Mean_Value=("Metric", "mean"),
                Std_Dev=("Metric", "std"),
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
    local_params = locals()
    probe_parameters_vectors = {}
    for par, value in local_params.items():
        if value is not None and type(value) is list:
            if len(value) == 1:
                probe_parameters_vectors[par] = value
            else:
                if isinstance(value[0], str):
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