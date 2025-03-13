from ..experiments import ExperimentParametrisation
from ..utils.data_format import configuration_format
from ..utils.io.yaml_functions import load_yaml
import os
import supramolsim
from . import metrics
import numpy as np

def nested_sweep(
    experiment=ExperimentParametrisation,
    structures=[
        None,
    ],
    probes=[
        None,
    ],
    probe_parameters=[
        None,
    ],
    virtual_samples=[
        None,
    ],
    modalities=[
        None,
    ],
    repetitions=1,
    **kwargs,
    # empty lists, but fill up with default
):
    """
    Series of nested loops to cover every parameter combination
    in parameters.

    Input:
    """
    sweep_params = dict()
    sweep_outputs = dict()
    for struct in structures:
        experiment.structure_id = struct
        experiment._build_structure()
        # probe_parameters = probes["target_info"]
        for rep in range(repetitions):
            probe_n = 0

            for probe in probes:
                combination = ""
                # combination += str(probe_n)
                experiment.remove_probes()
                print(f"probe: {probe}")
                probe_param_n = 0
                for p_param in probe_parameters:
                    # combination += str(probe_param_n)

                    experiment.structure_label = probe
                    if p_param is not None:
                        if "fluorophpre_id" not in p_param.keys():
                            experiment.fluorophore_id = "AF647"
                        experiment.probe_parameters[probe] = p_param
                    experiment._build_particle(keep=True)
                    if experiment.generators_status("particle"):
                        if len(experiment.particle.emitters) == 0:
                            print(f"Skipping {probe}. No emitters were generated")
                            break
                    else:
                        break
                    vsample_n = 0
                    for vsample in virtual_samples:
                        _exported_field = None
                        # combination += str(vsample_n)
                        _exported_field = experiment._build_coordinate_field(
                            keep=False, use_self_particle=True
                        )
                        for (
                            modality,
                            acquisition,
                        ) in modalities.items():  # load all modalities
                            if acquisition is None:
                                experiment.selected_mods[modality] = (
                                    configuration_format.format_modality_acquisition_params()
                                )
                            else:
                                experiment.selected_mods[modality] = acquisition
                        experiment._build_imager(use_local_field=False)
                        experiment.imager.import_field(**_exported_field)
                        # iteration_name = combination
                        iteration_output = experiment.run_simulation(
                            name="", save=False
                        )
                        mod_n = 0
                        combination_n = (
                            str(probe_n) + str(probe_param_n) + str(vsample_n)
                        )
                        for mod in modalities.keys():
                            mod_comb = combination_n + str(mod_n)
                            _parameters = [
                                struct,
                                probe,
                                p_param,
                                vsample,
                                mod,
                                modalities[mod],
                            ]
                            if mod_comb not in sweep_params.keys():
                                sweep_params[mod_comb] = _parameters
                                sweep_outputs[mod_comb] = []
                            sweep_outputs[mod_comb].append(iteration_output[mod])
                            mod_n += 1
                        vsample_n += 1
                    probe_param_n += 1
                probe_n += 1  # changes more due to repetitions
        return sweep_outputs, sweep_params


def sweep_vasmples(
    experiment: ExperimentParametrisation = None,
    structures=None,
    probes=None,
    probe_parameters=None,
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
        probe_parameters = [
            default_params,
        ]
    if virtual_samples is None:
        # vprobe_filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        # default_vsample =  load_yaml(vprobe_filepath)
        virtual_samples = [
            None,
        ]
    vsample_params = dict()
    vsample_outputs = dict()

    for struct in structures:
        experiment.structure_id = struct
        experiment._build_structure()
        for rep in range(repetitions):
            probe_n = 0
            for probe in probes:
                print(f"probe: {probe}")
                probe_param_n = 0
                for p_param in probe_parameters:
                    experiment.remove_probes()
                    if p_param is None:
                        experiment.probe_parameters[probe] = default_params
                    else:
                        p_param["fluorophore_id"] = default_fluorophore
                        experiment.structure_label = probe
                        experiment.probe_parameters[probe] = p_param
                    experiment._build_particle(keep=True)
                    if experiment.generators_status("particle"):
                        if len(experiment.particle.emitters) == 0:
                            print(f"Skipping {probe}. No emitters were generated")
                            break
                    else:
                        break
                    vsample_n = 0
                    for vsample in virtual_samples:
                        _exported_field = None
                        # combination += str(vsample_n)
                        _exported_field = experiment._build_coordinate_field(
                            keep=False, use_self_particle=True
                        )
                        combination_n = (
                            str(probe_n) + str(probe_param_n) + str(vsample_n)
                        )
                        _parameters = [struct, probe, p_param, vsample]
                        if combination_n not in vsample_params.keys():
                            vsample_params[combination_n] = _parameters
                            vsample_outputs[combination_n] = []
                        vsample_outputs[combination_n].append(_exported_field)
                        #
                        #
                        vsample_n += 1
                    probe_param_n += 1
                probe_n += 1
    return experiment, vsample_outputs, vsample_params


def sweep_modalities(
    experiment: ExperimentParametrisation = None,
    vsample_outputs=None,
    vsampl_pars=None,
    modalities=None,
):
    default_mod = "Confocal"
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
    if modalities is None:
        modalities = dict()
        modalities[default_mod] = None
    for modality, acquisition in modalities.items():  # load all modalities
        if acquisition is None:
            experiment.selected_mods[modality] = (
                configuration_format.format_modality_acquisition_params()
            )
        else:
            experiment.selected_mods[modality] = acquisition
    experiment._build_imager(use_local_field=False)
    for vsmpl_id in vsampl_pars.keys():
        for virtualsample in vsample_outputs[vsmpl_id]:
            experiment.imager.import_field(**virtualsample)
            # iteration_name = combination
            iteration_output = experiment.run_simulation(name="", save=False)
            mod_outputs
            mod_n = 0
            for mod, acq in modalities.items():
                mod_comb = vsmpl_id + str(mod_n)
                mod_parameters = [vsampl_pars[vsmpl_id], mod, acq]
                if mod_comb not in mod_params.keys():
                    mod_params[mod_comb] = mod_parameters
                    mod_outputs[mod_comb] = []
                mod_outputs[mod_comb].append(iteration_output[mod])
                mod_n += 1
    return experiment, mod_outputs, mod_params


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
        experiment.selected_mods[modality] = (
            configuration_format.format_modality_acquisition_params()
        )
    else:
        experiment.selected_mods[modality] = modality_acquisition
    experiment._build_imager(use_local_field=False)
    experiment.imager.import_field(**reference_vsample)
    reference_parameters = dict()
    reference_parameters["Vector"] = [reference_vsample_params, modality, modality_acquisition]
    reference_output = experiment.run_simulation(name="", save=False)
    imager_scale = experiment.imager.roi_params["scale"]
    scalefactor = np.ceil(imager_scale / 1e-9)  # resulting pixel size in nanometers
    reference_parameters["ref_pixelsize"] = experiment.imager.modalities["Reference"]["detector"]["pixelsize"]*scalefactor 
    return reference_output[modality], reference_parameters


def analyse_image_sweep(img_outputs, img_params, reference, analysis_case_params=None):
    measurement_vectors = []
    #ref_pixelsize = analysis_case_params["ref_pixelsize"]
    for params_id in img_params.keys():
        rep_number = 0
        mod_name = img_params[params_id][1]
        for img_r in img_outputs[params_id]:
            item_vector = []
            im1 = img_r[0]
            im_ref = reference[0]
            rep_measurement, ref_used, qry_used = metrics.img_compare(im1, im_ref, **analysis_case_params[mod_name])
            item_vector.append(params_id)
            item_vector.append(rep_number)
            item_vector.append(rep_measurement)
            measurement_vectors.append(item_vector)
            # inputs[params_id][rep_number] = [qry_used, im1]
            rep_number += 1
    return measurement_vectors