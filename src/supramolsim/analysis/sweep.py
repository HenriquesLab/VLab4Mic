from ..experiments import ExperimentParametrisation
from ..utils.data_format import configuration_format
from ..utils.io.yaml_functions import load_yaml
import os
import supramolsim

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
    repetitions = 1,
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
                #combination += str(probe_n)
                experiment.remove_probes()
                print(f"probe: {probe}")
                probe_param_n = 0
                for p_param in probe_parameters:
                    #combination += str(probe_param_n)
                    
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
                        #combination += str(vsample_n)
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
                        #iteration_name = combination
                        iteration_output = experiment.run_simulation(
                            name="", save=False
                        )
                        mod_n = 0
                        combination_n = str(probe_n) + str(probe_param_n) + str(vsample_n) 
                        for mod in modalities.keys():
                            mod_comb = combination_n + str(mod_n)
                            _parameters = [struct, probe, p_param, vsample, mod, modalities[mod]]
                            if mod_comb not in sweep_params.keys():
                                sweep_params[mod_comb] = _parameters
                                sweep_outputs[mod_comb] = []
                            sweep_outputs[mod_comb].append(iteration_output[mod])
                            mod_n += 1
                        vsample_n += 1
                    probe_param_n += 1
                probe_n += 1 # changes more due to repetitions
        return sweep_outputs, sweep_params
        


def sweep_vasmples(
    experiment: ExperimentParametrisation = None,
    structures=None,
    probes=None,
    probe_parameters=None,
    virtual_samples=None,
    modalities=None,
    repetitions = 1,
    **kwargs,
):
    # empty lists, but fill up with default
    default_probe = "NHS_ester"
    default_fluorophore = "AF647"
    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    local_dir = os.path.join(pck_dir, "configs")
    if experiment is None:
        experiment = ExperimentParametrisation()
    if structures is None:
        structures = ["1XI5",]
    if probes is None:
        probes = [default_probe,]
    if probe_parameters is None:
        filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        default_params = load_yaml(filepath) 
        probe_parameters = [default_params, ]
    if virtual_samples is None:
        filepath = os.path.join(local_dir, "probes", default_probe + ".yaml")
        #default_vsample =  load_yaml(filepath) 
        virtual_samples = [None,]

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
                    p_param["fluorophore_id"] = default_fluorophore
                    experiment.remove_probes()
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
                        #combination += str(vsample_n)
                        _exported_field = experiment._build_coordinate_field(
                            keep=False, use_self_particle=True
                        )
                        combination_n = str(probe_n) + str(probe_param_n) + str(vsample_n)
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
        vsample_outputs = None
    ):
    if vsample_outputs is None:
        pass # default sample can be a minimal field with a single emitter