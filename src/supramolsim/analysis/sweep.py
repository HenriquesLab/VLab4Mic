from ..experiments import ExperimentParametrisation
from ..utils.data_format import configuration_format


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
        probe_n = 0
        for rep in range(repetitions):
            for probe in probes:
                combination = ""
                combination += str(probe_n)
                # _probe
                print(f"probe: {probe}")
                probe_param_n = 0
                for p_param in probe_parameters:
                    combination += str(probe_param_n)
                    experiment.structure_label = None
                    experiment.particle = None
                    # particle defects should be here
                    # print(p_param)
                    experiment.structure_label = probe
                    if p_param is not None:
                        if "fluorophpre_id" not in probe_parameters.keys():
                            experiment.fluorophore_id = "AF647"
                    # experiment.label_parameters[probe] = p_param
                    # experiment._build_label()
                    experiment._build_particle(keep=True)
                    if len(experiment.particle.emitters) == 0:
                        print(f"Skipping {probe}. No emitters were generated")
                        break
                    vsample_n = 0
                    for vsample in virtual_samples:
                        combination += str(vsample_n)
                        print(combination)
                        _exported_field = experiment._build_coordinate_field(
                            keep=False, use_self_particle=True
                        )
                        # print(_exported_field)
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
                        iteration_name = combination
                        iteration_output = experiment.run_simulation(
                            name=iteration_name, save=False
                        )
                        mod_n = 0
                        for mod in modalities.keys():
                            mod_comb = combination + str(mod_n)
                            _parameters = [struct, probe, p_param, vsample, mod, modalities[mod]]
                            if mod_comb not in sweep_params.keys():
                                sweep_params[mod_comb] = _parameters
                                sweep_outputs[mod_comb] = []
                            sweep_outputs[mod_comb].append(iteration_output[mod])
                            mod_n += 1
                    vsample_n += 1
                probe_param_n += 1
            probe_n += 1
        return sweep_outputs, sweep_params
        
