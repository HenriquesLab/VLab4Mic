from ..io.yaml_functions import load_yaml


def compile_modality_parameters(
    modality_id: list, congifuration_directory, fluo_emissions=None
):
    """
    args:
        modality_id: (string) ID of modality
        congifuration_directory: (string) Path of configuration settings
        emission: (string) Emission type, either constant, or blinking.

    returns:
        modality_optical_pars: dictionary of parameters per modality
    """
    # load modality configuration file
    modality_config_path = (
        congifuration_directory + "/modalities/" + modality_id + ".yaml"
    )
    modality_params = load_yaml(modality_config_path)
    mod_emission = modality_params["emission_type"]
    # print(modality_params)
    # get PSF params
    psf_id = modality_params["psf"]
    psf_config_path = congifuration_directory + "/psfs/" + psf_id + ".yaml"
    psf_params = load_yaml(psf_config_path)
    # get detector params
    detector_id = modality_params["detector"]
    detector_config_path = (
        congifuration_directory + "/detector/" + detector_id + ".yaml"
    )
    detector_params = load_yaml(detector_config_path)
    # get filters default should be one per channel
    if fluo_emissions is None:
        filter_dictionary = None
        emission_behaviour = None
        modality_name = modality_params["id"]
        modality_params = dict(
            filters=filter_dictionary,
            psf_params=psf_params,
            detector=detector_params,
            emission=emission_behaviour,  # blinking or constant
            modality=modality_name,
        )
        return modality_params
    else:
        if (
            "filters" not in modality_params.keys()
        ):  # each fluorophore will be assign and independent channel
            print("Creating channel for each fluorophore")
            ch = 0
            filter_dictionary = dict()
            for fluoname in list(fluo_emissions.keys()):
                print(fluoname)
                channel_name = "ch" + str(ch)
                fluorophores_in_chanel = []
                fluorophores_in_chanel.append(
                    fluoname
                )  # this is because a channel can have several fluorophores
                filter_dictionary[channel_name] = fluorophores_in_chanel
                ch += 1
        else:
            filter_dictionary = modality_params["filters"]
        # get emission type from fluorophore
        emission_behaviour = define_emission_behaviour(mod_emission, fluo_emissions)
        print(f"emission for {str(modality_id)} set to {emission_behaviour}")
        # assign modality name
        modality_name = modality_params["id"]
        modality_params = dict(
            filters=filter_dictionary,
            psf_params=psf_params,
            detector=detector_params,
            emission=emission_behaviour,  # blinking or constant
            modality=modality_name,
        )
        return modality_params


def define_emission_behaviour(mod_emission, fluo_emissions):
    """
    Function that checks if the modality emission and fluorophores emission match
    """
    emission = None
    if mod_emission == "constant":
        emission = "constant"
    elif mod_emission == "blinking":
        emission = "blinking"
        for key, value in fluo_emissions.items():
            if value == "constant":
                emission = "constant"

    return emission


def format_modality_acquisition_params(
    exp_time=0.001,
    noise=True,
    save=True,
    nframes=10,
    channels=[
        "ch0",
    ],
    **kwargs,
):
    """
    Create expected format for acquisition parameters
    """
    mod_acquisition_params = dict(
        exp_time=exp_time, noise=noise, save=save, nframes=nframes, channels=channels
    )
    return mod_acquisition_params
