from ..io.yaml_functions import load_yaml
import os

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
    mode_file = modality_id + ".yaml"
    modality_config_path = os.path.join(congifuration_directory, "modalities", mode_file)
    mod_pars = load_yaml(modality_config_path)
    # approximate psf stack size
    psfx_pixels = int(mod_pars["detector"]["FOV_size"]["X"] / mod_pars["PSF"]["resolution"]["X"])
    psfy_pixels = int(mod_pars["detector"]["FOV_size"]["Y"] / mod_pars["PSF"]["resolution"]["Y"])
    if mod_pars["depth"] < mod_pars["PSF"]["resolution"]["Z"]:
        psfz_pixels = int(mod_pars["depth"] / mod_pars["PSF"]["resolution"]["Z"])
    else:
        psfz_pixels = int(mod_pars["PSF"]["resolution"]["Z"])
    psf_params = dict(
        stack_source = mod_pars["PSF"]["source"],
        scale = mod_pars["scale"],
        shape=[  # arbitrary 
            2*psfx_pixels,
            2*psfy_pixels,
            4*psfz_pixels
        ],
        std_devs=[
            mod_pars["PSF"]["resolution"]["X"] / mod_pars["PSF"]["voxelsize"],
            mod_pars["PSF"]["resolution"]["Y"] / mod_pars["PSF"]["voxelsize"],
            mod_pars["PSF"]["resolution"]["Z"] / mod_pars["PSF"]["voxelsize"]
            ],
        voxelsize=[
            mod_pars["PSF"]["voxelsize"],
            mod_pars["PSF"]["voxelsize"],
            mod_pars["PSF"]["voxelsize"]
        ],
        depth=int(mod_pars["depth"]/mod_pars["PSF"]["voxelsize"])
    )
    factor_ = 1000  # to match detector params to imager, needs fix
    detector_params = dict(
        scale=1.0e-06,
        pixelsize = mod_pars["detector"]["pixelsize"]/factor_,
        noise_model=dict(
                binomial= {"p":mod_pars["detector"]["noise"]["binomial"]},
                gamma={"g":1},
                baselevel={"bl":0},
                gaussian={"sigma":0},
                conversion={"adu":1}
        ),
        noise_order= [
            "binomial",
            "gamma",
            "baselevel",
            "gaussian",
            "conversion"
        ],
        bits_pixel=32
    )
    # get filters default should be one per channel
    if fluo_emissions is None:
        filter_dictionary = None
        emission_behaviour = None
        modality_name = mod_pars["name"]
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
            "filters" not in mod_pars.keys()
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
            filter_dictionary = mod_pars["filters"]
        # get emission type from fluorophore
        emission_behaviour = define_emission_behaviour("constant", fluo_emissions)
        print(f"emission for {str(modality_id)} set to {emission_behaviour}")
        # assign modality name
        modality_name = mod_pars["name"]
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
    save=False,
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
