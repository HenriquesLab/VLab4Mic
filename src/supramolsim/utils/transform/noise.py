import numpy as np


def add_poisson_noise(imagestack):
    # print("adding poisson noise")
    noisy = np.random.poisson(imagestack)
    return noisy


def add_gaussian_noise(imagestack, sigma):
    # print(f"adding gaussian noise: sigma = {sigma}")
    gnoise = np.random.normal(0, sigma, size=imagestack.shape)
    noisy = imagestack + gnoise
    return noisy


def add_binomial_noise(imagestack, p=1):
    # print(f"adding binomial noise: p = {p}")
    int_imagestack = np.floor(imagestack).astype("int64")
    return np.random.binomial(int_imagestack, p)


def add_gamma_noise(imagestack, g=1.0):
    # print(f"adding gamma noise: g = {g}")
    amplified = np.random.gamma(imagestack, scale=g)
    return amplified


def add_conversion_factor(imagestack, adu=1):
    """
    Emulates the digitalisation process for generating
    pixel values as integers
    """
    # print(f"Using conversion factor: ADU = {adu}")
    adu_stack = np.floor(imagestack / adu).astype("int64")
    return adu_stack


def add_integer_baselevel(imagestack, bl=0):
    """
    Adds a constant value as integer.
    As the input data should be integers as well it is
    explicitly casted as integer
    """
    stack_with_bl = (imagestack.astype("int64")) + bl
    return stack_with_bl


def add_image_noise(noise_type: str, stack, **kwargs):
    """
    Add noise on a image-based scheme. This models a EMCCD
    or any noise that is not pixel-dependent
    """
    if noise_type == "binomial":
        noisy = add_binomial_noise(stack, **kwargs)
    elif noise_type == "gaussian":
        noisy = add_gaussian_noise(stack, **kwargs)
    elif noise_type == "gamma":
        noisy = add_gamma_noise(stack, **kwargs)
    elif noise_type == "poisson":
        noisy = add_poisson_noise(stack)
    elif noise_type == "conversion":
        noisy = add_conversion_factor(stack, **kwargs)
    elif noise_type == "baselevel":
        noisy = add_integer_baselevel(stack, **kwargs)
    else:
        noisy = None
    return noisy
