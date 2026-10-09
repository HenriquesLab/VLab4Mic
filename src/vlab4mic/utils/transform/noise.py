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


def stochastic_round(imagestack):
    """
    Round expected photon counts to integers without bias.

    Each value x becomes floor(x) + 1 with probability x - floor(x), and
    floor(x) otherwise, so the mean is kept (flooring would remove the
    signal of dim pixels, such as PSF tails). Negative values become 0.
    """
    values = np.clip(np.asarray(imagestack, dtype=float), 0, None)
    low = np.floor(values)
    return (low + (np.random.random(values.shape) < (values - low))).astype(np.int64)


def add_binomial_noise(imagestack, p=1.0):
    """
    Photon detection with quantum efficiency p.

    Expected photon counts are rounded without bias (stochastic_round) and
    each photon is detected with probability p.
    """
    imstack_binom = np.random.binomial(stochastic_round(imagestack), p)
    return imstack_binom


def add_gamma_noise(imagestack, g=1.0, em_gain=False):
    """
    Signal amplification with gain g.

    With em_gain=True, the amplification is gamma distributed (shape =
    number of electrons, scale = g), as in an EMCCD, which doubles the
    variance (excess noise factor 2). Otherwise the gain is deterministic,
    as in an sCMOS or CCD camera.
    """
    if em_gain:
        return np.random.gamma(imagestack, scale=g)
    return np.asarray(imagestack, dtype=float) * g


def add_conversion_factor(imagestack, adu: int = 1):
    """
    Emulates the digitalisation process for generating
    pixel values as integers
    """
    # print(f"Using conversion factor: ADU = {adu}")
    adu_stack = np.floor(imagestack / adu).astype(np.int32)
    return adu_stack


def add_integer_baselevel(imagestack, bl: int = 0):
    """
    Adds a constant value as integer.
    As the input data should be integers as well it is
    explicitly casted as integer
    """
    stack_with_bl = (imagestack.astype(np.int32)) + bl
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
