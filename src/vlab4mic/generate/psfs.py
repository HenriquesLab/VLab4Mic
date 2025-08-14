import numpy as np


def elliptical_gaussian_3sigmas(shape, std_devs):
    """
    Create a 3D elliptical Gaussian array.

    Parameters:
    - shape: Tuple, shape of the array (depth, height, width).
    - std_dev: Tuple, standard deviations along each axis (z, y, x).

    Returns:
    - ndarray: 3D NumPy array representing the elliptical Gaussian.
    """
    center = tuple(dim // 2 for dim in shape)
    z, y, x = np.ogrid[: shape[0], : shape[1], : shape[2]]
    elliptical_gaussian = np.exp(
        -(
            (z - center[0]) ** 2 / (2 * std_devs[0] ** 2)
            + (y - center[1]) ** 2 / (2 * std_devs[1] ** 2)
            + (x - center[2]) ** 2 / (2 * std_devs[2] ** 2)
        )
    )
    elliptical_gaussian = elliptical_gaussian / np.sum(elliptical_gaussian)
    return elliptical_gaussian
