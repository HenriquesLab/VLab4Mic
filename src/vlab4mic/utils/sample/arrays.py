import numpy as np


def sample_spherical_normalised(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.reshape(-1)


def boolean_epitope_selection(epitopes, selected_epitopes, new_epitope, min_distance):
    """
    Check whether an epitope is free of steric hindrance.

    Parameters
    ----------
    epitopes : numpy.ndarray
        Nx3 array of epitope coordinates.
    selected_epitopes : list of int
        Indices (rows of epitopes) of the epitopes already bound.
    new_epitope : int
        Index of the epitope to check.
    min_distance : float
        Minimum distance between bound epitopes, in the units of epitopes.
        0 disables the check.

    Returns
    -------
    int
        1 if new_epitope is at least min_distance away from every selected
        epitope (or if min_distance is 0 or nothing is selected yet),
        0 otherwise.
    """
    if min_distance == 0 or len(selected_epitopes) == 0:
        return 1
    distances = np.linalg.norm(
        epitopes[selected_epitopes] - epitopes[new_epitope], axis=1
    )
    return int(np.all(distances >= min_distance))


def sample_epitopes_sterically(epitopes, min_distance, p=1):
    """
    Select the epitopes bound by a probe.

    Epitopes are visited in random order. Each one is bound with probability p
    (labelling efficiency), unless it lies closer than min_distance to an
    epitope that is already bound (steric hindrance). Epitopes that fail the
    efficiency trial stay unbound and do not block their neighbours.

    Random draws use the global numpy random state, so results are
    reproducible after np.random.seed.

    Parameters
    ----------
    epitopes : numpy.ndarray
        Nx3 array of epitope coordinates.
    min_distance : float
        Minimum distance between bound epitopes, in the units of epitopes.
        0 or None disables steric hindrance.
    p : float, optional
        Probability that a probe binds each epitope (labelling efficiency),
        between 0 and 1. Default is 1.

    Returns
    -------
    list of int
        Indices of the bound epitopes, in the random order they were bound.
    """
    n_epitopes = epitopes.shape[0]
    randomized_indices = np.random.permutation(n_epitopes)
    # bernoulli trials for probe to bind each epitope
    binding_trials = np.random.rand(n_epitopes) < p
    candidates = randomized_indices[binding_trials]
    if not min_distance:
        return candidates.tolist()
    selected_indices = []
    for i in candidates:
        # i here contains the index of the next epitope
        # verify is the next epitope is within the minimum distance
        if boolean_epitope_selection(epitopes, selected_indices, i, min_distance):
            selected_indices.append(i)
    return selected_indices


def sample_array(array, fraction=1):
    # assumes array has shape Nx3, so it will sample rows
    if fraction < 1:
        size = array.shape[0]
        ids = np.random.choice(np.arange(0, size), int(size * fraction), replace=False)
        subset = array[ids, :]
        return subset
    elif fraction == 1:
        return array


def binomial_epitope_sampling(epitopes, p=1, normals=None, min_distance=0.0):
    """
    Sample the epitopes bound by a probe, with their normals.

    Applies the labelling efficiency and steric hindrance model of
    sample_epitopes_sterically and returns the bound epitopes.

    Parameters
    ----------
    epitopes : numpy.ndarray
        Nx3 array of epitope coordinates.
    p : float, optional
        Probability that a probe binds each epitope (labelling efficiency),
        between 0 and 1. Default is 1.
    normals : numpy.ndarray, optional
        Nx3 array of normals, one per epitope. If given, the normals of the
        bound epitopes are returned too.
    min_distance : float, optional
        Minimum distance between bound epitopes, in the units of epitopes.
        Default is 0 (no steric hindrance).

    Returns
    -------
    subset_epitopes : numpy.ndarray
        Mx3 array with the coordinates of the M bound epitopes.
    n_epitopes : int
        Number of bound epitopes, M.
    subset_normals : numpy.ndarray or None
        Mx3 array with the normals of the bound epitopes, or None if
        normals is None.
    """
    ids_selected = sample_epitopes_sterically(
        epitopes=epitopes,
        min_distance=min_distance,
        p=p)
    n_epitopes = len(ids_selected)
    # create a subset of only the epitopes that were selected
    subset_epitopes = epitopes[ids_selected, :]
    if normals is None:
        return subset_epitopes, n_epitopes,  None
    else:
        subset_normals = normals[ids_selected, :]
        return subset_epitopes, n_epitopes, subset_normals

def get_random_pixels(binary_image, num_pixels=3, min_distance=2):
    # Find indices of all positive pixels (value = 1)
    positive_pixels = np.argwhere(binary_image > 0)
    selected_pixels = []
    np.random.shuffle(positive_pixels)
    for pixel in positive_pixels:
        selected = True
        for i in range(len(selected_pixels)):
          if np.linalg.norm(pixel - selected_pixels[i]) < min_distance:
             selected = False
             break
        if selected:
            selected_pixels.append(pixel)
        if len(selected_pixels) >=  num_pixels:
           break
    return selected_pixels