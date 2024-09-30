import numpy as np


def sample_spherical_normalised(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.reshape(-1)


def boolean_epitope_selection(epitopes, selected_epitopes, new_epitope, min_distance):
    # epitopes is the numpy array containing epitopes coordinates
    # selected epitopes is a list of indices of the selected epitopes
    # min_distance is the minimum distance between the epitopes
    # function should return indices of the selected epitopes
    n_epi = len(selected_epitopes)
    for i in range(n_epi):
        is_available = 1
        if np.linalg.norm(epitopes[selected_epitopes[i]] - epitopes[new_epitope]) < min_distance:
                is_available = 0
                break
    return is_available


def sample_epitopes_sterically(epitopes, min_distance):
    # epitopes is the coordinates array for the epitopes
    # min_distance is the minimum eculidean distance between the epitopes
    # step 1: select an epitope randomly
    selected_indices = []
    # randomize indices to have a guaranteed stop condition
    randomized_indices = np.random.choice(range(epitopes.shape[0]),epitopes.shape[0], replace=False )
    selected_indices.append(randomized_indices[0])
    for i in randomized_indices[1:]:
        # i here contains the index of the next epitope
        # verify is the next epitope is within the minimum distance
        if boolean_epitope_selection(epitopes, selected_indices, i, min_distance):
            selected_indices.append(i)
    # function should return indices of the selected epitopes
    return selected_indices


def sample_array(array, fraction=1):
    # assumes array has shape Nx3, so it will sample rows
    if fraction < 1:
        size = array.shape[0]
        ids = np.random.choice(np.arange(0, size), int(size*fraction), replace=False)
        subset = array[ids, :]
        return subset
    elif fraction == 1:
        return array


def binomial_epitope_sampling(epitopes, p=1, normals=None):
    n_epitopes = epitopes.shape[0]
    total_selected = np.random.binomial(n_epitopes, p)
    ids_selected = np.random.choice(np.arange(0, n_epitopes), int(total_selected), replace=False)
    subset_epitopes = epitopes[ids_selected, :]
    if normals is None:
        return subset_epitopes
    else:
        subset_normals = normals[ids_selected, :]
        return subset_epitopes, subset_normals
