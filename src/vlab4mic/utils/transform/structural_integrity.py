import warnings

import numpy as np
from sklearn.cluster import DBSCAN
from collections import Counter
import scipy


def remove_items_fromlist(test_list, item):
    # using list comprehension to perform the task
    res = [i for i in test_list if i != item]
    return res


def ids2delete2(xmer_id, tree, neigh, upbound):
    # this function generates an initial proposal for the IDs to drop
    # by considering neighbors given an initial id
    """
    Propose subunits to remove around one subunit.

    Queries the neighbours of xmer_id within upbound (at most neigh) and
    removes a random number of them, always including xmer_id itself.

    Returns
    -------
    list of int
        Subunit IDs proposed for removal.
    """
    point2query = tree.data[xmer_id, :]
    # print(f"query with parameters: {xmer_id},{point2query}, {neigh}, {upbound}")
    distances, index_in_data = tree.query(
        point2query, k=neigh, distance_upper_bound=upbound
    )
    # print(f"Result of tree.query distance: {distances}")
    # print(f"Result of tree.query index: {index_in_data}")
    # print(f"type of index: {type(index_in_data)}")

    # Handle the case when k=1, tree.query returns scalars instead of arrays
    if neigh == 1:
        if np.isscalar(index_in_data):
            index_in_data = np.array([index_in_data])
        if np.isscalar(distances):
            distances = np.array([distances])

    # neighbours beyond upbound are returned with infinite distance and an
    # out-of-range index; keep only real neighbours
    index_in_data = np.asarray(index_in_data)[np.isfinite(distances)]
    availabeids = index_in_data.shape[0]
    if availabeids <= 1:
        # no neighbour within upbound: remove only this subunit
        return [xmer_id]
    todelete = np.random.choice(np.arange(1, availabeids))
    ids2remove = np.random.choice(index_in_data, todelete, replace=False)
    xmers_removed = [int(i) for i in ids2remove]
    if xmer_id not in xmers_removed:
        xmers_removed.append(xmer_id)
    return xmers_removed  # is a list


def xmer_ids_remove(xmer0_id, tree, percentageoff, neighbors, upbound):
    totalpoints = tree.data.shape[0]
    print("total emitters", totalpoints)
    total2remove = np.floor(totalpoints * percentageoff)
    # total2remove = np.random.poisson(np.floor(totalpoints*percentageoff))
    print("removing: ", total2remove)
    removed = list([xmer0_id])
    id_query = xmer0_id
    max_iter = 100000
    iterat = 0
    while len(removed) < int(total2remove):
        if iterat > max_iter:
            print("exit due to overiteration")
            break
        tmp = ids2delete2(id_query, tree, neighbors, upbound)
        if totalpoints in tmp:  # total points would be an out of bound index
            # and it is added when quering the kdTree if the nearest neighbor is
            # outside the dist_upper_bound # https://github.com/scipy/scipy/issues/3210
            # print("before removing in tmp: ", tmp)
            tmp = remove_items_fromlist(
                tmp, totalpoints
            )  # this exception can occur multiple times
            # so we make sure to eliminate all occurrences when at leas one is detected
            # print("after removing in tmp: ", tmp)
        removed.extend(tmp)
        asarray = np.array(removed)
        uniq = np.unique(asarray)
        id_query = np.random.choice(uniq)
        removed = uniq.tolist()
        # print(f" in loop: {len(removed)}, {int(total2remove)}")
        iterat = iterat + 1
    return removed  ## also a list


def notin_logical_list(V, integers_to_check):
    # for V, we want to exclude all occurrences of integers_to_check
    # and it returns a logical vector indicating "False" for indices
    # in V where an element of integers_to_check appear
    # Create a list comprehension to generate the logical vector
    # returns a logical list that identify all occurrences of
    # each index of integers_to_check in V
    return [x not in integers_to_check for x in V]


def singlecluster_verification(
    xmer_centers, xmer_ids_all, ids_todelete, max_dist, min_samples
):
    """
    Complete a removal proposal so that the remaining subunits are connected.

    The subunits left after removing ids_todelete are clustered with
    DBSCAN (eps = max_dist). If they form several clusters, every cluster
    other than the largest is also removed, so no floating fragments
    remain.

    Parameters
    ----------
    xmer_centers : numpy.ndarray
        Coordinates of the subunit centres.
    xmer_ids_all : numpy.ndarray
        IDs of all subunits.
    ids_todelete : list of int
        IDs of the subunits proposed for removal.
    max_dist : float
        Distance below which remaining subunits are connected.
    min_samples : int
        DBSCAN min_samples.

    Returns
    -------
    numpy.ndarray or None
        IDs of all subunits to remove, or None if the proposal removes
        every subunit.
    """
    xmerids_logical = notin_logical_list(xmer_ids_all, ids_todelete)
    xmer_to_remain = xmer_ids_all[xmerids_logical]
    # up to here we are generating the complement list of xmer centers
    # that we have after ids_todelete
    xmer_subset = xmer_centers[xmerids_logical,]
    # print(xmer_subset, xmer_subset.shape)
    if xmer_subset.shape[0] == 0:
        return None
    # verify single cluster on this subset
    db_xmers = DBSCAN(eps=max_dist, min_samples=min_samples).fit(xmer_subset)
    disassembled_xmers_clusters = db_xmers.labels_
    nclusters_xmers = np.unique(disassembled_xmers_clusters).shape[0]
    if nclusters_xmers > 1:
        # keep the largest connected cluster; the other ones are floating
        # fragments and are removed too
        labels, counts = np.unique(disassembled_xmers_clusters, return_counts=True)
        largest = labels[np.argmax(counts)]
        floating = xmer_to_remain[disassembled_xmers_clusters != largest]
        return np.unique(np.concatenate([np.asarray(ids_todelete), floating]))
    return np.asarray(ids_todelete)


def xmersubset_byclustering(
    epitopes_coords,
    d_cluster_params,
    fracture=-24,
    deg_dissasembly=0.5,
    xmer_neigh_distance=100,
    return_ids=False,
):
    """
    Model structural integrity by nearest neighbors of the emitters and clustering

    """
    default_true = [True] * epitopes_coords.shape[0]
    if deg_dissasembly == 0:
        if return_ids:
            return default_true
        else:
            return epitopes_coords
    if deg_dissasembly == 1:
        if return_ids:
            default_false = [False] * epitopes_coords.shape[0]
            return default_false
        else:
            return np.array([])
    clusters1 = DBSCAN(
        eps=d_cluster_params["eps1"],
        min_samples=d_cluster_params["minsamples1"],
    ).fit(epitopes_coords)
    label_p_epitope = clusters1.labels_
    # obtain a center for each xmer
    xmer_ids_all = np.unique(label_p_epitope)
    n_xmer = len(xmer_ids_all)

    # generate a center point to refer to a xmer
    sums = np.zeros((n_xmer, 3))
    for d, l in zip(epitopes_coords, label_p_epitope):
        sums[l, :] += d
    total_labels_dictionary = Counter(
        label_p_epitope
    )  # this is a dictionary with the number of the cluster as key
    # its value are the amount of points
    center_xmers = np.zeros((n_xmer, 3))
    for i in total_labels_dictionary:
        center_xmers[i, :] = sums[i, :] / total_labels_dictionary[i]
    # define point of fracture randomly or defined
    # Define a percentage of the total objects to be removed from the whole data
    percentageoff = deg_dissasembly
    # create a list with the first proposal of xmers ids to delete
    # create the kdTree for the center_xmers
    xmer_tree = scipy.spatial.cKDTree(
        center_xmers, leafsize=10, copy_data=True
    )
    # print(total_labels_dictionary.values())
    # print(f"max num of elements on clusters: {neighbors}")
    upbound = xmer_neigh_distance  # in angstroms

    # print(f"neighbors for initial breakpoint: {len(neighbors)}")
    total_number_epitopes = epitopes_coords.shape[0]
    n_epitopes_to_keep = np.floor(total_number_epitopes * (1-percentageoff))
    lower_bound = np.floor(n_epitopes_to_keep - (total_number_epitopes*0.05))
    upper_bound = np.ceil(n_epitopes_to_keep + (total_number_epitopes*0.05))
    expected_number_reached = False
    i = 0
    epitopes_ids = None
    best_ids = None
    best_error = np.inf
    while i < 50: # maximum number of trials before returning empty selection
        # sample starting point each time if no fracture was specified
        if fracture == -24:
            xmer_fracture = np.random.choice(np.arange(0, n_xmer))
        else:
            xmer_fracture = fracture
        fracture_coord = xmer_tree.data[xmer_fracture, :]
        neighbors = xmer_tree.query_ball_point(fracture_coord, upbound)
        todelete = xmer_ids_remove(
            xmer_fracture, xmer_tree, percentageoff, len(neighbors), upbound
        )
        ## verifify that the resulting subset does not contain isolated entities
        ids_validated = singlecluster_verification(
            center_xmers,
            xmer_ids_all,
            todelete,
            d_cluster_params["eps2"],
            d_cluster_params["minsamples2"],
        )
        if ids_validated is None:
            i+=1
            continue
        # ids_validated are the IDs of the subunits to remove; keep the
        # epitopes of every other subunit
        epitopes_ids = notin_logical_list(label_p_epitope, ids_validated)
        n_kept = sum(epitopes_ids)
        if abs(n_kept - n_epitopes_to_keep) < best_error:
            best_error = abs(n_kept - n_epitopes_to_keep)
            best_ids = epitopes_ids
        if n_kept >= lower_bound and n_kept <= upper_bound:
            expected_number_reached = True
            break
        i+=1
    if not expected_number_reached and best_ids is not None:
        warnings.warn(
            "Structural integrity: no removal pattern within 5% of the "
            f"requested fraction after 50 trials; using the closest one "
            f"({sum(best_ids)} of {total_number_epitopes} epitopes kept, "
            f"{int(n_epitopes_to_keep)} requested)."
        )
        epitopes_ids = best_ids
        expected_number_reached = True
    if return_ids:
        if expected_number_reached:
            return epitopes_ids
        else:
            default_false = [False] * epitopes_coords.shape[0]
            return default_false
    else:
        if expected_number_reached:
            subset = epitopes_coords[epitopes_ids,]
            return subset
        else:
            return np.array([])
