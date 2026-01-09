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

    availabeids = index_in_data.shape[0]
    todelete = np.random.choice(np.arange(1, availabeids))
    ids2remove = np.random.choice(index_in_data, todelete, replace=False)
    xmers_removed = list(ids2remove)
    if xmer_id not in xmers_removed:
        xmers_removed.append(xmer_id)
    return xmers_removed  # is a list


def xmer_ids_remove(xmer0_id, tree, percentageoff, neighbors, upbound):
    totalpoints = tree.data.shape[0]
    # print("total emitters", totalpoints)
    total2remove = np.floor(totalpoints * percentageoff)
    # total2remove = np.random.poisson(np.floor(totalpoints*percentageoff))
    # print("removing: ", total2remove)
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
    xmerids_logical = notin_logical_list(xmer_ids_all, ids_todelete)
    xmer_to_remain = xmer_ids_all[xmerids_logical]
    # up to here we are generating the complement list of xmer centers
    # that we have after ids_todelete
    xmer_subset = xmer_centers[xmerids_logical,]
    # print(xmer_subset, xmer_subset.shape)
    if xmer_subset.shape[0] == 0:
        print("Verification step yield no result, try new parameters")
        return None
    # verify single cluster on this subset
    db_xmers = DBSCAN(eps=max_dist, min_samples=min_samples).fit(xmer_subset)
    disassembled_xmers_clusters = db_xmers.labels_
    nclusters_xmers = np.unique(disassembled_xmers_clusters).shape[0]
    if nclusters_xmers > 1:
        # running this line will get us the indexes in the second clustering
        # that are bigger than 0, which should be the ones "floating"
        # because by default the bigger cluster is termed 0th
        invalid_xmers_ids = [
            i for i, x in enumerate(disassembled_xmers_clusters) if x > 0
        ]
        cleanup_ids2 = xmer_to_remain[
            invalid_xmers_ids
        ]  # from the indices we wanted to preserve we remove the ones invalid
        final_logical = notin_logical_list(xmer_to_remain, cleanup_ids2)
        valid_xmer_ids = xmer_to_remain[final_logical,]
    else:
        valid_xmer_ids = ids_todelete
    return valid_xmer_ids


def xmersubset_byclustering(
    epitopes_coords,
    d_cluster_params,
    fracture=-24,
    deg_dissasembly=0.5,
    xmer_neigh_distance=100,
    return_ids=False,
):
    """
    Model Structureal defects by nearest neighbors of the emitters and clustering

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
    if fracture == -24:
        xmer_fracture = np.random.choice(np.arange(0, n_xmer))
    else:
        xmer_fracture = fracture
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

    fracture_coord = xmer_tree.data[xmer_fracture, :]
    neighbors = xmer_tree.query_ball_point(fracture_coord, upbound)
    # print(f"neighbors for initial breakpoint: {len(neighbors)}")
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
        print("error while simulating defects, returning No emitters")
        if return_ids:
            default_false = [False] * epitopes_coords.shape[0]
            return default_false
        else:
            return np.array([])
    # ids_validated are the ids of the center of each xmer that we want to preserve
    # we only need to then retrieve the appropriate indices of the epitopes themselves
    # that correspond to these labels
    epitopes_ids = notin_logical_list(
        label_p_epitope, ids_validated
    )  # this function retrieves false
    # for each time a value ids_validated appears in label_p_epitope
    final_epitopes_ids = epitopes_ids
    subset = epitopes_coords[final_epitopes_ids,]
    if return_ids:
        return final_epitopes_ids
    else:
        return subset
