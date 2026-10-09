import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import (
    breadth_first_order,
    connected_components,
    minimum_spanning_tree,
)
from scipy.spatial import cKDTree
from .points_transforms import rotate_point


def normals_by_scaling(epitope_locs, scale=0.95):
    """
    Estimate normals as the direction from the centroid of the epitopes.

    Each normal is the displacement of an epitope when the set of epitopes
    is shrunk by scale about its centroid, so it points away from the
    centroid with length (1 - scale) times the distance to it. This is the
    surface normal only for convex shapes centred on the centroid (e.g. a
    sphere); on flat or elongated surfaces use normals_by_local_plane.

    Parameters
    ----------
    epitope_locs : numpy.ndarray
        Nx3 array of epitope coordinates.
    scale : float, optional
        Shrinking factor, between 0 and 1. Default is 0.95.

    Returns
    -------
    numpy.ndarray
        Nx3 array of normals (not unit length).
    """
    scaled_epitopes_locs = coordinates_scaling(epitope_locs, scale)
    # calculate the normals of the epitopes
    # note that this approach will work for convex shapes
    # this assumes as well the scaled version is smaller
    normals = np.zeros((len(epitope_locs), 3))
    for i in range(len(epitope_locs)):
        normals[i, :] = epitope_locs[i] - scaled_epitopes_locs[i]
    # a site at the centroid (e.g. a single site) has no defined direction;
    # use +z instead of a zero vector
    zero = np.linalg.norm(normals, axis=1) == 0
    normals[zero] = [0, 0, 1]
    # then the points from which the normals are traced are the epitopes themselves
    # normals_ft_epitopes = [normals, epitope_locs]
    return normals
    # return normals_ft_epitopes

def normals_by_local_plane(
    epitope_locs, n_neighbours=10, reference_vector=None, flatness_tolerance=0.2
):
    """
    Estimate surface normals from a plane fitted around each epitope.

    For each epitope, the normal is the direction of least variance of its
    n_neighbours nearest epitopes (itself included). Unlike
    normals_by_scaling, this gives the surface normal on flat and
    non-convex surfaces as well as on convex ones.

    Normals are oriented consistently across the surface by propagating the
    side along a minimum spanning tree of neighbouring epitopes (Hoppe et
    al., 1992). In each connected patch, the side is set at the epitope whose
    normal is most aligned with the direction from the centroid of all
    epitopes, so normals point away from it. If no epitope in the patch is
    aligned better than flatness_tolerance (|cos|), as on a flat surface,
    normals point to the same side as reference_vector.

    Parameters
    ----------
    epitope_locs : numpy.ndarray
        Nx3 array of epitope coordinates.
    n_neighbours : int, optional
        Number of nearest epitopes (including itself) used to fit the plane
        at each epitope. Default is 10. Too few neighbours fail on epitopes
        that come in small clusters.
    reference_vector : numpy.ndarray, optional
        Side normals point to on flat surfaces. Default is [0, 0, 1].
    flatness_tolerance : float, optional
        Threshold on |cos| between a normal and the direction from the
        centroid, below which reference_vector decides the side. Default
        is 0.2.

    Returns
    -------
    numpy.ndarray
        Nx3 array of unit normals. With fewer than 3 epitopes, every normal
        is reference_vector (normalised).

    References
    ----------
    Hoppe, H. et al. Surface reconstruction from unorganized points.
    SIGGRAPH Comput. Graph. 26, 71-78 (1992).
    """
    locs = np.asarray(epitope_locs, dtype=float)
    n_epitopes = locs.shape[0]
    if reference_vector is None:
        reference_vector = np.array([0, 0, 1])
    reference_vector = np.asarray(reference_vector, dtype=float)
    reference_vector = reference_vector / np.linalg.norm(reference_vector)
    if n_epitopes < 3:
        # not enough points to fit a plane
        return np.tile(reference_vector, (n_epitopes, 1))
    k = min(n_neighbours, n_epitopes)
    _, neighbour_ids = cKDTree(locs).query(locs, k=k)
    neighbourhoods = locs[neighbour_ids]
    centred = neighbourhoods - neighbourhoods.mean(axis=1, keepdims=True)
    covariances = np.einsum("nki,nkj->nij", centred, centred)
    # eigenvalues come in ascending order: first eigenvector is the normal
    _, eigenvectors = np.linalg.eigh(covariances)
    normals = eigenvectors[:, :, 0]
    # direction from the centroid, used to decide the side of each patch
    outward = locs - locs.mean(axis=0)
    outward_norm = np.linalg.norm(outward, axis=1)
    outward_norm[outward_norm == 0] = 1
    cos_outward = np.einsum("ni,ni->n", normals, outward) / outward_norm
    # graph of neighbouring epitopes, weighted so that the spanning tree
    # follows pairs with nearly parallel normals
    rows = np.repeat(np.arange(n_epitopes), k - 1)
    cols = neighbour_ids[:, 1:].ravel()
    weights = 1 - np.abs(np.einsum("ni,ni->n", normals[rows], normals[cols]))
    # zero weights would be read as missing edges
    graph = csr_matrix((weights + 1e-9, (rows, cols)), shape=(n_epitopes,) * 2)
    tree = minimum_spanning_tree(graph.maximum(graph.T))
    tree = tree + tree.T
    n_patches, patch_ids = connected_components(tree, directed=False)
    for patch in range(n_patches):
        members = np.flatnonzero(patch_ids == patch)
        seed = members[np.argmax(np.abs(cos_outward[members]))]
        if abs(cos_outward[seed]) >= flatness_tolerance:
            flip = cos_outward[seed] < 0
        else:
            flip = normals[seed] @ reference_vector < 0
        if flip:
            normals[seed] *= -1
        order, predecessors = breadth_first_order(tree, seed, directed=False)
        for node in order[1:]:
            if normals[node] @ normals[predecessors[node]] < 0:
                normals[node] *= -1
    return normals


def global_normal_direction(epitope_locs, normal_vector = np.array([0,0,1])):
    """
    Assign the same normal to every epitope.

    Parameters
    ----------
    epitope_locs : numpy.ndarray
        Nx3 array of epitope coordinates (only its length is used).
    normal_vector : numpy.ndarray, optional
        Normal assigned to every epitope, used as given (not normalised).
        Default is [0, 0, 1].

    Returns
    -------
    numpy.ndarray
        Nx3 array with normal_vector in every row.
    """
    normals = np.zeros((len(epitope_locs), 3))
    for i in range(len(epitope_locs)):
        normals[i, :] = normal_vector
    return normals

def coordinates_scaling(
    epitope_array, scaling_factor
):  ## formerly named epitopes_scaling
    """
    Scale coordinates about their centroid.

    Parameters
    ----------
    epitope_array : numpy.ndarray
        Nx3 array of coordinates.
    scaling_factor : float
        Scaling factor; values below 1 shrink the set towards its centroid.

    Returns
    -------
    numpy.ndarray
        Nx3 array of scaled coordinates, with the same centroid as the
        input.
    """
    # epitope_array is a numpy array epitope locations
    # first get the centroid of the epitope array
    # (this input can be the averale location of the atoms in the epitopes)
    centroid0 = np.mean(epitope_array, axis=0)
    # scale down epitopes
    naive_scaled = epitope_array * scaling_factor  # for now this number can be fixed
    # get the centroid of the scaled epitope array
    centroid1 = np.mean(naive_scaled, axis=0)  # this is not necesary for now
    naive_scaled_centered = np.zeros(naive_scaled.shape)
    displacement = centroid1 - centroid0
    for i in range(naive_scaled.shape[0]):
        naive_scaled_centered[i] = naive_scaled[i] - displacement
    return naive_scaled_centered  # output is centered at the centroid of the input


def tilt_direction(direction, tilt_deg, rng=None):
    """
    Tilt a direction by a fixed angle with a random azimuth.

    Parameters
    ----------
    direction : numpy.ndarray
        Direction to tilt (e.g. a surface normal); need not be unit length.
    tilt_deg : float
        Angle between the input and the output direction, in degrees.
    rng : numpy.random.Generator, optional
        Generator for the azimuth. Default is the global numpy random state.

    Returns
    -------
    numpy.ndarray
        Unit vector at tilt_deg from direction, with an azimuth drawn
        uniformly in [0, 2 pi).
    """
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    if not tilt_deg:
        return d
    helper = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(d, helper)
    u = u / np.linalg.norm(u)
    v = np.cross(d, u)
    phi = rng.uniform(0, 2 * np.pi) if rng is not None else np.random.uniform(0, 2 * np.pi)
    t = np.radians(tilt_deg)
    return np.cos(t) * d + np.sin(t) * (np.cos(phi) * u + np.sin(phi) * v)


def add_wobble(pivot, direction, cone_angle_deg=30, length=1):
    """
    Adds a random wobble to the given vector within a cone defined by cone_angle_deg,
    with uniform distribution over the cone's surface area.
    The wobble is relative to the pivot and the direction of the vector.

    :param pivot: (np.array) The pivot point in 3D space
    :param direction: (np.array) The unit vector representing the initial direction from the pivot.
    :param cone_angle_deg: (float) The maximum angle of the cone (in degrees).
    :param length: (float) The length of the translation vector to preserve after wobbling.

    :return: (np.array) The new wobbling direction vector.
    """
    # Ensure direction is a unit vector
    direction = direction / np.linalg.norm(direction)

    # Convert the cone angle from degrees to radians
    cone_angle_rad = np.radians(cone_angle_deg)

    # Sample azimuthal angle uniformly
    phi = np.random.uniform(0, 2 * np.pi)

    # The cone defines a boundary for theta angles to choose from
    cos_theta_max = np.cos(cone_angle_rad)
    # this max theta is the angle between the originall vector
    # and the new one, which is to be place within a cone
    # then this max theta is actually a lower bound
    cos_theta = np.random.uniform(cos_theta_max, 1)
    # but this is the arc in radians, so we get its arcosine
    theta = np.arccos(cos_theta)

    # Random wobble vector in spherical coordinates
    # note that this is the relative wobble to a unitary vector
    x_wobble = np.sin(theta) * np.cos(phi)
    y_wobble = np.sin(theta) * np.sin(phi)
    z_wobble = np.cos(theta)

    # Build the unitary wobble vector
    wobble_vector = np.array([x_wobble, y_wobble, z_wobble])

    # Now, align the wobble vector to the cone defined by the initial direction
    # Find the axis of rotation (cross product of the initial direction and the z-axis)
    rotation_axis = np.cross(np.array([0, 0, 1]), direction)

    # If the direction is already aligned with the z-axis, pick a rotation axis perpendicular to it
    if np.linalg.norm(rotation_axis) < 1e-6:
        rotation_axis = np.cross(np.array([1, 0, 0]), direction)

    # Rotate the wobble vector around the axis defined by the initial direction
    angle = np.arccos(
        np.dot(np.array([0, 0, 1]), direction)
    )  # Angle between the z-axis and the direction vector
    rotated_wobble = rotate_point(wobble_vector, rotation_axis, angle)

    # Adjust the resulting vector to the desired length (this is the wobble vector)
    wobble_translation = rotated_wobble * length
    # Add the pivot and the translation to get the endpoint
    wobble_endpoint = pivot + wobble_translation
    # Return the wobble as an endpoint relative to the pivot
    return wobble_endpoint
