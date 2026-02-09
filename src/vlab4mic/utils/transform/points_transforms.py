import numpy as np
import copy
from scipy.spatial.transform import Rotation as R  

# # Rodrigues Rotation
def rotate_point(point, axis, angle):
    ## implements the Rodrigues' rotation formula
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    axis = axis / np.linalg.norm(axis)
    cross = np.cross(axis, point)
    dot = np.dot(axis, point)
    return cos_theta * point + sin_theta * cross + (1 - cos_theta) * dot * axis


def rotate_set(points, axi, angl):
    rotated_pts = []
    for pt in range(len(points)):
        rotated_pts.append(rotate_point(points[pt], axi, angl))
    rotated_pts = np.array(rotated_pts)
    return rotated_pts


# following function was reorient_set
def labeling_reorient_set(label_points, piv_id, intern_axis_id, direction, end_point):
    # pts are the points to reorient
    # piv_id is the pivot index that will be use to reorient pts
    # intern_axis_id is an index for the other point
    # in pts that will define its internal axix of rotation
    # this also dictates the orientation of the internal axis so that
    # direction sets the new direction all the points need to
    # follow relative to its internal axis
    # ref is the reference point for overall traslation of pts
    # Translate towards the origin
    # translated_origin, displacement = displace_set2(pts,piv_id,[0,0,0])
    translated_origin, displacement = transform_displace_set(
        label_points, label_points[piv_id], np.array([0, 0, 0])
    )
    # define an internal axis of orientation
    int_axis = translated_origin[intern_axis_id] - translated_origin[piv_id]
    # get rotation axis
    V = np.cross(int_axis, direction)
    if np.linalg.norm(V) == 0:
        V = direction
    # get angle or rotation in radians
    theta = np.arccos(
        np.dot(int_axis, direction)
        / (np.linalg.norm(int_axis) * np.linalg.norm(direction))
    )
    # Rotate pts
    rotated_at_origin = rotate_set(translated_origin, V, theta)
    reoriented, _ = transform_displace_set(
        rotated_at_origin, rotated_at_origin[piv_id], end_point
    )
    return reoriented, V


def transform_displace_set(
    points, external_reference, end_point
):  # previously displace_set_3
    # points are the set of points to displace
    # external_reference is the reference point
    # that usually is inside the cloud of points
    ## end_point is where the external_reference will be relocated, so that
    # a translation vector for all points is defined
    transl_vec = end_point - external_reference  # specify the translation vector
    translated = points + transl_vec
    return translated, transl_vec


def rotate_pts_by_vector(pts, current_vector, new_vector, pts_referece_center):
    translated_origin, displacement = transform_displace_set(
        pts, pts_referece_center, np.array([0, 0, 0])
    )
    # get rotation axis by calculating the cross product of the current axis with new
    ax_of_rot = np.cross(current_vector, new_vector)
    if np.linalg.norm(ax_of_rot) == 0:
        ax_of_rot = new_vector
    # get angle or rotation in radians
    theta = np.arccos(
        np.dot(current_vector, new_vector)
        / (np.linalg.norm(current_vector) * np.linalg.norm(new_vector))
    )
    # Rotate pts
    rotated = rotate_set(translated_origin, ax_of_rot, theta)
    reoriented, _ = transform_displace_set(
        rotated, np.array([0, 0, 0]), pts_referece_center
    )
    return reoriented


# method that takes the epitopes with their corresponding normals
# and takes the labeling entity generated with gen_labeling_entity
def decorate_epitopes_normals(normals_ft_epitopes, labeling_entity, dol:int = None, **kwargs):
    """
    labeling_entity is assumed to have the pivot and the axis orientation
    as the first two points, only after that are defined the actual emitters
    they will be removed after decorating the epitopes

    output: the new positioned emitters without pivots
    """
    labeling_entity_copy = copy.copy(labeling_entity)
    pivot_reference_index = 0
    axis_defining_index = 1
    list_reoriented_points = []
    list_reoriented_points_normals = []
    for repl in range(np.shape(normals_ft_epitopes)[1]):
        new_points, rot_axis = labeling_reorient_set(
            labeling_entity_copy,
            pivot_reference_index,
            axis_defining_index,
            normals_ft_epitopes[0][repl],
            normals_ft_epitopes[1][repl],
        )
        # model degree of labelling
        if dol is not None:
            int_dol = np.random.poisson(lam=dol)
            print(f"DOL = {int_dol}")
            max_emitters = len(new_points)
            if int_dol != 0:
                if int_dol > max_emitters - 2 :
                    # use max value
                    list_reoriented_points.append(new_points)
                else:
                    # take the first two points to keep pivot and axis
                    new_points_dol = copy.copy(new_points[0:2])
                    # then randomly select the rest of the emitters
                    emitter_indices = np.arange(2, max_emitters)
                    selected_emitters = np.random.choice(
                        emitter_indices, size=int_dol, replace=False
                    )
                    for se in selected_emitters:
                        new_points_dol = np.vstack((new_points_dol, new_points[se]))
                    list_reoriented_points.append(new_points_dol)
                    #new_points = new_points[0:(2+int_dol)]
                    #list_reoriented_points.append(new_points)
        else:
            # print(new_points, rot_axis)
            list_reoriented_points.append(new_points)
            list_reoriented_points_normals.append(normals_ft_epitopes[0][repl])
    #print(f"before cleaning: {len(list_reoriented_points)}")
    # cleanup label entities and leave only the true emitters
    cleaned_up_labels = cleanup_labeling_entities(list_reoriented_points)
    return cleaned_up_labels, list_reoriented_points_normals


def cleanup_labeling_entities(labels_on_epitopes):
    """
    after repositioning and reotienting labels the list of labels
    is assumed to contain pivot and axis point on positions 0 and 1
    THis method removes those leaving only the true emitters
    """
    n_epitopes = len(labels_on_epitopes)
    # print(f"cleaning {n_epitopes}")
    pivot_reference_index = 0
    # prepare the list with the first entity's emitters
    fluorophores_xyz = labels_on_epitopes[pivot_reference_index][
        2:
    ]  # the first two are to be removed
    for i in range(n_epitopes - 1):
        # joining them like this allows for different number of emitters per antibody
        fluorophores_xyz = np.vstack(
            (fluorophores_xyz, labels_on_epitopes[i + 1][2:])
        )  # from i+1 since we already added the first
    return fluorophores_xyz


def planar_rotation(points, points_reference, unitary_vector_of_rotation, theta):
    """ 
    Rotate a set of points around a given rotation axis.

    Parameters
    ----------
    points : Nx3 array
        Array containing the N points to rotate
    points_reference : 1x3 array
        reference point for the set 
    unitary_vector_of_rotation: 1x3 array
        Vector at the origin to define the rotation axis
    theta: float
        Angle for rotation, in radians

    Returns
    -------
    numpy.ndarray or None
        Emitter coordinates or None if not found.
    """
    translated_origin, displacement = transform_displace_set(
        points, points_reference, np.array([0, 0, 0])
    )
    rotated_at_origin = rotate_set(translated_origin, unitary_vector_of_rotation, theta)
    rotated, _ = transform_displace_set(
        rotated_at_origin, np.array([0, 0, 0]), points_reference
    )
    return rotated
  
def apply_euler_rotation(vector, phi=0, theta=0, psi=0, order = "zyx", reset_orientation=False, **kwargs):
    """
    Rotations are extrinsic.

    
    theta: rotation around x axis in degrees
    phi: rotation around z axis in degrees
    psi: rotation around y axis in degrees
    """
    combined_R = R.from_euler(order, [phi,theta, psi], degrees=True)
    new_vector = combined_R.apply(vector)
    return new_vector