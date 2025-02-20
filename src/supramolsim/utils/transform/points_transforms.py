import numpy as np


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
def decorate_epitopes_normals(normals_ft_epitopes, labeling_entity):
    """
    labeling_entity is assumed to have the pivot and the axis orientation
    as the first two points, only after that are defined the actual emitters
    they will be removed after decorating the epitopes

    output: the new positioned emitters without pivots
    """
    #
    pivot_reference_index = 0
    axis_defining_index = 1
    list_reoriented_points = []
    list_reoriented_points_normals = []
    for repl in range(np.shape(normals_ft_epitopes)[1]):
        new_points, rot_axis = labeling_reorient_set(
            labeling_entity,
            pivot_reference_index,
            axis_defining_index,
            normals_ft_epitopes[0][repl],
            normals_ft_epitopes[1][repl],
        )
        # print(new_points, rot_axis)
        list_reoriented_points.append(new_points)
        list_reoriented_points_normals.append(rot_axis)
    print(f"before cleaning: {len(list_reoriented_points)}")
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
    print(f"cleaning {n_epitopes}")
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
