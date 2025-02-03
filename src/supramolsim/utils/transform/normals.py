import numpy as np
from .points_transforms import rotate_point


def normals_by_scaling(epitope_locs, scale=0.95):
    scaled_epitopes_locs = coordinates_scaling(epitope_locs, scale)
    # calculate the normals of the epitopes
    # note that this approach will work for convex shapes
    # this assumes as well the scaled version is smaller
    normals = np.zeros((len(epitope_locs), 3))
    for i in range(len(epitope_locs)):
        normals[i, :] = epitope_locs[i] - scaled_epitopes_locs[i]
    # then the points from which the normals are traced are the epitopes themselves
    # normals_ft_epitopes = [normals, epitope_locs]
    return normals
    # return normals_ft_epitopes


def coordinates_scaling(
    epitope_array, scaling_factor
):  ## formerly named epitopes_scaling
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


def add_wobble(pivot, direction, cone_angle_deg=30, length=1):
    """
    Adds a random wobble to the given vector within a cone defined by cone_angle_deg,
    with uniform distribution over the cone's surface area.
    The wobble is relative to the pivot and the direction of the vector.

    :param pivot: (np.array) The pivot point in 3D space
    :param direction: (np.array) The unit vector representing the initial direction from the pivot.
    :param cone_angle_deg: (float) The maximum angle of the cone (in degrees).
    :param length: (float) The length of the translation vector to preserve after wobbling.

    :return: (np.array) The new wobbling direction vector (translation vector).
    """
    # Ensure direction is a unit vector
    direction = direction / np.linalg.norm(direction)

    # Convert the cone angle from degrees to radians
    cone_angle_rad = np.radians(cone_angle_deg)

    # Sample azimuthal angle niformly
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
