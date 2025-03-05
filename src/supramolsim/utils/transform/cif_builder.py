import numpy as np
from .points_transforms import decorate_epitopes_normals
from ..sample.arrays import binomial_epitope_sampling


def do_sym_operation(pts_set, sym_operation):
    # pts_set contains the coordinates of al epitopes/atoms of an asymetric unit
    # sym is a list of 2 objects, the first is a matrix that establishes a rotation
    rot = sym_operation[0]
    translation = sym_operation[1]
    symmetry_transformed = []
    for epi in pts_set:
        rotated = np.matmul(rot, epi)  ## matrix * vector
        rot_trans = rotated + translation
        symmetry_transformed.append(rot_trans)
    pts_set_symm = np.array(symmetry_transformed)
    return pts_set_symm


def sym_transforms_numeric_only(
    assembly_transf_dictionary,
    point_coords,
    by_sym_unit=False,
    show_transforms_used=False,
):
    # point coords is a numpy array with dimensions nx3,
    # where n is the number of points to transform
    # assembly_transf_dictionary is a dicitonary
    # where each key is the name of the transform,
    # and the value of each key a list
    # the first element of the list is the rotation matrix,
    # the second element is the translation vector
    new_points = []
    for op in assembly_transf_dictionary.keys():
        # ignore transformations that have letters
        flag = True
        try:
            # try converting to integer
            int(op)
        except ValueError:
            flag = False
        if show_transforms_used:
            print(op, flag)
        if flag:
            new_points.append(
                do_sym_operation(point_coords, assembly_transf_dictionary[op])
            )
    new_points_per_sym = np.array(new_points)
    if by_sym_unit:
        return new_points_per_sym
    else:
        allatoms = new_points_per_sym.reshape(
            (np.shape(new_points_per_sym)[0]) * (np.shape(new_points_per_sym)[1]), 3
        )
        return allatoms


def summarize_epitope_atoms(atoms_per_epitope, method="average", residue_pos=0):
    # this function takes the output list from get_epitopes_by_sequence()
    # Each element of the list is the colection of atoms that correspond to that epitope
    nepitopes = len(atoms_per_epitope)
    epitopes_summary = np.zeros((nepitopes, 3))
    if method == "average":
        for i in range(nepitopes):
            epitopes_summary[i, :] = np.mean(atoms_per_epitope[i], axis=0)
    if method == "last":
        for i in range(nepitopes):
            epitopes_summary[i, :] = atoms_per_epitope[i][nepitopes - 1,]
    if method == "first":
        for i in range(nepitopes):
            epitopes_summary[i, :] = atoms_per_epitope[i][0,]
    if method == "residue":
        pass
    return epitopes_summary


def array_coords_subset(
    array, percentage
):  # TODO: use sample array function from utils
    # show a percentage of all the atoms contained in the full assembly of a CIF
    # this assumes the shape of array is N x 3, where N is the total number of atoms
    if percentage < 1:
        size = array.shape[0]
        ids = np.random.choice(
            np.arange(0, size), int(size * percentage), replace=False
        )
        subset = array[ids, :]
        return subset
    elif percentage == 1:
        return array


def create_instance_label(
    coords_normals: dict, target_type: str, label_data: dict, **kwargs
):
    # receive points + normals for the target
    # receive the labeling parameters
    """
    This function is intended to be called for each type of label
    in the source of labelling instance

    coords_normals: dictionary as follows
        "coordinates" = numpy.array
        "normals" = None or numpy array
        ...
    """
    if target_type == "direct":
        # it could be that the corresponding label is indirect
        # but would not be able to compute since there are no normals
        return (direct_labelling(coords_normals, label_data), None)
    elif target_type == "indirect":
        return indirect_labelling(coords_normals, label_data)


def direct_labelling(coords_normals, label_data, **kwargs):
    coordinates = coords_normals["coordinates"]
    efficiency = label_data["labelling_efficiency"]
    label_data["minimal_distance"]
    return binomial_epitope_sampling(
        epitopes=coordinates,
        p=efficiency,
        normals=None,
        min_distance=label_data["minimal_distance"]
    )


def indirect_labelling(coords_nomrals, label_data, **kwargs):
    """
    coords_nomrals is a dictionary containing both coordinates and
    the normals asociated with each
    """
    # label data contains emitter coordinates, labeling efficiency and more parameters
    # needed to create variable or static labeling

    # step1: verify that the target type and its label type correspond
    #       otherwise rise and exemption and create a direct labelling
    # step2: if they match sample the available epitopes sterically by
    #       considering the distance imposed by the label
    # step3: then, populate each of those sampled epitopes with the label
    #       coordinates following the orientation
    if label_data["emitters"] is not None:
        # first, if there is a size constraint, sample the maximum number
        # that will fit these constraint

        # second: from those, we take the labeling efficiency into account
        # this is assuming the labeling efficiency is a reflection of
        # the binding energy needed not the steric impediiment of nearby fluorophroes

        # select randomly, and according to labelling_efficiency
        # the places that will be labelled
        #efficiency = label_data["labelling_efficiency"]

        epitopes, normals = binomial_epitope_sampling(
            epitopes=coords_nomrals["coordinates"],
            p=label_data["labelling_efficiency"],
            normals=coords_nomrals["normals"],
            min_distance=label_data["minimal_distance"]
        )  #####################
        normals_ft_epitopes = [normals, epitopes]
        indirect_realisation, list_reoriented_points_normals = decorate_epitopes_normals(
            normals_ft_epitopes, label_data["emitters"]
        )
    else:
        print("No emitters defined in label. Using direct labelling")
        indirect_realisation = direct_labelling(coords_nomrals, label_data)
        list_reoriented_points_normals = None
    return indirect_realisation, list_reoriented_points_normals
