import numpy as np


def normals_by_scaling(epitope_locs, scale = 0.95):
    scaled_epitopes_locs = coordinates_scaling(epitope_locs, scale)
    # calculate the normals of the epitopes
    # note that this approach will work for convex shapes
    # this assumes as well the scaled version is smaller
    normals = np.zeros((len(epitope_locs), 3))
    for i in range(len(epitope_locs)):
        normals[i, :] = epitope_locs[i] - scaled_epitopes_locs[i]
    # then the points from which the normals are traced are the epitopes themselves
    #normals_ft_epitopes = [normals, epitope_locs]
    return normals
    #return normals_ft_epitopes


def coordinates_scaling(epitope_array, scaling_factor):  ## formerly named epitopes_scaling
    # epitope_array is a numpy array of the final localtions of the final epitope locations
    # first get the centroid of the epitope array (this input can be the averale location of the atoms in the epitopes)
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