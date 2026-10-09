from scipy.signal import convolve, convolve2d
import numpy as np
from tqdm import tqdm
import math
import skimage
from scipy.ndimage import zoom, gaussian_filter
import copy


# keep
def convolve3D(voxelated_field, psf, method="fft"):
    """
    psf is the 3D superpixelated psf.
    The pixelsize of both psf and voxelated field must be equal
    psf is also assumed to be centered
    """
    result = convolve(voxelated_field, psf, mode="same", method=method)
    return result


"""
Section: convolve in 3D, use zrange to get emitter that are to be considered
"""


# keep
def voxelate_points_withrange(coordinates, bin_pixelsize, ranges):
    """
    ideally, this calculation is done in nm scale
    all variables are assumed to have the same dimensions
    pixelsize refers to the detection
    imsize, in pixels, are the number of pixel per dimensions
    bin_pixelsize is the voxel size of the data after binning
    """

    # voxels of exactly bin_pixelsize from the lower limit of each range;
    # the grid is extended to cover the whole range when the range is not
    # a whole number of voxels
    nbins = [
        int(np.ceil(np.diff(r)[0] / bin_pixelsize - 1e-9)) for r in ranges
    ]
    binrange = [
        (r[0], r[0] + n * bin_pixelsize) for r, n in zip(ranges, nbins)
    ]
    hist3D, edges = np.histogramdd(coordinates, bins=nbins, range=binrange)
    return hist3D


# keep
def voxelize_active_fluorophores_withrange(
    coordinates, photons_in_frame, ranges, bin_pixelsize
):
    """
    - photons_in_frame: is a vector with correspondance 1-1 with the emitters
        coordinates and contains the amount of photons emited in the time frame
    - Coordinates: is a numpy array of shape Nx3, being N the number of
        emitters
    - range: is a list of 3 elemtns, one per dimension (XYZ in that order),
        each entry is a tuple of 2 elements stating
        the minimum and maximum absolute value that will be used to create
        the discrete field
    """
    activated = photons_in_frame > 0
    activated_coordinates = coordinates[activated, :]
    photons_per_emitter = photons_in_frame[activated]
    coordinates_intensity = []
    coordinates_intensity = np.repeat(
        activated_coordinates, repeats=photons_per_emitter, axis=0
    )
    voxel = voxelate_points_withrange(coordinates_intensity, bin_pixelsize, ranges)
    return voxel


## keep
def frame_by_volume_convolution(
    coordinates,
    photon_vector,
    ranges,
    bin_size,
    kernel3D,
    zfocus_slice,
    projection_depth=2,
    activation_only=False,
    asframe=True,
    as_mask=False,
    render_sigma_xy=None,
):
    """
    Generate a 2D image from the 3D field of active emitters.

    By default the voxelised emitters are convolved with the 3D PSF
    (kernel3D) and the slices within projection_depth of the focal slice
    are summed. If render_sigma_xy is given, the points are localisations
    that already carry their localisation error: they are not convolved
    with the PSF, but rendered with a lateral Gaussian of standard
    deviation render_sigma_xy (in voxels; 0 renders a histogram), and the
    slices within projection_depth of the focal slice are summed.

    projection_depth is in units of slices. All units are in the same scale,
    except image dimensions.
    """
    projection_depth_half = int(projection_depth / 2)
    intensity_voxel = voxelize_active_fluorophores_withrange(
        coordinates, photon_vector, ranges, bin_pixelsize=bin_size
    )
    if as_mask:
        frame = np.sum(np.array(intensity_voxel),axis=2)
        return frame
    if activation_only:
        return intensity_voxel
    if render_sigma_xy is not None:
        # localisations: render without PSF convolution
        if render_sigma_xy > 0 and np.sum(photon_vector) > 0:
            convolved_intensity = gaussian_filter(
                np.asarray(intensity_voxel, dtype=float),
                sigma=(render_sigma_xy, render_sigma_xy, 0),
            )
        else:
            convolved_intensity = intensity_voxel
    elif np.sum(photon_vector) > 0:
        convolved_intensity = convolve3D(intensity_voxel, kernel3D)
    else:
        # there's no need to calculate the convolution if there are no photons
        convolved_intensity = intensity_voxel
    subset_of_frames = (
        zfocus_slice - projection_depth_half,
        zfocus_slice + projection_depth_half,
    )
    if projection_depth_half > int(convolved_intensity.shape[2] / 2):
        projection_depth_half = int(convolved_intensity.shape[2] / 2)
    if asframe:
        frame = np.sum(
            np.array(
                convolved_intensity[
                    :, :, zfocus_slice - projection_depth_half : zfocus_slice + projection_depth_half
                ]
            ),
            axis=2,
        )
        return frame
    else:
        return convolved_intensity
    

# keep
def generate_frames_volume_convolution(
    field_scale: float,
    field_coordinates: np.array,
    field_size: tuple,
    field_pixelsizeXY: float,
    field_zfocus: float,
    photons_frames: np.array,
    psf_scale: float,
    psf_array: np.array,
    psf_zstep: float,
    psf_focus_slice: int,
    psf_pixelsizeXY: float,
    asframes=True,
    activation_only=False,
    psf_projection_depth=1,
    as_mask=False,
    render_sigma_xy=None,
    **kwargs,
):
    """
    Generate frames by voxelising emitters and convolving them with the PSF.

    Inputs are taken from the method _homogenise_scales4convolution_modality.
    If render_sigma_xy is given (in PSF voxels), the coordinates are
    localisations and are rendered with that lateral Gaussian instead of
    being convolved with the PSF (see frame_by_volume_convolution).
    """
    if asframes:
        cframes = []
        nframes = np.shape(photons_frames)[1]
        for f in tqdm(range(int(nframes))):
            photons_in_frame_vect = photons_frames[
                :, f
            ]  # extract the vector that corresponds to the current frame
            # calculate depth. This value considers the PSF dimenstions
            coordinates_scaled = field_coordinates
            pixelsize_scaled = float(field_pixelsizeXY)
            zplane = field_zfocus
            depth = (
                psf_array.shape[2] / 2
            ) * psf_zstep  # Assumes the focus plane is at the middle
            ranges = [
                (0, field_size[0] * pixelsize_scaled),
                (0, field_size[1] * pixelsize_scaled),
                (zplane - depth, zplane + depth),
            ]
            frame_n = frame_by_volume_convolution(
                coordinates_scaled,
                photons_in_frame_vect,
                ranges,
                psf_zstep,
                psf_array,
                psf_focus_slice,
                projection_depth=psf_projection_depth,
                as_mask=as_mask,
                render_sigma_xy=render_sigma_xy,
            )
            cframes.append(frame_n)

        stack_zxy = np.array(cframes)
        currentpixelsizes = (psf_pixelsizeXY, psf_pixelsizeXY)
        newpixelsizes = (float(field_pixelsizeXY), float(field_pixelsizeXY))
        binned_sequence = lateral_binning_stack(
            stack_zxy, currentpixelsizes, newpixelsizes, field_size
        )
        return binned_sequence
    else:
        volume_frames = []
        nframes = np.shape(photons_frames)[1]
        for f in tqdm(range(int(nframes))):
            photons_in_frame_vect = photons_frames[
                :, f
            ]  # extract the vector that corresponds to the current frame
            # calculate depth. This value considers the PSF dimenstions
            coordinates_scaled = field_coordinates
            pixelsize_scaled = float(field_pixelsizeXY)
            zplane = field_zfocus
            depth = (
                psf_array.shape[2] / 2
            ) * psf_zstep  # Assumes the focus plane is at the middle
            ranges = [
                (0, field_size[0] * pixelsize_scaled),
                (0, field_size[1] * pixelsize_scaled),
                (zplane - depth, zplane + depth),
            ]
            v_frame = frame_by_volume_convolution(
                coordinates_scaled,
                photons_in_frame_vect,
                ranges,
                psf_zstep,
                psf_array,
                psf_focus_slice,
                projection_depth=psf_projection_depth,
                asframe=asframes,
                activation_only=activation_only
            )
            volume_frames.append(v_frame)
        return volume_frames


# keep
def _overlap_matrix(n_in, size_in, n_out, size_out):
    # W[o, i] is the fraction of input pixel i that falls in output pixel o
    edges_in = np.arange(n_in + 1) * size_in
    edges_out = np.arange(n_out + 1) * size_out
    low = np.maximum(edges_out[:-1, None], edges_in[None, :-1])
    high = np.minimum(edges_out[1:, None], edges_in[None, 1:])
    return np.clip(high - low, 0, None) / size_in


def lateral_binning_stack(stack, currentpixelsizes, newpixelsizes, field_size):
    """
    Bin a stack of images to the detector pixel size, conserving intensity.

    Each output pixel receives the fraction of every input pixel that it
    overlaps, so the detector pixel size does not need to be a multiple of
    the input (PSF voxel) size, and both can be non-integer.

    Parameters
    ----------
    stack : numpy.ndarray
        Images with shape (frames, x, y), at currentpixelsizes.
    currentpixelsizes : tuple of float
        Pixel size of the input along x and y.
    newpixelsizes : tuple of float
        Detector pixel size along x and y.
    field_size : tuple of int
        Number of detector pixels along x and y.

    Returns
    -------
    numpy.ndarray
        Binned images with shape (frames, field_size[0], field_size[1]).
    """
    print("Binning image stack")
    stack = np.asarray(stack)
    nframes = stack.shape[0]
    nx, ny = int(field_size[0]), int(field_size[1])
    weights_x = _overlap_matrix(stack.shape[1], currentpixelsizes[0], nx, newpixelsizes[0])
    weights_y = _overlap_matrix(stack.shape[2], currentpixelsizes[1], ny, newpixelsizes[1])
    binned_images = np.zeros((nframes, nx, ny))
    for i in range(int(nframes)):
        binned_images[i] = weights_x @ stack[i] @ weights_y.T
    return binned_images


# keep
def discretise_field_xyz(coordinates: np.ndarray, binsizes: list, ranges: list):
    """
    Generate a 3D histogram from a set of coordinates.

    Parameters:
    - coordinates (np.ndarray): The coordinates to be discretised.
    - binsizes (list): The size of the bins in each dimension.
    - ranges (list): The range of values in each dimension.

    Returns:
    - hist3D (np.ndarray): The 3D histogram of the discretised field.
    - edges (list): The bin edges for each dimension.
    """
    binrange = [ranges[0], ranges[1], ranges[2]]
    nbins = [
        int(np.diff(ranges[0]) / binsizes[0]),
        int(np.diff(ranges[1]) / binsizes[1]),
        math.ceil(np.diff(ranges[2]) / binsizes[2]),
    ]
    hist3D, edges = np.histogramdd(coordinates, bins=nbins, range=binrange)
    return hist3D, edges


"""
Convolution functions as 2D slices. Convolve plane-wise with 2D PSF.
Methods for optimising lookup convolution
"""


# keep
def _prepare_data_4convolutions(
    field_coordinates: np.array,
    field_size: tuple,
    field_pixelsizeXY: float,
    field_zfocus: float,
    photons_frames: np.array,
    psf_array: np.array,
    psf_zstep: float,
    psf_focus_slice: int,
    psf_pixelsizeXY: float,
    **kwargs,
):
    psf_half_range = (
        min((psf_array.shape[2] - psf_focus_slice), psf_focus_slice) * psf_zstep
    )
    zrange = psf_half_range * 2
    nframes = np.shape(photons_frames)[1]
    n_emitters = np.shape(photons_frames)[0]
    # pixel dimensions for binning
    currentpixelsizes = (psf_pixelsizeXY, psf_pixelsizeXY)
    newpixelsizes = (float(field_pixelsizeXY), float(field_pixelsizeXY))
    # limits to discretise
    ranges_xyz = [
        (0, field_size[0] * field_pixelsizeXY),
        (0, field_size[1] * field_pixelsizeXY),
        (field_zfocus - zrange / 2, field_zfocus + zrange / 2),
    ]
    # #
    bin_pixelsizes_xyz = [psf_pixelsizeXY, psf_pixelsizeXY, psf_zstep]
    data4conv = dict(
        field_coordinates=field_coordinates,
        field_size=field_size,
        nframes=nframes,
        n_emitters=n_emitters,
        bining_lateral_pixelsize=psf_pixelsizeXY,
        zrange=zrange,
        currentpixelsizes=currentpixelsizes,  # used at binning
        newpixelsizes=newpixelsizes,  # used at binning
        field_zfocus=field_zfocus,
        ranges_xyz=ranges_xyz,
        bin_pixelsizes_xyz=bin_pixelsizes_xyz,
        photons_frames=photons_frames,
        psf_array=psf_array,
        psf_focus_slice=psf_focus_slice,
    )
    return data4conv


# keep
def generate_frames_projection_conv_optimised2(
    field_scale: float,
    field_coordinates: np.array,
    field_size: tuple,
    field_pixelsizeXY: float,
    field_zfocus: float,
    photons_frames: np.array,
    psf_scale: float,
    psf_array: np.array,
    psf_zstep: float,
    psf_focus_slice: int,
    psf_pixelsizeXY: float,
    **kwargs,
):
    conv_data = _prepare_data_4convolutions(
        field_coordinates,
        field_size,
        field_pixelsizeXY,
        field_zfocus,
        photons_frames,
        psf_array,
        psf_zstep,
        psf_focus_slice,
        psf_pixelsizeXY,
    )

    convolved_emitters_list = create_convolved_emitters_list(**conv_data)
    #  Then construct every frame by adding their corresponding emitters
    frames = construct_frames(convolved_emitters_list, **conv_data)
    binned_frames = bin_timeseries(frames, **conv_data)
    return np.array(convolved_emitters_list), frames, binned_frames


# keep
def construct_frames(convolved_emitters_list, nframes, photons_frames, **kwargs):
    print("constructing frames")
    timeserie = []
    for f in tqdm(range(int(nframes))):
        # extract the vector that corresponds to the current frame
        photons_in_frame_vect = photons_frames[:, f]
        timeserie.append(
            merge_emissions_lookup(photons_in_frame_vect, convolved_emitters_list)
        )
    return np.array(timeserie)


# keep
def bin_timeseries(timeseires, currentpixelsizes, newpixelsizes, field_size, **kwargs):
    # field_size is equivalent to image size
    # assumes the stack has dimensions tXY
    print(
        f"current pixels size for bining: {currentpixelsizes}; new ones {newpixelsizes}"
    )
    print("Binning images")
    binned_timeserie = lateral_binning_stack(
        timeseires, currentpixelsizes, newpixelsizes, field_size
    )
    return binned_timeserie


# keep
def create_convolved_emitters_list(
    field_coordinates,
    n_emitters,
    psf_array,
    ranges_xyz,
    bin_pixelsizes_xyz,
    psf_focus_slice,
    field_zfocus,
    **kwargs,
):
    """ """
    print("Calculating single-emitter profile")
    convolved_emitters_list = []
    for em in tqdm(range(n_emitters)):
        emitter_coordinate = field_coordinates[em, :]
        row_vector_coordinate = emitter_coordinate[
            np.newaxis, :
        ]  # to preserve dimension of 1X3
        convolved_emitters_list.append(
            single_emitter_convolved_lookup(
                row_vector_coordinate,
                ranges_xyz,
                bin_pixelsizes_xyz,
                psf_array,
                psf_focus_slice,
                field_zfocus,
            )
        )
    return convolved_emitters_list


# keep
def single_emitter_convolved_lookup(
    coordinate,
    ranges_xyz,
    bin_pixelsizes_xyz,
    psf_array,
    psf_focus_slice,
    field_zfocus,
    **kwargs,
):
    """
    coordinate has to have shape (1,3)

    """
    discrete_emitter_volume, edges = discretise_field_xyz(
        coordinate, bin_pixelsizes_xyz, ranges_xyz
    )
    # wchich slice contains the emitter
    where_ = np.where(discrete_emitter_volume > 0)
    slice_of_emiter = where_[2][0]
    # which slice on the discrete version corresponds to the focus
    where_discrete_focus = np.where(np.isclose(edges[2], field_zfocus))
    discrete_focus_idx = where_discrete_focus[0][0]
    # according to this get which would be the corresponding slice of PSF
    psf_slice_id = psf_focus_slice - (discrete_focus_idx - slice_of_emiter)
    discrete_2D_slice = discrete_emitter_volume[:, :, slice_of_emiter]
    psf_slice2D = psf_array[:, :, psf_slice_id]
    # convolve these 2D images
    convolved_emitter = convolve2d(discrete_2D_slice, psf_slice2D, mode="same")
    # Profit
    return convolved_emitter


# keep
def merge_emissions_lookup(photons_in_frame_vect, convolved_emitters_list):
    """
    This method will take the pre-calculated emittsion profile
    from each isolated emitter in the field
    and according to the photon vector will combine the
    emission from the active emitters and its corresponding
    photon flux at that frame
    """
    # exception if there are no activated emitters in frame
    total_photons_frame = np.sum(photons_in_frame_vect)
    frame = np.zeros(convolved_emitters_list[0].shape)
    if total_photons_frame > 0:
        for pos in range(len(photons_in_frame_vect)):
            if photons_in_frame_vect[pos] > 0:
                intensity = convolved_emitters_list[pos] * photons_in_frame_vect[pos]
                frame = np.sum([frame, intensity], axis=0)
    # shall return a single image
    return frame
