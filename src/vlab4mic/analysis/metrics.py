from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max
from scipy.stats import pearsonr
import cv2

def img_compare(ref, query, metric=["ssim",], force_match=False, zoom_in=0, **kwargs):
    """
    Compare two images using specified similarity metrics.

    Parameters
    ----------
    ref : numpy.ndarray
        Reference image.
    query : numpy.ndarray
        Image to compare against the reference.
    metric : list of str, optional
        List of metrics to use for comparison. Supported: "ssim", "pearson". Default is ["ssim"].
    force_match : bool, optional
        If True, resize images to match pixel sizes before comparison. Default is False.
    zoom_in : float, optional
        Zoom factor for cropping the images before comparison. Default is 0.
    **kwargs
        Additional keyword arguments, e.g., 'ref_pixelsize', 'modality_pixelsize'.

    Returns
    -------
    similarity_vector : list
        List of similarity values for each metric.
    ref : numpy.ndarray
        (Possibly resized) reference image.
    query : numpy.ndarray
        (Possibly resized) query image.
    """
    if force_match:
        if 'ref_pixelsize' in kwargs and 'modality_pixelsize' in kwargs:
            ref, query = resize_images_interpolation(
                img1=ref,
                img2=query,
                px_size_im1=kwargs['ref_pixelsize'],
                px_size_im2=kwargs['modality_pixelsize'],
                zoom_in=zoom_in
            )
        else:
            ref, query = resize_images_interpolation(
                img1=ref,
                img2=query,
                zoom_in=zoom_in
            )
    similarity_vector = []
    for method in metric:
        if method == "ssim":
            similarity = ssim(ref, query, data_range=query.max() - query.min())
            similarity_vector.append(similarity)
        elif method == "pearson":
            similarity, pval = pearsonr(ref.flatten(), query.flatten())
            similarity_vector.append(similarity)
    return similarity_vector, ref, query, 


def _padding(img1, img2, zoom_in=0, **kwargs):
    """
    Pad two images with zeros so they have the same shape, centering the originals.

    Parameters
    ----------
    img1 : numpy.ndarray
        First image.
    img2 : numpy.ndarray
        Second image.
    zoom_in : float, optional
        Zoom factor for cropping after padding. Default is 0.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    tuple of numpy.ndarray
        The padded (and possibly zoomed) images.
    """
    if img1.shape == img2.shape:
        if zoom_in:
            img1 = zoom_img(img1, zoom_in=zoom_in)
            img2 = zoom_img(img2, zoom_in=zoom_in)
        return img1, img2
    else:
        height1, width1 = img1.shape
        height2, width2 = img2.shape
        max_width = max(width1, width2)
        max_height = max(height1, height2)

        # Calculate the padding needed for both images
        padding1_top = (max_height - height1) // 2
        padding1_left = (max_width - width1) // 2

        padding2_top = (max_height - height2) // 2
        padding2_left = (max_width - width2) // 2

        # Create padded arrays (with zero padding)
        img1_padded = np.zeros((max_height, max_width), dtype=np.float32)
        img2_padded = np.zeros((max_height, max_width), dtype=np.float32)

        # Place the original images into the center of the padded arrays
        img1_padded[
            padding1_top : padding1_top + height1,
            padding1_left : padding1_left + width1,
        ] = img1
        img2_padded[
            padding2_top : padding2_top + height2,
            padding2_left : padding2_left + width2,
        ] = img2
        if zoom_in:
            img1_padded = zoom_img(img1_padded, zoom_in=zoom_in)
            img2_padded = zoom_img(img2_padded, zoom_in=zoom_in)
        return img1_padded, img2_padded


def resize_images_interpolation(
    img1, img2, px_size_im1=1, px_size_im2=1, interpolation_order=3, zoom_in=0
):
    """
    Resize and interpolate images to match pixel sizes, then pad as needed.

    Parameters
    ----------
    img1 : numpy.ndarray
        Reference image.
    img2 : numpy.ndarray
        Image to compare against reference.
    px_size_im1 : float, optional
        Pixel size of image 1. Default is 1.
    px_size_im2 : float, optional
        Pixel size of image 2. Default is 1.
    interpolation_order : int, optional
        Order of interpolation (passed to scipy.ndimage.zoom). Default is 3.
    zoom_in : float, optional
        Zoom factor for cropping after resizing. Default is 0.

    Returns
    -------
    tuple of numpy.ndarray
        The resized and padded images.
    """
    pixel_size_ratio = px_size_im2 / px_size_im1

    resized_img2 = zoom(img2, pixel_size_ratio, order=interpolation_order)

    return _padding(img1, resized_img2, zoom_in)

def zoom_img(img, zoom_in=0, **kwargs):
    """
    Crop the center of the image by a zoom factor.

    Parameters
    ----------
    img : numpy.ndarray
        Image to crop.
    zoom_in : float, optional
        Fraction of the image to crop from each side. Default is 0.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    numpy.ndarray
        Cropped image.
    """
    if zoom_in >= 0.9 or zoom_in < 0:
        zoom_in = 0
    if zoom_in:
        imshape = img.shape
        patchsize = np.ceil(imshape[0] - (imshape[0] * zoom_in))
        center_x = int(imshape[0] / 2)
        center_y = int(imshape[1] / 2)
        half_patch = int(patchsize / 2)
        img = img[center_y-half_patch:center_y+half_patch, center_x-half_patch:center_x+half_patch]
    return img


def image_preprocess(img, background=None, sigma=None, **kwargs):
    """
    Preprocess an image by subtracting background and applying Gaussian filter.

    Parameters
    ----------
    img : numpy.ndarray
        Image to preprocess.
    background : float or numpy.ndarray, optional
        Value or image to subtract as background.
    sigma : float, optional
        Standard deviation for Gaussian filter.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    numpy.ndarray
        Preprocessed image.
    """
    if background:
        img = img - background
    if sigma:
        img = gaussian_filter(img, sigma=sigma)
    return img

def local_maxima_positions(img, min_distance=1, threshold=None, **kwargs):
    """
    Find local maxima positions in an image.

    Parameters
    ----------
    img : numpy.ndarray
        Input image.
    min_distance : int, optional
        Minimum number of pixels separating peaks. Default is 1.
    threshold : float, optional
        Minimum intensity of peaks. Default is None.
    **kwargs
        Additional keyword arguments for preprocessing.

    Returns
    -------
    xy : numpy.ndarray
        Array of (row, col) indices of local maxima.
    img_pre : numpy.ndarray
        Preprocessed image used for peak finding.
    """
    if len(img.shape) > 2:
        img = np.mean(img, axis=-1)
    # remove background as offset value
    img_pre = image_preprocess(img, **kwargs)
    xy = peak_local_max(img_pre, min_distance=min_distance, threshold_abs=threshold)
    return xy, img_pre

def pixel_positions_to_relative(indices, image_sizes, pixelsize):
    """
    Convert pixel indices to relative positions in the image.

    Parameters
    ----------
    indices : list or numpy.ndarray
        List of pixel indices.
    image_sizes : tuple or list
        Size of the image (height, width).
    pixelsize : float
        Size of a pixel.

    Returns
    -------
    xyz_relative : list of numpy.ndarray
        List of relative positions (x, y, z=0).
    """
    image_relative_positions=[(np.array(p)*pixelsize)/image_sizes for p in indices]
    xyz_relative = [ np.append(xypos, 0)  for xypos in image_relative_positions]
    return xyz_relative

def get_circles(img,
                blur_px=1,
                dp=0.1,
                minDist=3, 
                param1=1,
                param2=12,
                minRadius=3,
                maxRadius=7,
                **kwargs):
    """
    Detect circles in an image using the Hough Circle Transform.

    Parameters
    ----------
    img : numpy.ndarray
        Input image.
    blur_px : int, optional
        Size of the blur kernel. Default is 1.
    dp : float, optional
        Inverse ratio of the accumulator resolution to the image resolution. Default is 0.1.
    minDist : int, optional
        Minimum distance between detected centers. Default is 3.
    param1 : int, optional
        First method-specific parameter for HoughCircles. Default is 1.
    param2 : int, optional
        Second method-specific parameter for HoughCircles. Default is 12.
    minRadius : int, optional
        Minimum circle radius. Default is 3.
    maxRadius : int, optional
        Maximum circle radius. Default is 7.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    circles : numpy.ndarray or None
        Detected circles.
    gray_blurred : numpy.ndarray
        Blurred grayscale image used for detection.
    cirlce_params : dict
        Parameters used for circle detection.
    """
    gray_blurred = img
    if gray_blurred.min() < 0:
        gray_blurred += -gray_blurred.min()
    if gray_blurred.dtype != np.uint8:
        if gray_blurred.dtype in [np.float32, np.float64]:
            #gray_blurred = np.uint8(gray_blurred * 255)
            gray_blurred = np.uint8(gray_blurred/gray_blurred.max() * 255)
        else:
            #gray_blurred = gray_blurred.astype(np.uint8)
            gray_blurred = (gray_blurred/256).astype('uint8')
    if blur_px:
        gray_blurred = cv2.blur(gray_blurred, (blur_px, blur_px), 0)
    else: 
        gray_blurred = gray_blurred
    circles = cv2.HoughCircles(
        gray_blurred,
        cv2.HOUGH_GRADIENT,
        dp=dp,
        minDist=minDist, 
        param1=param1,
        param2=param2,
        minRadius=minRadius,
        maxRadius=maxRadius,
    )
    cirlce_params=dict(
        dp=dp,
        minDist=minDist, 
        param1=param1,
        param2=param2,
        minRadius=minRadius,
        maxRadius=maxRadius,
    )
    return circles, gray_blurred, cirlce_params