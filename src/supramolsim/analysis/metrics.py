from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max
from scipy.stats import pearsonr
import cv2

def img_compare(ref, query, metric="ssim", force_match=False, zoom_in=0, **kwargs):
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
    if metric == "ssim":
        similarity = ssim(ref, query, data_range=query.max() - query.min())
    elif metric == "pearson":
        similarity, pval = pearsonr(ref.flatten(), query.flatten())
    return similarity, ref, query


def _padding(img1, img2, zoom_in=0, **kwargs):
    """
    Padding the images with zeroes.
    Padded images will have same sizes, makes no assumption
    on the final size.
    Padding is done by placing the old image in the center
    of a new array with final dimensions.

    Args:
        img1 (numpy array): First image.
        img2 (numpy array): Second image.

    Returns:
        tuple: Resized images by padding zeroes.
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
    Interpolate image with bigger pixel size by using cubic interpolation.
    Pad images if needed after interpolating image 2

    Args:
        img1 (numpy array): Reference image
        img2 (numpy array): Image to compare against reference
        px_size_im1 (float): pixel size of image 1
        px_size_im2 (float): pixel size of image 2

    Returns:
        tuple: Resized images.
    """
    pixel_size_ratio = px_size_im2 / px_size_im1

    resized_img2 = zoom(img2, pixel_size_ratio, order=interpolation_order)

    return _padding(img1, resized_img2, zoom_in)

def zoom_img(img, zoom_in=0, **kwargs):
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
    if background:
        img = img - background
    if sigma:
        img = gaussian_filter(img, sigma=sigma)
    return img

def local_maxima_positions(img, min_distance=1, threshold=None, **kwargs):
    """
    Find local maximas. Assumes single image, non mask.

    Returns:
    (list): list of pairs of indices for each local maxima found
    """
    if len(img.shape) > 2:
        img = np.mean(img, axis=-1)
    # remove background as offset value
    img_pre = image_preprocess(img, **kwargs)
    xy = peak_local_max(img_pre, min_distance=min_distance, threshold_abs=threshold)
    return xy, img_pre

def pixel_positions_to_relative(indices, image_sizes, pixelsize):
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
    if blur_px:
        gray_blurred = cv2.blur(img, (blur_px, blur_px), 0)
    else: 
        gray_blurred = img
    if gray_blurred.min() < 0:
        gray_blurred += -gray_blurred.min()
    if gray_blurred.dtype != np.uint8:
        if gray_blurred.dtype in [np.float32, np.float64]:
            gray_blurred = np.uint8(gray_blurred * 255)
        else:
            gray_blurred = gray_blurred.astype(np.uint8)
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