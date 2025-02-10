from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom
import numpy as np


def img_compare(ref, query, metric="ssim", force_match=False, **kwargs):
    if force_match:
        ref, query = resize_images_interpolation(ref, query)
    if metric == "ssim":
        similarity = ssim(ref, query, data_range=query.max() - query.min())
    return similarity, ref, query


def _padding(img1, img2):
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
        # print(f"Same shapes; {img1.shape}, {img2.shape}")
        return img1, img2
    else:
        # print(f"Different shapes; {img1.shape}, {img2.shape}")
        height1, width1 = img1.shape
        height2, width2 = img2.shape
        max_width = max(width1, width2)
        max_height = max(height1, height2)

        # Calculate the padding needed for both images
        padding1_top = (max_height - height1) // 2
        padding1_bottom = max_height - height1 - padding1_top
        padding1_left = (max_width - width1) // 2
        padding1_right = max_width - width1 - padding1_left

        padding2_top = (max_height - height2) // 2
        padding2_bottom = max_height - height2 - padding2_top
        padding2_left = (max_width - width2) // 2
        padding2_right = max_width - width2 - padding2_left

        # Create padded arrays (with zero padding, representing black)
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

        # print(f"new shapes; {img1_padded.shape}, {img2_padded.shape}")
        return img1_padded, img2_padded


def resize_images_interpolation(
    img1, img2, px_size_im1=1, px_size_im2=1, interpolation_order=3
):
    """
    Interpolate image with bigger pixel size by using cubic interpolation.
    Pad images if needed after interpolating image 2

    Args:
        img1 (numpy array): First image.
        img2 (numpy array): Second image.
        px_size_im1 (float): pixel size of image 1
        px_size_im2 (float): pixel size of image 2

    Returns:
        tuple: Resized images.
    """
    pixel_size_ratio = px_size_im2 / px_size_im1
    print(pixel_size_ratio)

    # Resize the smaller image using cubic interpolation
    resized_img2 = zoom(img2, pixel_size_ratio, order=interpolation_order)

    # check if paddin needed
    return _padding(img1, resized_img2)
