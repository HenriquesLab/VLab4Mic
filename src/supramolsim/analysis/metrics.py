from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom


def img_compare(ref, query, metric="ssim", force_match=False, **kwargs):
    if force_match:
        ref, query = resize_images_interpolation(ref, query)
    if metric == "ssim":
        similarity = ssim(ref, query, data_range=query.max() - query.min())
    return similarity, ref, query


def _padding(img1, img2):
    # check if padding needed
    return img1, img2


def resize_images_interpolation(img1, img2, px_size_im1=1, px_size_im2=1, interpolation_order=3):
    """
    Resize the smaller image to match the dimensions of the larger image using cubic interpolation.
    
    Args:
        img1 (numpy array): First image.
        img2 (numpy array): Second image.
        px_size_im1 (float): pixel size of image 1
        px_size_im2 (float): pixel size of image 2
    
    Returns:
        tuple: Resized images.
    """
    pixel_size_ratio = px_size_im2 / px_size_im1
    
    # Resize the smaller image using cubic interpolation
    resized_smaller_img = zoom(img2, pixel_size_ratio, order=interpolation_order)
    
    # check if paddin needed
    return _padding(img1, img2)
    
