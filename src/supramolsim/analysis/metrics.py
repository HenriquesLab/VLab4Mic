from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom


def img_compare(ref, query, metric="ssim", force_match=False, **kwargs):
    if force_match:
        ref, query = resize_images_interpolation(ref, query)
    if metric == "ssim":
        similarity = ssim(ref, query, data_range=query.max() - query.min())
        return similarity


def resize_images_interpolation(img1, img2, interpolation_order=3):
    """
    Resize the smaller image to match the dimensions of the larger image using cubic interpolation.
    
    Args:
        img1 (numpy array): First image.
        img2 (numpy array): Second image.
    
    Returns:
        tuple: Resized images.
    """
    # Determine which image is smaller
    if img1.shape[0] * img1.shape[1] > img2.shape[0] * img2.shape[1]:
        larger_img = img1
        smaller_img = img2
    else:
        larger_img = img2
        smaller_img = img1
    
    # Calculate the zoom factors
    zoom_factors = (larger_img.shape[0] / smaller_img.shape[0], larger_img.shape[1] / smaller_img.shape[1])
    
    # Resize the smaller image using cubic interpolation
    resized_smaller_img = zoom(smaller_img, zoom_factors, order=interpolation_order)
    
    # Ensure both images are now the same size
    if larger_img is img1:
        return larger_img, resized_smaller_img
    else:
        return resized_smaller_img, larger_img