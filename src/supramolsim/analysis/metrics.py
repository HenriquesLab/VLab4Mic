from skimage.metrics import structural_similarity as ssim


def img_compare(ref, query, metric="ssim", **kwargs):
    if metric == "ssim":
        similarity = ssim(ref, query, data_range=query.max() - query.min())
        return similarity
