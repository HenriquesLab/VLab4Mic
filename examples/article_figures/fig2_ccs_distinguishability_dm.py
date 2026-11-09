from vlab4mic.analysis import distinguishability
from vlab4mic.sweep_generator import run_parameter_sweep
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import ndimage
import requests
from pathlib import Path
import tempfile

labelling_efficiencies_for_sweep = (0,1,0.25)
exp_time_for_sweep = [0.001,]
sweep_repetitions = 10
output_name="vlab_example_sweep_"
modalities=["STED", "SMLM", "AiryScan", "Confocal", "Widefield"]

dome_model = "https://zenodo.org/records/20719843/files/Dome_model_v2.0.cif"
flat_model = "https://zenodo.org/records/20719843/files/Flat_lattice_model_v2.0.cif"

primary = dict(
    probe_templates = ["Antibody",],
    probe_target_value = ["NNRIA",], 
)
def download(url):
    path = Path(tempfile.gettempdir()) / Path(url).name
    if not path.exists():
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path.write_bytes(response.content)
    return str(path)

dome_model_path = download(dome_model)
flat_model_path = download(flat_model)

FOV_dimensions_XY_nm = 600
smlm_pixelsize = 10
mask_image_dimension = int(FOV_dimensions_XY_nm / smlm_pixelsize)
binary_mask = np.zeros(shape=(mask_image_dimension, mask_image_dimension))

#manually placing them
im_centre = int(binary_mask.shape[0]/2)
binary_mask[im_centre, im_centre] = 1

image_parameters = dict(
        mode="mask",
        pixelsize = smlm_pixelsize, 
        min_distance = 10*smlm_pixelsize,
        nparticles = 1,
    )
plt.imshow(binary_mask)


def _as_2d(image):
    image = np.asarray(image, dtype=float)
    if image.ndim == 3:
        # frames, x, y: integrate the frames
        image = image.sum(axis=0)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2-D image or a stack of frames, got {image.shape}")
    return image

def aspect_ratio_features_2(image, threshold_rel=0.2, smoothing_px=1.0):
    """
    Aspect ratio of a single-particle image (size and rotation invariant).

    The image is smoothed with a Gaussian of smoothing_px pixels (so that
    isolated noise pixels do not count), the background (median) is
    subtracted and pixels below threshold_rel of the maximum are ignored;
    the aspect ratio is the square root of the
    ratio of the largest to the smallest eigenvalue of the intensity
    second-moment matrix, i.e. the ratio of the long to the short principal
    axis of the particle seen in projection (1 for a circular particle).

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, or a stack of frames that is summed.
    threshold_rel : float, optional
        Fraction of the maximum below which pixels are ignored. Default 0.2.
    smoothing_px : float, optional
        Standard deviation of the smoothing Gaussian, in pixels. Default 1;
        0 disables smoothing.

    Returns
    -------
    numpy.ndarray
        One-element feature vector [aspect ratio].
    """
    if image.sum() == 0:
        return np.array([0])
    image = _as_2d(image)
    if smoothing_px:
        image = ndimage.gaussian_filter(image, smoothing_px)
    image = np.clip(image - np.median(image), 0, None)
    if image.max() <= 0:
        return np.array([1.0])
    image = np.where(image >= threshold_rel * image.max(), image, 0)
    weights = image / image.sum()
    x, y = np.indices(image.shape)
    cx, cy = (weights * x).sum(), (weights * y).sum()
    dx, dy = x - cx, y - cy
    cov = np.array(
        [
            [(weights * dx * dx).sum(), (weights * dx * dy).sum()],
            [(weights * dx * dy).sum(), (weights * dy * dy).sum()],
        ]
    )
    low, high = np.sort(np.linalg.eigvalsh(cov))
    if low <= 0:
        return np.array([np.inf if high > 0 else 1.0])
    return np.array([np.sqrt(high / low)])


def image_aspect_ratio_2(
        reference_image=None,
        reference_image_pixelsize_nm=None,
        simulated_image=None,
        simulated_image_pixelsize_nm=None,
        image_mask=None,
        resized_reference_image=None,
        resized_simulated_image=None,
        *args,
        **kwargs,
    ):
        feature = aspect_ratio_features_2(simulated_image)
        return feature[0]

sweep_gen_model1 = run_parameter_sweep(
    structures=[dome_model_path,],
    structure_is_path=True,
    use_reference=False,
    **primary,
    modalities=modalities,
    sweep_repetitions=sweep_repetitions,
    # parameters for sweep
    labelling_efficiency=labelling_efficiencies_for_sweep,
    # reference image
    image4vsample=binary_mask,  # 
    image4vsample_parameters=image_parameters,
    exp_time=exp_time_for_sweep,  # 2 values
    # output and analysis
    output_name=output_name,
    return_generator=True,
    run_analysis=True,
    save_sweep_images=False,  #
    save_analysis_results=False,
)

sweep_gen_model1.add_custom_analysis_metrics(
    custom_metrics=[
        image_aspect_ratio_2,])
sweep_gen_model1.run_analysis()

sweep_gen_model2 = run_parameter_sweep(
    structures=[flat_model_path,],
    structure_is_path=True,
    use_reference=False,
    **primary,
    modalities=modalities,
    sweep_repetitions=sweep_repetitions,
    # parameters for sweep
    labelling_efficiency=labelling_efficiencies_for_sweep,
    # reference image
    image4vsample=binary_mask,  # 
    image4vsample_parameters=image_parameters,
    exp_time=exp_time_for_sweep,  # 2 values
    # output and analysis
    output_name=output_name,
    return_generator=True,
    run_analysis=True,
    save_sweep_images=False,  #
    save_analysis_results=False,
)
sweep_gen_model2.add_custom_analysis_metrics(
    custom_metrics=[
        image_aspect_ratio_2,])
sweep_gen_model2.run_analysis()

df1 = sweep_gen_model1.analysis["dataframes"].copy()
df2 = sweep_gen_model2.analysis["dataframes"].copy()

logistic_reg_dataframe__aspectratio = distinguishability.group_dataframes_1value(
    dataframe1=df1,
    dataframe2=df2,
    main_category="modality_name",
    param_name_1="labelling_efficiency",
    metric_name="image_aspect_ratio_2"
    )

figure = distinguishability.plot_distinguishability_sweep(
    logistic_reg_dataframe__aspectratio, "labelling_efficiency", metric="accuracy")

out = sweep_gen_model2.output_directory
prefix = os.path.join(out, sweep_gen_model2.experiment.date_as_string + "fig2_clathrin_")
figure.savefig( prefix + "distinguishability.png", dpi=300, bbox_inches="tight")
plt.close()