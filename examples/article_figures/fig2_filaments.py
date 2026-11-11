import numpy as np
import matplotlib.pyplot as plt
from vlab4mic import experiments
import copy
import seaborn as sns
import os

structure = "5SYF"
probe_model = "1IGT"          # used for every primary and secondary antibody
labelling_efficiency = 1
random_seed = 24  
primary_1 = dict(
    probe_template="Antibody",
    probe_name="ATubulin_epitopes",
    probe_distance_to_epitope = 0,
    probe_model=probe_model,
    probe_target_type="Sequence",
    probe_target_value="EEVGV",
    labelling_efficiency=1,
)


sample_dimensions = [1000, 1000, 10]
common = dict(
    structure=structure,
    structure_global_normal_orientation="local_plane",   # normals from a local plane fit
    sample_dimensions=sample_dimensions,                    # nm
    number_of_particles=1,                               # a single labelled particle
    modality="SMLM",
    SMLM=dict(
        lateral_precision=5,        # nm
        axial_precision=5,          # nm
        n_localisations=30,
        precision_model="fixed",    # not photon-limited
    ),
    clear_experiment=True,
    random_seed=random_seed,
)

def line_of_points2(
    x_dimension,
    y_dimension,
    spacing,
    npoints=None,
    line_position=(0.5, 0.5),
    orientation=0,
    length=None,
):
    if (npoints is None) == (length is None):
        raise ValueError("Specify exactly one of npoints or length.")
    if length is not None:
        npoints = int(np.floor(length / spacing + 1e-9)) + 1
    angle = np.deg2rad(orientation)
    # signed distances from the centre, symmetric on both sides
    distances = (np.arange(npoints) - (npoints - 1) / 2) * spacing
    points = np.zeros((npoints, 3))
    points[:, 0] = line_position[0] + distances * np.cos(angle) / x_dimension
    points[:, 1] = line_position[1] + distances * np.sin(angle) / y_dimension
    return points

pts = line_of_points2(sample_dimensions[0], sample_dimensions[1], spacing=7, line_position=(0.5, 0.5), orientation=0, length=500)

__, noiseless_, my_experiment = experiments.image_vsample(
    probe_list=[primary_1],
    **common,
    random_orientations=True,
    xz_orientations=[90,],
)
my_experiment.clear_virtual_sample()
my_experiment.set_virtualsample_params(
    particle_positions=pts)
my_experiment.build(modules=["coordinate_field","imager"])

fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(1,2,1, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], emitters_plotsize=1, axis_object=ax1)
ax2 = fig.add_subplot(1,2,2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[0,0,0], emitters_plotsize=1, axis_object=ax2)

n_molecules = len(my_experiment.coordinate_field.molecules)
epitopes_new = np.empty(shape=(0,3))
normals_new = np.empty(shape=(0,3))
for particle_number in range(n_molecules):
    epitopes_new = np.vstack((epitopes_new, my_experiment.coordinate_field.molecules[particle_number].source["targets"]["ATubulin_epitopes"]["coordinates"] ))
    normals_new = np.vstack((normals_new, my_experiment.coordinate_field.molecules[particle_number].source["targets"]["ATubulin_epitopes"]["normals"] ))

assert epitopes_new.shape == (1872,3)

absolute_positions = epitopes_new.copy()
normals = normals_new.copy()
particle_generator_data = dict(
    targets={
        "ATubulin_epitopes":dict(
            coordinates = absolute_positions,
            normals=normals
        )
    },
    scale=1e-9,
    reference_point = my_experiment.coordinate_field.get_field_param("absolute_reference_point").copy(),
    axis=dict(pivot=np.array([0.0, 0.0, 0.0]), direction=np.array([0.0, 0.0, 1.0]))
    )
probe_parameters_custom = copy.deepcopy(my_experiment.particle.labels["ATubulin_epitopes"])

### Filament diameter
n_replicates = 200
labelling_efficiency_primary = 1
labelling_efficiency_secondary = 1

primary_1 = dict(
    probe_template="Antibody",
    probe_name="ATubulin_epitopes",
    probe_distance_to_epitope = 0,
    probe_model=probe_model,
    probe_target_type="Sequence",
    probe_target_value="EEVGV",
    labelling_efficiency=1,
    probe_paratope= "TITCH",
    probe_secondary_epitope = "EKKNW"
)
secondary_1 = dict(
    probe_template="Antibody",
    probe_name="Secondary_L1",
    probe_model=probe_model,
    probe_target_type="Primary",
    probe_target_value="ATubulin_epitopes",
    probe_conjugation_target_info=dict(
        type="Atom_residue",
        value=dict(atoms=["CA"], residues=["LYS"]),
    ),
    probe_paratope= "TITCH", 
    #probe_DoL=2,
    labelling_efficiency=labelling_efficiency_secondary,
)
secondary_2 = dict(
    probe_template="Antibody",
    probe_name="Secondary_L2",
    probe_model=probe_model,
    probe_target_type="Primary",
    probe_target_value="ATubulin_epitopes",
    probe_conjugation_target_info=dict(
        type="Atom_residue",
        value=dict(
            atoms= ["CA", ],
            residues=["GLN", ],
            position=305
            ),
    ),
    probe_DoL=None,
    labelling_efficiency=labelling_efficiency,
)
### First strategy
__, noiseless_, my_experiment_probe1 = experiments.image_vsample(
    probe_list=[primary_1],
    **common
)
my_experiment_probe1.particle.source = dict()
my_experiment_probe1.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe1.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe1.particle.generate_instance()
my_experiment_probe1.build(modules=["coordinate_field", "imager"])
replicates_1 = my_experiment_probe1.run_replicates(n_replicates=n_replicates, modality="SMLM")
## Second strategy
__, noiseless_, my_experiment_probe2 = experiments.image_vsample(
    primary_probe=primary_1,
    secondary_probe=secondary_1,
    **common
)
my_experiment_probe2.particle.source = dict()
my_experiment_probe2.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe2.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe2.particle.generate_instance()
my_experiment_probe2.build(modules=["coordinate_field", "imager"])
replicates_2 = my_experiment_probe2.run_replicates(n_replicates=n_replicates, modality="SMLM")
## Third
__, noiseless_, my_experiment_probe3 = experiments.image_vsample(
    primary_probe=primary_1,
    secondary_probe=secondary_2,
    **common
)
my_experiment_probe3.particle.source = dict()
my_experiment_probe3.particle.load_source(**copy.deepcopy(particle_generator_data))
my_experiment_probe3.set_virtualsample_params(
    sample_dimensions=sample_dimensions,
    update_mode = False 
)
my_experiment_probe3.particle.generate_instance()
my_experiment_probe3.build(modules=["coordinate_field", "imager"])
replicates_3 = my_experiment_probe3.run_replicates(n_replicates=n_replicates, modality="SMLM")


import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation


# Global methods return a single threshold value
_GLOBAL_METHODS = {
    "otsu": filters.threshold_otsu,
    "li": filters.threshold_li,
    "yen": filters.threshold_yen,
    "triangle": filters.threshold_triangle,
    "isodata": filters.threshold_isodata,
    "mean": filters.threshold_mean,
    "minimum": filters.threshold_minimum,
}

# Local (adaptive) methods return a threshold image with the same shape as the input
_LOCAL_METHODS = {
    "local": filters.threshold_local,       # e.g. block_size=35, offset=0
    "niblack": filters.threshold_niblack,   # e.g. window_size=15, k=0.2
    "sauvola": filters.threshold_sauvola,   # e.g. window_size=15, k=0.2
}

AVAILABLE_METHODS = sorted([*_GLOBAL_METHODS, *_LOCAL_METHODS, "multiotsu"])


def compute_threshold(image, method="otsu", method_kwargs=None, multiotsu_index=-1):
    """
    Compute a threshold for a 2D image.

    Parameters
    ----------
    image : 2D array
    method : str or callable
        Name of a skimage threshold method (see AVAILABLE_METHODS) or a callable
        ``f(image, **method_kwargs)`` returning a scalar or an array of image shape.
    method_kwargs : dict, optional
        Extra keyword arguments passed to the threshold function.
    multiotsu_index : int
        For ``method="multiotsu"``: which of the returned thresholds to use
        (-1 = highest, i.e. keep only the brightest class).

    Returns
    -------
    threshold : float or np.ndarray
    """
    kwargs = dict(method_kwargs or {})

    if callable(method):
        return method(image, **kwargs)

    name = method.lower()
    if name in _GLOBAL_METHODS:
        return _GLOBAL_METHODS[name](image, **kwargs)
    if name in _LOCAL_METHODS:
        if name == "local":
            kwargs.setdefault("block_size", 35)  # must be odd
        return _LOCAL_METHODS[name](image, **kwargs)
    if name == "multiotsu":
        thresholds = filters.threshold_multiotsu(image, **kwargs)
        return thresholds[multiotsu_index]

    raise ValueError(
        f"Unknown threshold method '{method}'. Choose from {AVAILABLE_METHODS} "
        f"or pass a callable."
    )


def threshold_mask_and_measure(
    image,
    *,
    method="otsu",           # threshold method name or callable
    method_kwargs=None,      # dict of kwargs for the threshold method
    multiotsu_index=-1,      # which multi-Otsu threshold to use
    threshold=None,          # explicit threshold; overrides `method`
    threshold_scale=1.0,     # multiply the computed threshold
    smooth_sigma=0.0,        # Gaussian pre-smoothing (0 = off)
    invert=False,            # True if objects are darker than background
    min_size=0,              # remove objects smaller than this (pixels)
    fill_holes=False,        # fill holes inside objects
    opening_radius=0,        # morphological opening radius (0 = off)
    clear_border=False,      # drop objects touching the image border
    connectivity=2,          # 1 = 4-connected, 2 = 8-connected
    pixel_size=1.0,          # physical size per pixel (e.g. nm/px)
    extra_properties=None,   # list of additional regionprops names
):
    """
    Create a binary mask from a 2D grey-value image using a selectable
    thresholding method, then label connected objects and measure each one.

    Returns
    -------
    mask : np.ndarray (bool)
    table : pd.DataFrame
        One row per object. Metadata in ``table.attrs``:
        ``method``, ``threshold`` (float, or None for local methods),
        ``n_objects``, ``pixel_size``.
    """
    img = np.asarray(image, dtype=float)
    if img.ndim != 2:
        raise ValueError(f"Expected a 2D array, got shape {img.shape}")

    # --- preprocessing ---
    work = filters.gaussian(img, sigma=smooth_sigma, preserve_range=True) \
        if smooth_sigma > 0 else img

    # --- threshold ---
    if threshold is None:
        threshold = compute_threshold(
            work, method=method, method_kwargs=method_kwargs,
            multiotsu_index=multiotsu_index,
        )
        method_used = method if isinstance(method, str) else getattr(method, "__name__", "custom")
    else:
        method_used = "manual"
    threshold = np.asarray(threshold, dtype=float) * threshold_scale

    mask = work < threshold if invert else work > threshold

    # --- mask cleanup ---
    if opening_radius > 0:
        mask = morphology.binary_opening(mask, morphology.disk(opening_radius))
    if fill_holes:
        mask = ndi.binary_fill_holes(mask)
    if min_size > 0:
        mask = morphology.remove_small_objects(
            mask, min_size=min_size, connectivity=connectivity
        )
    if clear_border:
        mask = segmentation.clear_border(mask)

    # --- labelling ---
    labels = measure.label(mask, connectivity=connectivity)

    # --- measurements ---
    properties = [
        "label",
        "area",
        "perimeter",
        "centroid",
        "bbox",
        "equivalent_diameter_area",
        "major_axis_length",
        "minor_axis_length",
        "orientation",
        "eccentricity",
        "solidity",
        "extent",
        "mean_intensity",
        "max_intensity",
        "min_intensity",
    ]
    if extra_properties:
        properties += [p for p in extra_properties if p not in properties]

    table = pd.DataFrame(
        measure.regionprops_table(labels, intensity_image=img, properties=properties)
    )

    if not table.empty:
        table["circularity"] = (
            4 * np.pi * table["area"] / np.square(table["perimeter"]).replace(0, np.nan)
        )
        table["integrated_intensity"] = table["mean_intensity"] * table["area"]

        table["area_physical"] = table["area"] * pixel_size**2
        table["perimeter_physical"] = table["perimeter"] * pixel_size
        table["equivalent_diameter_physical"] = (
            table["equivalent_diameter_area"] * pixel_size
        )
        table["major_axis_physical"] = table["major_axis_length"] * pixel_size
        table["minor_axis_physical"] = table["minor_axis_length"] * pixel_size

        table = table.rename(columns={"centroid-0": "centroid_y",
                                      "centroid-1": "centroid_x"})

    table.attrs["method"] = method_used
    table.attrs["threshold"] = float(threshold) if threshold.ndim == 0 else None
    table.attrs["n_objects"] = int(labels.max())
    table.attrs["pixel_size"] = pixel_size

    return mask, table, work


def measure_filament_width(
    image,
    pixel_size=1.0,
    search_halfwidth=None,
    row_binning=1,
    background_width=5,
    min_snr=2.0,
    mode="single",
    min_sum = 0.001,
):
    """
    Measure the apparent width (FWHM) of a vertical filament centred in a 2D image.

    Each row (or bin of rows) is treated as a transversal line profile. The FWHM
    is measured per profile with sub-pixel linear interpolation. The FWHM of the
    centre-aligned average profile along the whole filament is also reported.

    Parameters
    ----------
    image : 2D array
        Grayscale image with a bright, roughly vertical filament near the centre column.
    pixel_size : float
        Physical size of a pixel. Widths are returned in these units.
    search_halfwidth : int or None
        Half-width (in px) of the window around the centre column where the peak
        is searched. Defaults to a quarter of the image width.
    row_binning : int
        Number of consecutive rows averaged into one profile (improves SNR).
    background_width : int
        Number of pixels at each edge of the search window used to estimate
        the local background.
    min_snr : float
        Profiles whose (peak - background) / background_std is below this value
        are discarded (NaN).

    Returns
    -------
    dict with keys
        "widths"          : 1D array, FWHM per profile (NaN where invalid)
        "centers"         : 1D array, sub-pixel filament centre per profile (px)
        "rows"            : 1D array, central row index of each profile
        "mean", "median", "std" : statistics of valid per-profile widths
        "n_valid"         : number of valid profiles
        "fwhm_avg_profile": FWHM of the centre-aligned mean profile
        "avg_profile"     : (x, profile) tuple of the aligned mean profile
    """
    img = np.asarray(image, dtype=float)
    if img.ndim != 2:
        raise ValueError("image must be 2D")
    n_rows, n_cols = img.shape

    if image.sum() < min_sum:
        return  {
            "widths": None,
            "centers": None,
            "rows": None,
            "mean": None,
            "median": None,
            "std": None,
            "n_valid": None,
            "fwhm_avg_profile": None,
            "avg_profile": None
        }
    if search_halfwidth is None:
        search_halfwidth = max(n_cols // 4, 3)
    c0 = n_cols // 2
    x_lo = max(c0 - search_halfwidth, 0)
    x_hi = min(c0 + search_halfwidth + 1, n_cols)
    x = np.arange(x_lo, x_hi, dtype=float)

    # Bin rows
    n_bins = n_rows // row_binning
    binned = img[: n_bins * row_binning, x_lo:x_hi]
    binned = binned.reshape(n_bins, row_binning, -1).mean(axis=1)
    rows = np.arange(n_bins) * row_binning + (row_binning - 1) / 2.0

    widths = np.full(n_bins, np.nan)
    centers = np.full(n_bins, np.nan)
    aligned = []

    for i, prof in enumerate(binned):
        res = _fwhm_1d(prof, x, background_width, min_snr, mode)
        if res is None:
            continue
        widths[i], centers[i], prof_bg = res
        aligned.append((x - centers[i], prof_bg))

    valid = widths[~np.isnan(widths)] * pixel_size
    widths *= pixel_size

    # Centre-aligned mean profile over the whole filament
    fwhm_avg = np.nan
    avg_profile = (None, None)
    if aligned:
        step = 0.25  # sub-pixel sampling grid
        x_al = np.arange(-search_halfwidth, search_halfwidth + step, step)
        try:
            stack = np.array(
                [np.interp(x_al, xs, p, left=np.nan, right=np.nan) for xs, p in aligned]
            )
            mean_prof = np.nanmean(stack, axis=0)
            ok = ~np.isnan(mean_prof)
            x_al, mean_prof = x_al[ok], mean_prof[ok]
            # Profile is already background-subtracted -> baseline 0
            res = _fwhm_1d(mean_prof, x_al, background_width=0, min_snr=0, mode=mode)
            if res is not None:
                fwhm_avg = res[0] * pixel_size
            avg_profile = (x_al * pixel_size, mean_prof)
        except:
            avg_profile = (None, None)


    return {
        "widths": widths,
        "centers": centers,
        "rows": rows,
        "mean": np.mean(valid) if valid.size else np.nan,
        "median": np.median(valid) if valid.size else np.nan,
        "std": np.std(valid) if valid.size else np.nan,
        "n_valid": int(valid.size),
        "fwhm_avg_profile": fwhm_avg,
        "avg_profile": avg_profile,
    }


def _fwhm_1d(profile, x, background_width, min_snr, mode="single"):
    """
    FWHM of a 1D profile using linear interpolation at half maximum.

    mode="single": width of the brightest peak (walk outward from the maximum
                   until the profile first drops below half max).
    mode="total" : total width of the structure, from the outermost left
                   half-max crossing to the outermost right one. Several
                   peaks are treated as one structure.

    Returns (width, center, background_subtracted_profile) or None if invalid.
    """
    if mode not in ("single", "total"):
        raise ValueError("mode must be 'single' or 'total'")

    p = np.asarray(profile, dtype=float)
    n = p.size

    if background_width > 0 and 2 * background_width < n:
        edges = np.concatenate([p[:background_width], p[-background_width:]])
        bg, bg_std = np.median(edges), np.std(edges)
    else:
        bg, bg_std = 0.0, 0.0

    p = p - bg
    i_max = int(np.argmax(p))
    amp = p[i_max]
    if amp <= 0:
        return None, None, None
    if bg_std > 0 and amp / bg_std < min_snr:
        return None, None, None

    half = amp / 2.0

    if mode == "single":
        # Walk outward from the peak until the profile drops below half max
        i = i_max
        while i > 0 and p[i] >= half:
            i -= 1
        if p[i] >= half:
            return None
        j = i_max
        while j < n - 1 and p[j] >= half:
            j += 1
        if p[j] >= half:
            return None
        i_in, i_out = i + 1, i      # left crossing between i and i+1
        j_in, j_out = j - 1, j      # right crossing between j-1 and j
    else:
        # Outermost samples at or above half max
        above = np.flatnonzero(p >= half)
        i_first, i_last = above[0], above[-1]
        if i_first == 0 or i_last == n - 1:
            return None  # structure touches the window edge
        i_in, i_out = i_first, i_first - 1
        j_in, j_out = i_last, i_last + 1

    x_left = x[i_out] + (half - p[i_out]) * (x[i_in] - x[i_out]) / (p[i_in] - p[i_out])
    x_right = x[j_in] + (half - p[j_in]) * (x[j_out] - x[j_in]) / (p[j_out] - p[j_in])

    width = x_right - x_left
    center = 0.5 * (x_left + x_right)
    return width, center, p
def simple_width(image):
    profile_ = image.sum(axis=0)
    x = np.arange(0, 500, dtype=float)
    width, center, p = _fwhm_1d(profile_, x, 5, 0, mode="total")
    return width

def measure_filaments(replicates, 
                      max_length_nm=500, 
                      pixel_size=1, 
                      modality = "SMLM", 
                      channel="ch0", 
                      fill_holes = True, 
                      smooth_sigma=1, 
                      row_binning=1, 
                      min_snr=2, 
                      mode="single", 
                      threshold_method = "otsu",
                      threshold_method_kwargs = None,):
    tables = []
    masks = []
    major_axis_lengths = []
    filament_metrics_list = []
    rows = []
    blurred_list = []
    for result in replicates:
        mask, table, blurred = threshold_mask_and_measure(result["images"][modality][channel][0], method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = measure_filament_width(result["images"][modality][channel][0], pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = simple_width(result["images"][modality][channel][0])
        filament_metrics_list.append(filament_metrics)
        blurred_list.append(blurred)
        tables.append(table)
        masks.append(mask)
        n_segments = table.shape[0]
        if n_segments == 0:
            min_eccentricity = 0
            max_eccentricity = 0
            mean_eccentricity = 0
            major_lenght = 0
            major_axis_lengths.append(major_lenght)
            apparent_continuity = 0
        else:
            min_eccentricity = np.min(table["eccentricity"])
            max_eccentricity = np.max(table["eccentricity"])
            if len(list(table["eccentricity"])) > 1:
                mean_eccentricity = np.mean(table["eccentricity"])
            else:
                    mean_eccentricity = table["eccentricity"][0]
            major_lenght = table.iloc[table['major_axis_length'].idxmax()]["major_axis_length"]
            major_axis_lengths.append(major_lenght)
            apparent_continuity = (major_lenght / max_length_nm)/ n_segments

        rows.append(dict(
            replicate=result["replicate"],
            mean_filament_width = filament_metrics["mean"],
            std_filament_width = filament_metrics["std"],
            median_filament_width = filament_metrics["median"],
            mask_continuity=apparent_continuity,
            mask_major_lenght = major_lenght,
            mask_number_of_segments = n_segments,
            mask_prefilter_smooth_sigma=smooth_sigma,
            pixelsize = pixel_size,
            row_binning=row_binning,
            s_width=s_width,
            min_eccentricity=min_eccentricity,
            max_eccentricity=max_eccentricity,
            mean_eccentricity=mean_eccentricity
        ))
    df = pd.DataFrame(rows)
    return df, filament_metrics_list, masks, blurred_list, tables

filaments_params = dict(
    pixel_size=2, 
    max_length_nm=500, 
    modality="SMLM", 
    channel="ch0", 
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=2,
    threshold_method = "otsu"
)

df1, filament_metrics_list1, masks1, blurred1, tables1 = measure_filaments(
    replicates=replicates_1, 
    **filaments_params
)
df2, filament_metrics_list2, masks2, blurred2, tables2 = measure_filaments(
    replicates=replicates_2, 
    **filaments_params
    )
df3, filament_metrics_list3, masks3, blurred3, tables3 = measure_filaments(
    replicates=replicates_3, 
    **filaments_params
    )

combined_df_sw = pd.DataFrame(data={
    "Primary": df1["s_width"],
    "Primary_Secondary_SiteSpec": df3["s_width"],
    "Primary_Secondary_NHS": df2["s_width"],
    }
)
custom_palette = {
    "Primary": "#F58B19",
    "Primary_Secondary_SiteSpec": "green",
    "Primary_Secondary_NHS": "#a34bb0",
}

filament_diameters = sns.histplot(data=combined_df_sw, palette=custom_palette, bins=30)
fig = filament_diameters.get_figure()
filename = my_experiment_probe3.date_as_string + 'vlab4mic_apparent_diameter_microtubule.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 