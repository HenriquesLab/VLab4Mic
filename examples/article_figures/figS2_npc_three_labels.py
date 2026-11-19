"""Fig. 2e and Supp. Fig. S2c: nuclear pore ring radius for three labels.

Nup96 (7R5K, C-terminal motif ELAVGSL) labelled with a SNAP-tag, an
anti-GFP nanobody (GFP + nanobody, 6XZF) and an antibody, imaged by SMLM
and STED seen from above (pore axis along z). For each label and modality,
N independent realisations (one pore each) give the ring radius, compared
with the published radii of the reference dataset (Thevathasan et al.
2019): 53.7 +/- 2.1, 55.0 +/- 1.9 and 64.3 +/- 2.6 nm.

The antibody label is modelled as in the experiment: GFP fused to the
Nup96 C-terminus (primary) and an anti-GFP antibody (Antibody template)
bound to an epitope on the GFP barrel (secondary). Nup96-GFP is a genetic
knock-in, so every Nup96 carries GFP (labelling efficiency 1); the antibody
binds with the labelling efficiency of the other labels.

SMLM: the radius and width come from a circle fit to the localisations.
Only pores whose localisations cover at least MIN_CORNERS of the eight
corners (around the pore centre) are fitted, as sparse partial arcs give
unstable circle fits; the number of pores used is reported.

STED (modality parameters of figS_NPC_quantification.py): not localisation
based, so the radius is measured on the image. The pixels within
STED_REGION_NM of the image centre (where the pore is) are thresholded
(Otsu) and a circle is fitted to the positions of the bright pixels
(fit_circle_geometric); the width is the standard deviation of their
distances to the centre.

For each modality, one simulated image and the average image over all
realisations per label (pores always at the image centre) are saved and
shown with a dotted circle of the pore diameter at the Nup96 epitopes
(NPC_DIAMETER_NM).

Set VLAB4MIC_N_REALISATIONS to change the number of realisations.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from vlab4mic import experiments
from vlab4mic.analysis.particle_measures import (
    compute_threshold,
    count_occupied_sectors,
    fit_circle,
    fit_circle_geometric,
    ring_measures,
)

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 44
LABELLING_EFFICIENCY = 0.4  # Supp. Note 4
MIN_CORNERS = 5
NPC_DIAMETER_NM = 107  # pore diameter at the Nup96 epitopes
SCALEBAR_NM = 100
MODALITIES = ["SMLM", "STED"]
# STED as in figS_NPC_quantification.py (Thevathasan et al. 2019)
STED_PARAMS = dict(lateral_resolution_nm=15, axial_resolution_nm=20, psf_voxel_nm=5,
                   depth_of_field_nm=100, pixelsize_nm=15)
STED_EXP_TIME = 0.0002
STED_BINOMIAL_P = 0.7
# STED pixels further than this from the image centre are background
STED_REGION_NM = 100
TARGET = dict(probe_target_type="Sequence", probe_target_value="ELAVGSL")
LABELS = {
    "SNAP-tag": [dict(probe_template="SNAP-tag", probe_name="SNAP-tag", **TARGET)],
    "anti-GFP nanobody": [dict(probe_template="GFP_w_nanobody", probe_name="GFP_w_nanobody", probe_DoL=2, **TARGET)],
    "anti-GFP antibody": [
        dict(probe_template="GFP", probe_name="Nup96-GFP", as_primary=True,
             labelling_efficiency=1.0, **TARGET),
        dict(probe_template="Antibody", probe_name="anti-GFP",
             probe_target_type="Primary", probe_target_value="Nup96-GFP"),
    ],
}
PUBLISHED_RADIUS_NM = {
    "SNAP-tag": (53.7, 2.1),
    "anti-GFP nanobody": (55.0, 1.9),
    "anti-GFP antibody": (64.3, 2.6),
}


def pixelsize_nm(experiment, modality):
    detector = experiment.imaging_modalities[modality]["detector"]
    return detector["pixelsize"] * detector["scale"] / 1e-9


def measure_smlm(result):
    """Ring radius and width from the localisations of one pore."""
    record = list(result["positions"]["SMLM"]["ch0"].values())[0]
    points = record["localisations"]
    if points is None or len(points) < 3:
        return dict(used=False, n_localisations=0)
    # pore centre from the emitters (known structure) to count corners
    centre = fit_circle(record["emitters"])[0]
    corners = count_occupied_sectors(points, 8, centre)
    if corners < MIN_CORNERS:
        return dict(used=False, corners=corners, n_localisations=len(points))
    measures = ring_measures(points)
    return dict(used=True, corners=corners, radius_nm=measures["radius"],
                width_nm=measures["width"], n_localisations=len(points))


def measure_sted(image, pixelsize):
    """Ring radius and width from a circle fit to the bright pixels around the image centre."""
    rows_px, cols_px = np.indices(image.shape)
    xy = np.c_[cols_px.ravel() + 0.5, rows_px.ravel() + 0.5] * pixelsize  # pixel centres, nm
    values = np.asarray(image, dtype=float).ravel()
    image_centre = np.array(image.shape[::-1]) * pixelsize / 2
    # the pore is at the image centre: ignore the background further out
    region = np.linalg.norm(xy - image_centre, axis=1) <= STED_REGION_NM
    xy, values = xy[region], values[region]
    points = xy[values > compute_threshold(values, method="otsu")]
    if len(points) < 3:
        return dict(used=False, n_pixels=len(points))
    centre, radius = fit_circle_geometric(points)
    width = np.linalg.norm(points - centre, axis=1).std()
    return dict(used=True, n_pixels=len(points), radius_nm=radius, width_nm=width)


rows = []
images = {}  # one image per modality and label, preferring a pore used in the fits
image_sums = {}  # sum of the images of all realisations per modality and label
pixelsizes = {}
for m, modality in enumerate(MODALITIES):
    for i, (label, probes) in enumerate(LABELS.items()):
        _, _, experiment = experiments.image_vsample(
            structure="7R5K",
            probe_list=[dict({"labelling_efficiency": LABELLING_EFFICIENCY}, **p) for p in probes],
            multimodal=[modality],
            number_of_particles=1,
            sample_dimensions=[600, 600, 100],
            sample_inital_orientation=[0, 0, 1],
            run_simulation=False,
            clear_experiment=True,
            random_seed=RANDOM_SEED + len(LABELS) * m + i,
        )
        if modality == "STED":
            experiment.update_modality(modality_name="STED", **STED_PARAMS)
            experiment.set_modality_acq(modality_name="STED", exp_time=STED_EXP_TIME)
            experiment.imager.modalities["STED"]["detector"]["noise_model"]["binomial"]["p"] = STED_BINOMIAL_P
        pixelsizes[modality] = pixelsize_nm(experiment, modality)
        experiment.set_virtualsample_params(random_orientations=False, random_rotations=True)
        key = (modality, label)
        for result in experiment.run_replicates(N_REALISATIONS, modality=modality):
            image = np.asarray(result["images"][modality]["ch0"])[0]
            images.setdefault(key, (image, False))
            image_sums[key] = image_sums.get(key, 0) + image.astype(float)
            if modality == "SMLM":
                measures = measure_smlm(result)
            else:
                measures = measure_sted(image, pixelsizes[modality])
            if measures["used"] and not images[key][1]:
                images[key] = (image, True)
            rows.append(dict(modality=modality, label=label, **measures))

results = pd.DataFrame(rows)
used = results[results.used]
summary = used.groupby(["modality", "label"]).agg(
    radius_mean=("radius_nm", "mean"), radius_sd=("radius_nm", "std"),
    width_mean=("width_nm", "mean"), n_used=("radius_nm", "size"),
)
index = pd.MultiIndex.from_product([MODALITIES, list(LABELS)], names=["modality", "label"])
summary = summary.reindex(index)
summary["n_total"] = results.groupby(["modality", "label"]).size().reindex(index)
summary["published_radius"] = [PUBLISHED_RADIUS_NM[k][0] for k in index.get_level_values("label")]
summary["published_sd"] = [PUBLISHED_RADIUS_NM[k][1] for k in index.get_level_values("label")]
out = experiment.output_directory
prefix = os.path.join(out, experiment.date_as_string + "figS2_npc_three_labels_")
results.to_csv(prefix + "particles.csv", index=False)
summary.to_csv(prefix + "summary.csv")
experiment.save_parameters(out, name="figS2_npc_three_labels")


def plot_images_with_npc_diameter(images_per_label, pixelsize, name):
    """Save each image as tif and plot them with the pore diameter at the epitopes."""
    fig, axes = plt.subplots(1, len(images_per_label), figsize=[5 * len(images_per_label), 5])
    for ax, (label, image) in zip(axes, images_per_label.items()):
        tifffile.imwrite(prefix + name + "_" + label.replace(" ", "_") + ".tif", image)
        height, width = image.shape
        ax.imshow(image, cmap="grey", extent=[0, width * pixelsize, height * pixelsize, 0])
        # the pore is at the centre of the field of view
        centre = (width * pixelsize / 2, height * pixelsize / 2)
        ax.add_patch(plt.Circle(centre, NPC_DIAMETER_NM / 2, fill=False, linestyle=":",
                                color="yellow", linewidth=1.5))
        ax.set_title(label)
        ax.set_xticks([])
        ax.set_yticks([])
    # axes are in nm, so the scale bar size is in nm
    axes[-1].add_artist(AnchoredSizeBar(axes[-1].transData, SCALEBAR_NM, f"{SCALEBAR_NM} nm",
                                        "lower right", color="white", frameon=False,
                                        size_vertical=2))
    fig.savefig(prefix + name + "s_npc_diameter.png", dpi=300, bbox_inches="tight")
    plt.close()


for modality in MODALITIES:
    modality_summary = summary.loc[modality]
    fig, ax = plt.subplots(figsize=[5, 5])
    ax.errorbar(modality_summary.published_radius, modality_summary.radius_mean,
                xerr=modality_summary.published_sd, yerr=modality_summary.radius_sd, fmt="o")
    for label, row in modality_summary.iterrows():
        ax.annotate(label, (row.published_radius, row.radius_mean))
    lims = [45, 75]
    ax.plot(lims, lims, "--", color="grey")
    ax.set_xlabel("published ring radius (nm)")
    ax.set_ylabel(f"predicted ring radius, {modality} (nm)")
    fig.savefig(prefix + modality + "_radius.png", dpi=300, bbox_inches="tight")
    plt.close()

    # one simulated image per label
    plot_images_with_npc_diameter({label: images[(modality, label)][0] for label in LABELS},
                                  pixelsizes[modality], modality + "_image")
    # average image per label over all realisations (pore always at the image centre)
    plot_images_with_npc_diameter(
        {label: (image_sums[(modality, label)] / modality_summary.n_total[label]).astype(np.float32)
         for label in LABELS},
        pixelsizes[modality], modality + "_average_image",
    )
print(summary)
