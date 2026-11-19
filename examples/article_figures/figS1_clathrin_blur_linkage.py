"""Supp. Fig. S1: dome vs flat clathrin lattices against PSF width and linkage error.

Plots the axial_blur, lateral_blur and linkage_error sweeps of the Fig. 2
analysis (fig2_clathrin_distinguishability.py): no new simulation.
Lattices in ventral geometry, 100% labelling, structural integrity 1,
N = 100 realisations per lattice, aspect ratio + brightness features.

a  STED, axial PSF width 50-600 nm at lateral 20 nm.
b  STED, lateral PSF width 10-130 nm at axial 300 nm; widefield (lateral
   130 nm, axial 400 nm, 100 nm pixels) shown for comparison.
c  Linkage error 2-20 nm for every imaging method (10 nm: reference).

Run:   python figS1_clathrin_blur_linkage.py <...fig2_clathrin_distinguishability.csv>
Outputs, next to the input: <prefix>figS1_blur_linkage.{pdf,svg,png} and
<prefix>figS1_blur_linkage_source_data.csv.
"""

import sys

import matplotlib.pyplot as plt
import pandas as pd

MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]
LABELS = {"AiryScan": "Airyscan"}
# categorical slots 1-5 of the validated reference palette, in fixed order,
# with a distinct marker per method so identity is not colour alone
COLOURS = dict(zip(MODALITIES, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]))
MARKERS = dict(zip(MODALITIES, ["o", "s", "^", "D", "v"]))
INK, MUTED, GRID = "#1a1a1a", "#666666", "#d9d9d9"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "text.color": INK, "svg.fonttype": "none", "pdf.fonttype": 42,
})


def errorbars(ax, x, data, colour, marker, label=None, line=True):
    ax.errorbar(x, data.auc, yerr=[data.auc - data.auc_low, data.auc_high - data.auc],
                fmt=marker + ("-" if line else ""), color=colour, markersize=4.5,
                linewidth=1.5, elinewidth=1, capsize=0, label=label,
                markeredgecolor="white", markeredgewidth=0.6, zorder=3)


def null(ax, x, data):
    ax.scatter(x, data.null_auc, marker="x", s=14, color=MUTED, linewidths=0.9, zorder=2,
               label="null (dome vs dome)")


def style(ax, xlabel, title):
    ax.axhline(0.5, color=MUTED, linestyle="--", linewidth=0.8, zorder=1)
    ax.set_ylim(0.3, 1.04)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def main(path):
    scores = pd.read_csv(path)
    reference = scores[scores.sweep == "reference"]
    sted = "STED"
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.3), sharey=True)

    axial = scores[scores.sweep == "axial_blur"].sort_values("axial_resolution_nm")
    errorbars(axes[0], axial.axial_resolution_nm, axial, COLOURS[sted], MARKERS[sted], "STED")
    null(axes[0], axial.axial_resolution_nm, axial)
    axes[0].set_xscale("log")
    axes[0].set_xticks([50, 100, 200, 300, 600], ["50", "100", "200", "300", "600"])
    axes[0].minorticks_off()
    style(axes[0], "axial PSF width, s.d. (nm)", "a  Axial width (STED, lateral 20 nm)")
    axes[0].set_ylabel("Distinguishability (ROC AUC)")

    lateral = scores[scores.sweep == "lateral_blur"].sort_values("lateral_resolution_nm")
    errorbars(axes[1], lateral.lateral_resolution_nm, lateral, COLOURS[sted], MARKERS[sted], "STED")
    null(axes[1], lateral.lateral_resolution_nm, lateral)
    widefield = reference[reference.modality_name == "Widefield"]
    errorbars(axes[1], widefield.lateral_resolution_nm, widefield, COLOURS["Widefield"],
              MARKERS["Widefield"], "Widefield (100 nm pixels)", line=False)
    axes[1].annotate("Widefield", (widefield.lateral_resolution_nm.iloc[0], widefield.auc.iloc[0]),
                     xytext=(-6, 0), textcoords="offset points", ha="right", va="center",
                     fontsize=6.5, color=INK)
    axes[1].annotate("STED", (lateral.lateral_resolution_nm.iloc[-1], lateral.auc.iloc[-1]),
                     xytext=(-6, 5), textcoords="offset points", ha="right", va="bottom",
                     fontsize=6.5, color=INK)
    axes[1].set_xscale("log")
    axes[1].set_xticks([10, 20, 40, 80, 130], ["10", "20", "40", "80", "130"])
    axes[1].minorticks_off()
    style(axes[1], "lateral PSF width, s.d. (nm)", "b  Lateral width (axial 300 nm)")

    linkage = pd.concat([scores[scores.sweep == "linkage_error"], reference])
    for k, modality in enumerate(MODALITIES):
        sub = linkage[linkage.modality_name == modality].sort_values("linkage_error_nm")
        dodge = (k - (len(MODALITIES) - 1) / 2) * 0.3  # nm; overlapping methods stay visible
        errorbars(axes[2], sub.linkage_error_nm + dodge, sub, COLOURS[modality],
                  MARKERS[modality], LABELS.get(modality, modality))
    axes[2].set_xticks([2, 5, 10, 20])
    style(axes[2], "linkage error (nm)", "c  Linkage error")
    axes[2].legend(frameon=False, loc="lower right", ncol=2, handlelength=1.2,
                   columnspacing=0.8, borderaxespad=0.2)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, frameon=False, loc="lower left", handlelength=1.2)

    fig.tight_layout(w_pad=1.0)
    prefix = path.replace("fig2_clathrin_distinguishability.csv", "")
    for extension in ("pdf", "svg", "png"):
        fig.savefig(f"{prefix}figS1_blur_linkage.{extension}", dpi=600, bbox_inches="tight")
    plt.close(fig)

    columns = ["sweep", "modality_name", "axial_resolution_nm", "lateral_resolution_nm",
               "linkage_error_nm", "n_a", "n_b", "auc", "auc_low", "auc_high", "null_auc",
               "p_value", "aspect_ratio_auc", "brightness_auc"]
    source = pd.concat([axial, lateral, widefield.assign(sweep="reference (widefield)"),
                        linkage.sort_values(["modality_name", "linkage_error_nm"])])
    source[columns].to_csv(f"{prefix}figS1_blur_linkage_source_data.csv", index=False)
    print(f"saved {prefix}figS1_blur_linkage.pdf / .svg / .png and source data")


if __name__ == "__main__":
    main(sys.argv[1])
