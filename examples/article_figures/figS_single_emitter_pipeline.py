"""
Trace a single emitter through the image generation pipeline.

One emitter is placed at the centre of a virtual sample of 500x500x10 nm
and imaged with each modality in MODALITIES, for several photon yields and exposure
times (0 s gives an image without signal). The pipeline is reproduced stage
by stage with the library functions, on a single photon draw:

    emitted photons -> discrete sample (voxelised photons) -> convolved
    volume -> z-projection -> binning to detector pixels -> clip < 0 ->
    each detector noise step -> bit depth clip

For every stage the dtype, value range, total signal and negative pixels are
recorded, plotted and written to a CSV file. The last stage is compared with
the output of Imager.generate_imaging run with the same random seed.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

from vlab4mic import experiments
from vlab4mic.utils.transform import image_convolution as conv
from vlab4mic.utils.transform.noise import add_image_noise

RANDOM_SEED = 24
MODALITIES = ["Widefield", "STED", "Confocal", "AiryScan"]
SAMPLE_DIMENSIONS = [500, 500, 10]  # nm
FLUOROPHORE = "AF647"
EXPOSURES_S = [0, 0.001, 0.01, 0.1]
# photons per second; the intermediate yield is the geometric mean of the
# AF647 and AF488 yields
PHOTON_YIELDS = {
    "AF647-like": 100_000,
    "intermediate": 31_623,
    "AF488-like": 10_000,
}
# negative values above this are floating point round-off (FFT convolution)
ROUNDOFF = -1e-9


def stage_stats(array):
    array = np.asarray(array)
    n_negative = int(np.count_nonzero(array < 0))
    n_roundoff = int(np.count_nonzero((array < 0) & (array >= ROUNDOFF)))
    return dict(
        dtype=str(array.dtype),
        shape="x".join(str(s) for s in array.shape),
        min=float(array.min()),
        max=float(array.max()),
        mean=float(array.mean()),
        sum=float(array.sum()),
        dynamic_range=float(array.max() - array.min()),
        n_unique=int(np.unique(array).size),
        n_negative=n_negative,
        frac_negative=n_negative / array.size,
        n_negative_roundoff=n_roundoff,
        frac_negative_roundoff=n_roundoff / array.size,
    )


def trace_pipeline(imager, modality, exp_time):
    """
    Reproduce Imager.generate_imaging for one channel, one fluorophore and
    one frame, keeping the output of every stage.

    Returns the emission record and a list of (stage name, array) tuples.
    Volumes are kept in 3D; images have shape (frames, x, y).
    """
    params = imager.modalities[modality]
    # emitters in the ROI, relative to the ROI origin (imaging.py:1096-1101)
    emitters = np.array(imager.get_emitters_in_ROI(FLUOROPHORE), dtype=float)
    n_emitters = emitters.shape[0]
    emitters[:, 0] -= imager.roi_params["ranges"][0][0]
    emitters[:, 1] -= imager.roi_params["ranges"][1][0]
    photons_frames, _ = imager.calculate_photons_per_frame(
        modality, FLUOROPHORE, n_emitters, 1, exp_time
    )
    emission = dict(
        n_emitters=n_emitters,
        n_active=int(np.count_nonzero(photons_frames[:, 0] > 0)),
        expected_photons=float(
            imager.get_photon_exposure(modality, exp_time)
            * imager.fluorophore_params[FLUOROPHORE]["photons_per_second"]
            * n_emitters
        ),
        emitted_photons=int(np.sum(photons_frames)),
    )
    field_data, psf_data = imager._homogenise_scales4convolution_modality(
        modality, emitters, photons_frames
    )
    stages = []
    discrete = conv.generate_frames_volume_convolution(
        **field_data, **psf_data, asframes=False, activation_only=True
    )[0]
    stages.append(("discrete sample", discrete))
    convolved = conv.generate_frames_volume_convolution(
        **field_data, **psf_data, asframes=False
    )[0]
    stages.append(("convolved", convolved))
    projected = conv._project_slices(
        convolved, psf_data["psf_focus_slice"], int(params["psf"]["depth"] / 2)
    )
    stages.append(("z-projected", projected))
    binned = conv.lateral_binning_stack(
        projected[np.newaxis],
        (psf_data["psf_pixelsizeXY"], psf_data["psf_pixelsizeXY"]),
        (float(field_data["field_pixelsizeXY"]),) * 2,
        field_data["field_size"],
    )
    stages.append(("binned", binned))
    stack = imager._crop_negative(binned)
    stages.append(("clip < 0", stack))
    noise_model = params["detector"]["noise_model"]
    for step in params["detector"]["noise_order"]:
        stack = add_image_noise(step, stack, **noise_model[step])
        stages.append((step, stack))
    stack = imager._adjust_to_pixel_depth(modality, stack)
    stages.append(("bit depth", stack))
    # the binned image must equal the frame the library renders
    library_frame = conv.generate_frames_volume_convolution(
        **field_data,
        **psf_data,
        asframes=True,
        psf_projection_depth=params["psf"]["depth"],
    )
    emission["binning_matches"] = bool(np.allclose(library_frame, binned))
    return emission, stages


def show_stage(ax, array, name):
    image = np.asarray(array, dtype=float)
    if image.ndim == 3 and name in ("discrete sample", "convolved"):
        image = image.sum(axis=2)  # volume: sum along z
    elif image.ndim == 3:
        image = image[0]  # first frame
    im = ax.imshow(image.T, origin="lower", cmap="magma")
    roundoff = (image < 0) & (image >= ROUNDOFF)
    for mask, colour in ((image < ROUNDOFF, "red"), (roundoff, "deepskyblue")):
        if mask.any():
            ax.imshow(
                np.ma.masked_where(~mask.T, np.ones_like(image.T)),
                origin="lower",
                cmap=ListedColormap([colour]),
                vmin=0,
                vmax=1,
            )
    cax = make_axes_locatable(ax).append_axes("right", size="6%", pad=0.03)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=5)
    ax.set_xticks([])
    ax.set_yticks([])


def plot_modality_yield(runs, modality, yield_name, photon_yield, filename):
    """
    Image strip (rows: exposure times, columns: stages) and stage statistics.
    """
    stage_names = [name for name, _ in runs[0]["stages"]]
    nrows, ncols = len(runs), len(stage_names)
    fig = plt.figure(figsize=(1.9 * ncols, 1.9 * nrows + 3.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.9 * nrows, 3.2], hspace=0.25)
    grid = outer[0].subgridspec(nrows, ncols, wspace=0.45, hspace=0.55)
    for r, run in enumerate(runs):
        for c, (name, array) in enumerate(run["stages"]):
            ax = fig.add_subplot(grid[r, c])
            show_stage(ax, array, name)
            st = run["stats"][c]
            ax.set_title(
                f"{name}\n{st['dtype']}\n[{st['min']:.3g}, {st['max']:.3g}]",
                fontsize=6,
            )
            if c == 0:
                em = run["emission"]
                ax.set_ylabel(
                    f"t = {run['exp_time'] * 1e3:g} ms\n"
                    f"expected {em['expected_photons']:.4g} ph\n"
                    f"emitted {em['emitted_photons']} ph\n"
                    f"emitters {em['n_emitters']} (active {em['n_active']})",
                    fontsize=6,
                )
    bottom = outer[1].subgridspec(1, 3, wspace=0.3)
    x = np.arange(ncols)
    tick_labels = [
        f"{name}\n{st['dtype']}" for name, st in zip(stage_names, runs[0]["stats"])
    ]
    colours = plt.cm.viridis(np.linspace(0, 0.9, nrows))
    ax_range = fig.add_subplot(bottom[0])
    ax_sum = fig.add_subplot(bottom[1])
    ax_neg = fig.add_subplot(bottom[2])
    for run, colour in zip(runs, colours):
        label = f"{run['exp_time'] * 1e3:g} ms"
        stats = run["stats"]
        for key, style in (("max", "-"), ("mean", "--"), ("min", ":")):
            ax_range.plot(
                x, [s[key] for s in stats], style, marker="o", markersize=3,
                color=colour, label=f"{label} {key}",
            )
        ax_sum.plot(
            x, [s["sum"] for s in stats], "-o", markersize=3, color=colour,
            label=label,
        )
        ax_sum.axhline(
            run["emission"]["emitted_photons"], color=colour, ls=":", lw=0.8
        )
        ax_neg.plot(
            x,
            [100 * (s["frac_negative"] - s["frac_negative_roundoff"]) for s in stats],
            "-o", markersize=3, color=colour, label=label,
        )
        ax_neg.plot(
            x, [100 * s["frac_negative_roundoff"] for s in stats], ":",
            marker="x", markersize=3, color=colour,
        )
    ax_range.axhline(0, color="k", ls="--", lw=0.8)
    ax_range.set_yscale("symlog", linthresh=1)
    ax_range.set_title("min / mean / max per stage", fontsize=8)
    ax_range.legend(fontsize=5, ncol=2, frameon=False)
    ax_sum.set_yscale("symlog", linthresh=1)
    ax_sum.set_title("total signal per stage (dotted: emitted photons)", fontsize=8)
    ax_sum.legend(fontsize=6, frameon=False)
    ax_neg.set_title(
        f"negative pixels (%); dotted: round-off ({ROUNDOFF:g} to 0)",
        fontsize=8,
    )
    ax_neg.legend(fontsize=6, frameon=False)
    for ax in (ax_range, ax_sum, ax_neg):
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=60, ha="right", fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        f"{modality}, single emitter, {yield_name} ({photon_yield:,} ph/s). "
        f"Red pixels: value < {ROUNDOFF:g}; blue: round-off ({ROUNDOFF:g} to 0)",
        fontsize=10,
    )
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def plot_summary(all_runs, filename):
    """Final image min, max and negative pixels against exposure time."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    x = np.arange(len(EXPOSURES_S))
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(PHOTON_YIELDS)))
    # one line style and marker per modality
    line_styles = ["-", "--", "-.", ":"]
    markers = ["o", "s", "^", "D", "v", "P"]
    styles = {
        modality: line_styles[i % len(line_styles)] + markers[i % len(markers)]
        for i, modality in enumerate(MODALITIES)
    }
    for modality in MODALITIES:
        for colour, yield_name in zip(colours, PHOTON_YIELDS):
            runs = all_runs[(modality, yield_name)]
            final = [run["stats"][-1] for run in runs]
            label = f"{modality}, {yield_name}"
            axes[0].plot(x, [s["min"] for s in final], styles[modality],
                         color=colour, label=label)
            axes[1].plot(x, [s["max"] for s in final], styles[modality],
                         color=colour, label=label)
            axes[2].plot(x, [100 * s["frac_negative"] for s in final],
                         styles[modality], color=colour, label=label)
    axes[0].axhline(0, color="k", ls="--", lw=0.8)
    axes[0].set_yscale("symlog", linthresh=1)
    axes[0].set_title("final image min")
    axes[1].set_yscale("symlog", linthresh=1)
    axes[1].set_title("final image max")
    axes[2].set_title("final image negative pixels (%)")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([f"{t * 1e3:g}" for t in EXPOSURES_S])
        ax.set_xlabel("exposure time (ms)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[2].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    _, _, experiment = experiments.image_vsample(
        multimodal=MODALITIES,
        sample_dimensions=SAMPLE_DIMENSIONS,
        number_of_particles=1,
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED,
    )
    imager = experiment.imager
    centre = np.array(imager.get_roi_params("ranges")).mean(axis=1)
    imager.emitters_by_fluorophore = {FLUOROPHORE: [centre.tolist()]}
    print(
        f"Single emitter at {centre} (um); "
        f"focus plane {imager.get_roi_params('focus_plane')}"
    )
    prefix = os.path.join(
        experiment.output_directory,
        experiment.date_as_string + "figS_single_emitter_pipeline_",
    )
    rows = []
    all_runs = {}
    seed = RANDOM_SEED
    for modality in MODALITIES:
        for yield_name, photon_yield in PHOTON_YIELDS.items():
            imager.fluorophore_params[FLUOROPHORE]["photons_per_second"] = photon_yield
            runs = []
            for exp_time in EXPOSURES_S:
                seed += 1
                np.random.seed(seed)
                emission, stages = trace_pipeline(imager, modality, exp_time)
                stats = [stage_stats(array) for _, array in stages]
                # the same seed must give the library image
                np.random.seed(seed)
                noisy, _, noiseless, _ = imager.generate_imaging(
                    modality=modality, exp_time=exp_time, noise=True
                )
                library_noisy = np.asarray(noisy["ch0"])
                library_noiseless = np.asarray(noiseless["ch0"])
                final = np.asarray(stages[-1][1])
                binned = dict(stages)["binned"]
                emission["library_noisy_matches"] = bool(
                    np.array_equal(library_noisy, final)
                )
                emission["library_noiseless_matches"] = bool(
                    np.allclose(library_noiseless, binned)
                )
                if not emission["library_noisy_matches"]:
                    print(
                        "  library noisy image differs, max abs difference: "
                        f"{np.max(np.abs(library_noisy - final))}"
                    )
                print(
                    f"{modality:9s} {yield_name:12s} t={exp_time * 1e3:6g} ms | "
                    f"emitted {emission['emitted_photons']:6d} ph | "
                    f"final {stats[-1]['dtype']} "
                    f"[{stats[-1]['min']:g}, {stats[-1]['max']:g}], "
                    f"{100 * stats[-1]['frac_negative']:.1f}% < 0 | "
                    f"library noisy match {emission['library_noisy_matches']}, "
                    f"noiseless match {emission['library_noiseless_matches']}, "
                    f"binning match {emission['binning_matches']}"
                )
                for i, ((name, _), st) in enumerate(zip(stages, stats)):
                    rows.append(dict(
                        modality=modality,
                        yield_name=yield_name,
                        photons_per_second=photon_yield,
                        exp_time_s=exp_time,
                        seed=seed,
                        **{k: emission[k] for k in (
                            "n_emitters", "n_active", "expected_photons",
                            "emitted_photons", "library_noisy_matches",
                        )},
                        stage_index=i,
                        stage=name,
                        **st,
                    ))
                runs.append(dict(
                    exp_time=exp_time, emission=emission, stages=stages,
                    stats=stats,
                ))
            all_runs[(modality, yield_name)] = runs
            plot_modality_yield(
                runs, modality, yield_name, photon_yield,
                f"{prefix}{modality}_{yield_name}.png",
            )
    plot_summary(all_runs, prefix + "summary.png")
    with open(prefix + "stages.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Outputs written with prefix {prefix}")


if __name__ == "__main__":
    main()
