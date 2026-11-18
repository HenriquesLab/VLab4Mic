import numpy as np
import matplotlib.pyplot as plt
from vlab4mic import experiments
import copy
import seaborn as sns
import os
from vlab4mic.analysis import particle_measures
import pandas as pd
from scipy import optimize, stats

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

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,2,1, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[90,0,0], emitters_plotsize=1, axis_object=ax1)
ax2 = fig.add_subplot(1,2,2, projection="3d")
my_experiment.coordinate_field.show_field(view_init=[0,0,0], emitters_plotsize=1, axis_object=ax2)
filename = my_experiment.date_as_string + 'vlab4mic_microtubule_filament.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 


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
n_replicates = 100
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

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,3,1)
ax1.imshow(replicates_1[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax1.set_title("Primary Antibody")
ax2 = fig.add_subplot(1,3,2)
ax2.imshow(replicates_2[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax2.set_title("Primary Antibody + Secondary Antibody (NHS)")
ax3 = fig.add_subplot(1,3,3)
ax3.imshow(replicates_3[0]["images"]["SMLM"]["ch0"][0], cmap="grey")
ax3.set_title("Primary Antibody + Secondary Antibody (site-specific)")
filename = my_experiment_probe3.date_as_string + 'vlab4mic_microtubules_different_probes.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight') 

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
    eccentricities = []
    major_axis_lengths = []
    filament_metrics_list = []
    rows = []
    blurred_list = []
    for result in replicates:
        mask, table, blurred = particle_measures.threshold_mask_and_measure(result["images"][modality][channel][0], method=threshold_method, smooth_sigma=smooth_sigma, fill_holes=fill_holes, method_kwargs=threshold_method_kwargs)
        filament_metrics = particle_measures.measure_filament_width(result["images"][modality][channel][0], pixel_size=pixel_size, row_binning=row_binning, min_snr=min_snr, mode=mode)
        s_width = particle_measures.simple_width(result["images"][modality][channel][0])
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

# simple_width is in pixels: convert to nm with the pixel size of the SMLM image
combined_df_sw = pd.DataFrame(data={
    "Primary": df1["s_width"],
    "Primary_Secondary_SiteSpec": df3["s_width"],
    "Primary_Secondary_NHS": df2["s_width"],
    }
) / filaments_params["pixel_size"]
custom_palette = {
    "Primary": "#F58B19",
    "Primary_Secondary_SiteSpec": "green",
    "Primary_Secondary_NHS": "#a34bb0",
}

# legend labels with the mean and standard deviation of each strategy
legend_labels = {
    strategy: f"{strategy} ({combined_df_sw[strategy].mean():.1f} ± {combined_df_sw[strategy].std():.1f} nm)"
    for strategy in combined_df_sw.columns
}
fig = plt.figure()
ax = fig.add_subplot(111)
filament_diameters = sns.histplot(
    data=combined_df_sw.rename(columns=legend_labels),
    palette={legend_labels[k]: v for k, v in custom_palette.items()},
    bins=30, ax=ax)
ax.set_xlabel("Apparent diameter (nm)")
filename = my_experiment_probe3.date_as_string + 'vlab4mic_apparent_diameter_microtubule.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()
filename = my_experiment_probe3.date_as_string + 'vlab4mic_apparent_diameter_microtubule.csv'
combined_df_sw.to_csv(os.path.join(my_experiment.output_directory, filename), index=False)

### Hollow cylinder fit of the 3D localisations
reference_diameter = 25    # nm, microtubule

def orthonormal_basis(direction):
    """Two unit vectors orthogonal to direction and to each other."""
    direction = direction / np.linalg.norm(direction)
    helper = np.array([0.0, 0.0, 1.0]) if abs(direction[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(direction, helper)
    u /= np.linalg.norm(u)
    v = np.cross(direction, u)
    return u, v

def cylinder_coordinates(points, pivot, direction):
    """Position along the axis and (u, v) position across it, for each point."""
    u, v = orthonormal_basis(direction)
    centred = points - pivot
    return centred @ direction, centred @ u, centred @ v

def fit_hollow_cylinder(points, reference_diameter=25, background_radius_factor=3, max_points=20000, random_seed=0):
    """
    Fit a hollow cylinder to 3D coordinates (nm).

    Model: points lie on a cylindrical shell of radius R, blurred with an
    isotropic Gaussian of width sigma (labelling offset, localisation
    precision). Their distance to the axis then follows a Rice distribution
    Rice(R, sigma). Points further than background_radius_factor *
    reference_diameter from the axis are ignored, and a uniform background
    in the cross-section accounts for localisations off the filament. The
    axis is first fitted geometrically (robust least squares of distance to
    axis minus R, starting from the principal axis of the points and the
    reference radius), then axis, R, sigma and the background fraction are
    refined by maximum likelihood.

    Returns a dict with the axis (pivot, direction), radius, sigma,
    background fraction, diameter and the coordinates of every point in the
    cylinder frame (along, u, v, radial).
    """
    points = np.asarray(points, dtype=float)
    rng = np.random.default_rng(random_seed)
    if points.shape[0] > max_points:
        subset = points[rng.choice(points.shape[0], max_points, replace=False)]
    else:
        subset = points
    centroid = subset.mean(axis=0)
    __, eigenvectors = np.linalg.eigh(np.cov((subset - centroid).T))
    direction0 = eigenvectors[:, -1]
    u0, v0 = orthonormal_basis(direction0)
    shift_limit = 2 * reference_diameter

    def axis_from_params(params):
        tilt_u, tilt_v, shift_u, shift_v = params[:4]
        direction = direction0 + tilt_u * u0 + tilt_v * v0
        direction /= np.linalg.norm(direction)
        pivot = centroid + shift_u * u0 + shift_v * v0
        return pivot, direction

    def radial_distances(params, coordinates):
        pivot, direction = axis_from_params(params)
        __, across_u, across_v = cylinder_coordinates(coordinates, pivot, direction)
        return np.hypot(across_u, across_v)

    # geometric fit: distance to the axis equal to the radius
    geometric = optimize.least_squares(
        lambda params: radial_distances(params, subset) - params[4],
        x0=[0, 0, 0, 0, reference_diameter / 2],
        bounds=([-0.5, -0.5, -shift_limit, -shift_limit, 0], [0.5, 0.5, shift_limit, shift_limit, np.inf]),
        loss="soft_l1", f_scale=5,
    )

    # maximum likelihood on the points within max_radial of the geometric
    # axis: distances to the axis follow Rice(R, sigma) plus a uniform
    # background in the cross-section
    max_radial = background_radius_factor * reference_diameter
    near_axis = subset[radial_distances(geometric.x, subset) < max_radial]
    residual_sd = np.std(radial_distances(geometric.x, near_axis) - geometric.x[4])

    def negative_log_likelihood(axis_params, radius, sigma, background):
        radial = radial_distances(axis_params, near_axis)
        shell = stats.rice.pdf(radial, radius / sigma, scale=sigma) / stats.rice.cdf(max_radial, radius / sigma, scale=sigma)
        uniform = 2 * radial / max_radial ** 2
        return -np.sum(np.log((1 - background) * shell + background * uniform + 1e-300))

    shell_bounds = [(0, max_radial), (0.5, max_radial), (0, 0.99)]
    # shell and background with the geometric axis, then all together
    shell_fit = optimize.minimize(
        lambda params: negative_log_likelihood(geometric.x[:4], *params),
        x0=[geometric.x[4], max(residual_sd, 1), 0.05],
        method="L-BFGS-B",
        bounds=shell_bounds,
    )
    likelihood = optimize.minimize(
        lambda params: negative_log_likelihood(params[:4], *params[4:]),
        x0=np.append(geometric.x[:4], shell_fit.x),
        method="L-BFGS-B",
        bounds=[(-0.5, 0.5), (-0.5, 0.5), (-shift_limit, shift_limit), (-shift_limit, shift_limit)] + shell_bounds,
    )
    pivot, direction = axis_from_params(likelihood.x)
    along, across_u, across_v = cylinder_coordinates(points, pivot, direction)
    radius, sigma, background = likelihood.x[4:]
    return dict(
        pivot=pivot,
        direction=direction,
        radius=radius,
        sigma=sigma,
        background=background,
        max_radial=max_radial,
        diameter=2 * radius,
        reference_diameter=reference_diameter,
        geometric_diameter=2 * geometric.x[4],
        along=along,
        across_u=across_u,
        across_v=across_v,
        radial=np.hypot(across_u, across_v),
        n_points=points.shape[0],
        success=bool(likelihood.success),
    )

def plot_hollow_cylinder_fit(fit, axs, color="#F58B19", title="", section_length=100, max_scatter=2000, random_seed=0):
    """
    Draw a hollow cylinder fit on three axes: a section of the fitted
    cylinder (section_length nm long, centred) with its points in 3D, in the
    frame of the cylinder (axs[0] must be 3D); the cross-section of all points
    with the fitted and reference circles; and the distances to the axis with
    the fitted model.
    """
    rng = np.random.default_rng(random_seed)
    radius, sigma, background = fit["radius"], fit["sigma"], fit["background"]
    reference_radius = fit["reference_diameter"] / 2
    limit = min(radius + 4 * sigma, fit["max_radial"])
    angles = np.linspace(0, 2 * np.pi, 100)
    # 3D section, in the frame of the cylinder
    along_centre = np.median(fit["along"])
    in_section = np.flatnonzero((np.abs(fit["along"] - along_centre) <= section_length / 2) & (fit["radial"] <= limit))
    shown = rng.choice(in_section, min(max_scatter, in_section.size), replace=False)
    along_grid, angle_grid = np.meshgrid(np.linspace(-section_length / 2, section_length / 2, 20), angles)
    axs[0].plot_surface(along_grid, radius * np.cos(angle_grid), radius * np.sin(angle_grid),
                        color="grey", alpha=0.25, linewidth=0)
    axs[0].scatter(fit["along"][shown] - along_centre, fit["across_u"][shown], fit["across_v"][shown], s=1, c=color, alpha=0.6, depthshade=False)
    axs[0].set_xlim(-section_length / 2, section_length / 2)
    axs[0].set_ylim(-limit, limit)
    axs[0].set_zlim(-limit, limit)
    axs[0].set_box_aspect((section_length, 2 * limit, 2 * limit))
    axs[0].view_init(elev=20, azim=-35)
    axs[0].set_xlabel("along axis (nm)")
    axs[0].set_ylabel("u (nm)")
    axs[0].set_zlabel("v (nm)")
    axs[0].set_title(title)
    # cross-section
    bins = np.linspace(-limit, limit, 60)
    cross_section, __, __ = np.histogram2d(fit["across_u"], fit["across_v"], bins=[bins, bins])
    axs[1].imshow(cross_section.T, cmap="grey", origin="lower", extent=[-limit, limit, -limit, limit])
    axs[1].plot(radius * np.cos(angles), radius * np.sin(angles), c=color,
                label=f"fit: {2 * radius:.1f} nm")
    axs[1].plot(reference_radius * np.cos(angles), reference_radius * np.sin(angles), "--", c="white",
                label=f"reference: {fit['reference_diameter']:.0f} nm")
    axs[1].set_xlabel("u (nm)")
    axs[1].set_ylabel("v (nm)")
    axs[1].legend(loc="upper right", fontsize="small")
    axs[1].set_title("Cross-section")
    # distances to the axis, with the fitted model
    max_radial = fit["max_radial"]
    radial_range = np.linspace(0, max_radial, 300)
    shell = stats.rice.pdf(radial_range, radius / sigma, scale=sigma) / stats.rice.cdf(max_radial, radius / sigma, scale=sigma)
    model = (1 - background) * shell + background * 2 * radial_range / max_radial ** 2
    axs[2].hist(fit["radial"][fit["radial"] < max_radial], bins=np.linspace(0, max_radial, 60), density=True, color=color, alpha=0.6)
    axs[2].plot(radial_range, model, c="black",
                label=f"fit: R = {radius:.1f} nm, $\\sigma$ = {sigma:.1f} nm, bg = {100 * background:.0f}%")
    axs[2].axvline(radius, c=color)
    axs[2].axvline(reference_radius, ls="--", c="grey", label=f"reference R = {reference_radius:.1f} nm")
    axs[2].set_xlabel("distance to axis (nm)")
    axs[2].set_ylabel("density")
    axs[2].legend(fontsize="small")
    axs[2].set_title("Distance to the cylinder axis")

def replicate_localisations(result, modality="SMLM", channel="ch0"):
    """3D localisations (nm) of all fluorophores of one channel in a replicate."""
    records = result["positions"][modality][channel].values()
    return np.vstack([record["localisations"] for record in records if record["localisations"] is not None])

replicates_per_strategy = {
    "Primary": replicates_1,
    "Primary_Secondary_SiteSpec": replicates_3,
    "Primary_Secondary_NHS": replicates_2,
}
cylinder_fits = dict()
cylinder_rows = []
for strategy, replicates in replicates_per_strategy.items():
    cylinder_fits[strategy] = []
    for result in replicates:
        fit = fit_hollow_cylinder(replicate_localisations(result), reference_diameter=reference_diameter)
        cylinder_fits[strategy].append(fit)
        cylinder_rows.append(dict(
            strategy=strategy,
            replicate=result["replicate"],
            diameter=fit["diameter"],
            sigma=fit["sigma"],
            background=fit["background"],
            success=fit["success"],
        ))
cylinder_df = pd.DataFrame(cylinder_rows)
print(cylinder_df.groupby("strategy")[["diameter", "sigma", "background"]].agg(["mean", "std"]))

# apparent hollow cylinder of the first replicate of each strategy
n_strategies = len(replicates_per_strategy)
fig = plt.figure(figsize=(21, 6 * n_strategies))
for i, strategy in enumerate(replicates_per_strategy):
    axs = [
        fig.add_subplot(n_strategies, 3, 3 * i + 1, projection="3d"),
        fig.add_subplot(n_strategies, 3, 3 * i + 2),
        fig.add_subplot(n_strategies, 3, 3 * i + 3),
    ]
    plot_hollow_cylinder_fit(cylinder_fits[strategy][0], axs, color=custom_palette[strategy], title=strategy)
filename = my_experiment_probe3.date_as_string + 'vlab4mic_hollow_cylinder_fit_microtubule.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.histplot(data=cylinder_df, x="diameter", hue="strategy", palette=custom_palette, bins=30, ax=ax)
ax.axvline(reference_diameter, ls="--", c="grey")
ax.set_xlabel("Hollow cylinder diameter (nm)")
filename = my_experiment_probe3.date_as_string + 'vlab4mic_hollow_cylinder_diameter_microtubule.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

### Apparent continuity

labelling_efficiencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
replicates_per_le = dict()
n_replicates_le = 100
for labeff in labelling_efficiencies:
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
        labelling_efficiency=labeff,
    )
    __, noiseless_, my_experiment_le = experiments.image_vsample(
        primary_probe=primary_1,
        secondary_probe=secondary_1,
        **common
    )
    my_experiment_le.particle.source = dict()
    my_experiment_le.particle.load_source(**copy.deepcopy(particle_generator_data))
    my_experiment_le.set_virtualsample_params(
        sample_dimensions=sample_dimensions,
        update_mode = False 
    )
    my_experiment_le.particle.generate_instance()
    my_experiment_le.build(modules=["coordinate_field", "imager"])
    replicates_le = my_experiment_le.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le[labeff] = copy.copy(replicates_le)

filament_measures_per_le = dict()
smooth_sigma = 0.0

filaments_params_le = dict(
    pixel_size=2, 
    max_length_nm=500, 
    modality="SMLM", 
    channel="ch0", 
    fill_holes=True,
    row_binning=7,
    min_snr=2,
    mode="total",
    smooth_sigma=smooth_sigma,
    threshold_method = "otsu"
)

for i in labelling_efficiencies:
    df_, filament_metrics_list_, masks_, blurred_, tables_ = measure_filaments(
    replicates=replicates_per_le[i], 
    **filaments_params_le
    )
    filament_measures_per_le[i] = dict()
    filament_measures_per_le[i]["dataframe"] = df_
    filament_measures_per_le[i]["filament_metrics_list"] = filament_metrics_list_
    filament_measures_per_le[i]["masks"] = masks_
    filament_measures_per_le[i]["blurred"] = blurred_
    filament_measures_per_le[i]["tables"] = tables_

dfs = dict() 
for lef in labelling_efficiencies:
    dfs[lef] = filament_measures_per_le[lef]["dataframe"]

merged = pd.concat(
    [d.assign(source=name) for name, d in dfs.items()],
    ignore_index=True,
)
filename = my_experiment_le.date_as_string + 'vlab4mic_filament_measures_per_LE.csv'
merged.to_csv(os.path.join(my_experiment_le.output_directory, filename), index=False)
fig = plt.figure()
ax = fig.add_subplot(111)
sns.swarmplot(data=merged, x="source", y="mean_eccentricity", ax=ax, s=1)
filename = my_experiment_le.date_as_string + 'vlab4mic_mean_eccentricities_per_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
plt.xlabel("Labelling efficiency")
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()

n_le = len(labelling_efficiencies)
fig, ax = plt.subplots(1, n_le, figsize=(30,10))
for i, lef in enumerate(labelling_efficiencies):
    ax[i].imshow(filament_measures_per_le[lef]["blurred"][0], cmap = "grey")
    ax[i].set_title(f"LE = {lef}")
filename = my_experiment_le.date_as_string + 'vlab4mic_microtubules_across_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
fig.savefig(filename2, dpi=300, bbox_inches='tight')
plt.close()