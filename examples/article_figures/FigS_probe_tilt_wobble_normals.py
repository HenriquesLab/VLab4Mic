"""
Probe tilt and wobble on epitopes on the z = 0 plane, normals along +z.

Panel a: the tilt step alone. At each site the probe axis is tilted from
the normal by the tilt angle with tilt_direction (random azimuth), so the
axes lie on the rim of the tilt cone.

Panel b: the combined result. Along x the tilt increases; along y the
wobble (binding.wobble_range.theta) increases. Probes are placed with
indirect_labelling, which tilts each axis (tilt_direction) and then wobbles
it within a cone about the tilted axis (add_wobble).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from vlab4mic.utils.transform.cif_builder import indirect_labelling
from vlab4mic.utils.transform.normals import tilt_direction

random_seed = 24
tilts = [0, 15, 30, 45, 60]  # degrees, x axis
wobbles = [0, 10, 20, 30]  # degrees, y axis of panel b
normal = np.array([0, 0, 1.0])
spacing = 40
probe_length = 12
n_tilted_axes = 8  # tilted axes drawn per site in panel a
n_probes = 200  # probes drawn per site in panel b
n_rim_points = 100  # tilted directions sampled to trace the tilt cone

normal_colour = "#222222"
probe_colour = "#377eb8"
cone_colour = "#9ecae1"


def place_probes(epitope, tilt, wobble):
    # labelling entity: anchor, axis point, one emitter along the axis
    entity = np.array([[0, 0, 0], [0, 0, 1], [0, 0, probe_length]], dtype=float)
    label_data = {
        "emitters": entity,
        "labelling_efficiency": 1,
        "minimal_distance": 0,
        "conjugation_sites": {"DoL": None},
        "binding": {"tilt": tilt, "wobble_range": {"theta": wobble or None}},
    }
    # n_probes copies of the same site, each with its own tilt azimuth and wobble
    sites = {
        "coordinates": np.tile(epitope, (n_probes, 1)),
        "normals": np.tile(normal, (n_probes, 1)),
    }
    emitters, _, _ = indirect_labelling(sites, label_data)
    return emitters


def cone_rim(epitope, tilt, rng):
    # sample tilted directions and order them by azimuth to close the circle
    dirs = np.array([tilt_direction(normal, tilt, rng=rng) for _ in range(n_rim_points)])
    dirs = dirs[np.argsort(np.arctan2(dirs[:, 1], dirs[:, 0]))]
    dirs = np.vstack([dirs, dirs[:1]])
    return epitope + probe_length * dirs


def draw_site(ax, epitope, tips, tilt, rng, add_tips=True):
    if add_tips:
        for tip in tips:
            ax.plot(*np.stack([epitope, tip]).T, color=probe_colour, lw=0.8, alpha=0.6, clip_on=False)
    ax.scatter(*tips.T, color=probe_colour, s=10, depthshade=False, clip_on=False)
    if tilt:
        ax.plot(*cone_rim(epitope, tilt, rng).T, color="Red", lw=1, clip_on=False)
    normal_tip = epitope + 1.3 * probe_length * normal
    #ax.plot(*np.stack([epitope, normal_tip]).T, color=normal_colour, lw=1.5)
    ax.scatter(*epitope, color=normal_colour, s=25, depthshade=False, clip_on=False)


def add_legend(ax):
    ax.plot([], [], color=normal_colour, lw=1.5, label="Epitope normal")
    ax.plot([], [], color=probe_colour, marker="o", ms=3, label="Probe axis and emitter")
    ax.plot([], [], color=cone_colour, label="Tilt cone")
    ax.legend(loc="upper left")


def fill_box(ax, zoom=1.3):
    # equal aspect from the data limits, with the box enlarged inside its axes
    ax.set_box_aspect(
        [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")], zoom=zoom
    )


def style_tilt_axis(ax):
    ax.set_xticks([j * spacing for j in range(len(tilts))])
    ax.set_xticklabels([f"{t}°" for t in tilts])
    ax.set_xlabel("Probe tilt from normal")
    ax.set_zticks([])
    fill_box(ax, zoom=1.3)
    ax.view_init(elev=30, azim=-60)


np.random.seed(random_seed)
rng = np.random.default_rng(random_seed)
fig = plt.figure(figsize=[20, 14])
top, bottom = fig.subfigures(2, 1, height_ratios=[0.55, 1])
top.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.05)
bottom.subplots_adjust(left=0.07, right=1, top=0.97, bottom=0, wspace=0.02, hspace=0.05)
ax_a = top.add_subplot(1, 2, 1, projection="3d")
ax_b = top.add_subplot(1, 2, 2, projection="3d")
axs_n = bottom.subplots(
    2, 5, subplot_kw={"projection": "3d"},
    gridspec_kw={"width_ratios": [1, 1, 0.65, 0.65, 1]},
)
for row, method in enumerate(["By local plane", "By scaling"]):
    axs_n[row, 0].text2D(-0.13, 0.5, method, transform=axs_n[row, 0].transAxes,
                         rotation=90, va="center", ha="center", fontsize=14)

# panel a: tilt only
for j, tilt in enumerate(tilts):
    epitope = np.array([j * spacing, 0.0, 0.0])
    axes = np.array([tilt_direction(normal, tilt, rng=rng) for _ in range(n_tilted_axes)])
    draw_site(ax_a, epitope, epitope + probe_length * axes, tilt, rng)
style_tilt_axis(ax_a)
ax_a.set_yticks([])
ax_a.set_title("Tilting the probe axis from the normal")
add_legend(ax_a)

# panel b: tilt combined with wobble about the tilted axis
for i, wobble in enumerate(wobbles):
    for j, tilt in enumerate(tilts):
        epitope = np.array([j * spacing, i * spacing, 0.0])
        draw_site(ax_b, epitope, place_probes(epitope, tilt, wobble), tilt, rng, add_tips=False)
style_tilt_axis(ax_b)
ax_b.set_yticks([i * spacing for i in range(len(wobbles))])
ax_b.set_yticklabels([f"{w}°" for w in wobbles])
ax_b.set_ylabel("Wobble about tilted axis")
ax_b.set_title("Tilt combined with wobble")
add_legend(ax_b)


import numpy as np

from vlab4mic.utils.transform.normals import (
    normals_by_local_plane,
    normals_by_scaling,
)


def flat_lattice(n_side=15, spacing=10.0, z=50.0):
    g = np.stack(
        np.meshgrid(np.arange(n_side) * spacing, np.arange(n_side) * spacing), -1
    ).reshape(-1, 2)
    return np.c_[g, np.full(len(g), z)]


def fibonacci_sphere(n=400, radius=50.0):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return radius * np.c_[
        np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    ]


def cylinder(radius=12.5, length=200.0, n_around=26, n_along=40):
    angles = np.linspace(0, 2 * np.pi, n_around, endpoint=False)
    zs = np.linspace(-length / 2, length / 2, n_along)
    a, z = np.meshgrid(angles, zs)
    return np.c_[radius * np.cos(a.ravel()), radius * np.sin(a.ravel()), z.ravel()]


def unit(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)

import matplotlib.pyplot as plt

def show(points, normals, title, ax, length=10, arrow_length_ratio = 0.3, view_init = [30,-60,0], zoom=1.15):
    ax.scatter(*points.T, s=4, c="grey", clip_on=False)
    ax.quiver(*points.T, *normals.T, length=length, arrow_length_ratio=arrow_length_ratio, clip_on=False)
    ax.set_title(title)
    # limits that include the arrow tips, so the box holds the whole panel
    ends = points + length * normals
    lo, hi = np.minimum(points, ends).min(axis=0), np.maximum(points, ends).max(axis=0)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], hi[2]))
    fill_box(ax, zoom=zoom)
    ax.set_axis_off()
    ax.set_anchor("N")
    ax.view_init(*view_init)


points_flat = flat_lattice()
up = normals_by_local_plane(points_flat, reference_vector=np.array([0, 0, 1]))
down = normals_by_local_plane(points_flat, reference_vector=np.array([0, 0, -1]))
show(points_flat, up, "Flat lattice (reference +z)", axs_n[0, 0], length=10, zoom=1.25)

points_flat2 = flat_lattice()
mask2 = np.any(points_flat2[:, :2] != 70, axis=1)
normals2 = unit(normals_by_scaling(points_flat2)[mask2])
show(points_flat2[mask2], normals2, "Flat lattice", axs_n[1, 0], arrow_length_ratio=0.0, zoom=1.25)

points_fib = fibonacci_sphere()

normals_fib_scaling = normals_by_scaling(points_fib)
show(points_fib, normals_fib_scaling, "Sphere", axs_n[1, 1], length=5, arrow_length_ratio=0.0, zoom=1.4)

normals_fib_localplane = normals_by_local_plane(points_fib)
show(points_fib, normals_fib_localplane, "Sphere", axs_n[0, 1], length=10, arrow_length_ratio=0.0, zoom=1.4)
points_cyl = cylinder(n_along=30)

normals_cyl_local = normals_by_local_plane(points_cyl, n_neighbours=8)
body_cyl = np.abs(normals_cyl_local[:, 2]) < 80
show(points_cyl[body_cyl], normals_cyl_local[body_cyl], "Cylinder", axs_n[0, 2], length=10, arrow_length_ratio=0, zoom=0.95)
show(points_cyl[body_cyl], normals_cyl_local[body_cyl], "Cylinder (side view)", axs_n[0, 3], length=10, arrow_length_ratio=0, view_init=[0,-60,0], zoom=0.9)

normals_cyl_scaling = normals_by_scaling(points_cyl)
body_cyl = np.abs(normals_cyl_scaling[:, 2]) < 80
show(points_cyl[body_cyl], normals_cyl_scaling[body_cyl], "Cylinder", axs_n[1, 2], length=10, arrow_length_ratio=0, zoom=0.95)
show(points_cyl[body_cyl], normals_cyl_scaling[body_cyl], "Cylinder (side view)", axs_n[1, 3], length=3, arrow_length_ratio=0, view_init=[0,-60,0], zoom=0.9)

from vlab4mic import experiments
random_seed = 24
structure = "7R5K"
probe_template = "GFP_w_nanobody"
modalities = ["Widefield", "Confocal", "AiryScan","STED", "SMLM"]
# Run simulation
# Prepare experiment object
#default
_1, _2, my_experiment = experiments.image_vsample(
    structure=structure,
    probe_template=probe_template,
    probe_DoL=2,
    labelling_efficiency=0.4,
    probe_target_type="Sequence",
    probe_target_value="ELAVGSL", # Nup96 C-terminal
    multimodal=modalities,
    random_rotation=True,
    clear_experiment=True,
    run_simulation=False,
    random_rotations=False,
    random_seed=random_seed
)
my_experiment.structure.show_target_labels(with_normals=True, with_assembly_atoms=True, assembly_fraction=0.001, view_init=[20,0,0], axesoff=True, show_axis=False, reference_point=False, axis_object=axs_n[1, 4])
axs_n[1, 4].set_title("NPC (7R5K)")
fill_box(axs_n[1, 4], zoom=1.5)
# by plane
n_neighbours_npc = 2
my_experiment.set_structure_normal_params(mode="local_plane", n_neighbours=n_neighbours_npc)
my_experiment.structure.assign_normals2targets()
my_experiment.structure.show_target_labels(with_normals=True, with_assembly_atoms=True, assembly_fraction=0.001, view_init=[20,0,0], axesoff=True, show_axis=False, reference_point=False, axis_object=axs_n[0, 4])
axs_n[0, 4].set_title("NPC (7R5K)")
fill_box(axs_n[0, 4], zoom=1.5)

filename = my_experiment.date_as_string + "figS_probe_tilt_wobble_normals"
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2 + ".png", dpi=300, bbox_inches="tight")
plt.close()

legend = (
    "Top. Probe orientation on epitopes on the z = 0 plane, with normals along +z. "
    f"(a) Tilt only: at each site, {n_tilted_axes} probe axes are tilted from the normal by "
    f"{', '.join(str(t) for t in tilts)}° with a random azimuth, so they lie on the tilt cone (red circle). "
    "(b) Tilt combined with wobble: the tilt increases along x and the wobble about the tilted axis "
    f"({', '.join(str(w) for w in wobbles)}°) increases along y. Each site shows the emitters (blue) "
    f"of {n_probes} probes with length {probe_length}; black dots mark the epitopes.\n\n"
    "Bottom. Normals estimated for target sites by fitting a local plane to neighbouring sites (top row) "
    "and by scaling about the centroid of the sites (bottom row), on a flat lattice, a sphere, "
    "a cylinder (oblique and side views) and the nuclear pore complex (PDB 7R5K, Nup96 C-terminal "
    f"sequence ELAVGSL; local plane fitted to {n_neighbours_npc} neighbouring sites). "
    "Scaling gives the surface normal only for convex shapes centred on the centroid: on the flat "
    "lattice its normals lie in the plane, and on the cylinder they bend towards the ends. "
    "Fitting a local plane gives the surface normal on each of these surfaces.\n"
)
with open(filename2 + "_legend.txt", "w") as f:
    f.write(legend)