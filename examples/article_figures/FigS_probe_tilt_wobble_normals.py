"""
Probe tilt and wobble on epitopes on the z = 0 plane, normals along +z.

Tilt and wobble figure (square): along x the tilt increases; along y the
wobble (binding.wobble_range.theta) increases. Probes are placed with
indirect_labelling, which tilts each axis (tilt_direction) and then wobbles
it within a cone about the tilted axis (add_wobble).

Normals figure: normals estimated by local plane fitting and by scaling on
a flat lattice, a sphere, a cylinder and the NPC.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from vlab4mic.utils.transform.cif_builder import indirect_labelling
from vlab4mic.utils.transform.normals import tilt_direction

plt.rcParams.update({
    "font.size": 18,
    "axes.titlesize": 22,
    "axes.labelsize": 20,
    "legend.fontsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
})

random_seed = 24
tilts = [0, 15, 30, 45, 60]  # degrees, x axis
wobbles = [0, 10, 20, 30]  # degrees, y axis
normal = np.array([0, 0, 1.0])
spacing = 40
probe_length = 12
n_probes = 200  # probes drawn per site
n_rim_points = 100  # tilted directions sampled to trace the tilt cone

normal_colour = "#222222"
probe_colour = "#377eb8"
cone_colour = "red"


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


def draw_site(ax, epitope, tips, tilt, rng):
    ax.scatter(*tips.T, color=probe_colour, s=10, depthshade=False, clip_on=False)
    if tilt:
        ax.plot(*cone_rim(epitope, tilt, rng).T, color=cone_colour, lw=1, clip_on=False)
    normal_tip = epitope + 1.3 * probe_length * normal
    #ax.plot(*np.stack([epitope, normal_tip]).T, color=normal_colour, lw=1.5)
    ax.scatter(*epitope, color=normal_colour, s=25, depthshade=False, clip_on=False)


def add_legend(ax):
    ax.plot([], [], color=normal_colour, marker="o", ls="none", ms=6, label="Epitope")
    ax.plot([], [], color=probe_colour, marker="o", ls="none", ms=4, label="Probe emitter")
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
    ax.set_xlabel("Probe tilt from normal", labelpad=55)
    ax.set_zticks([])
    fill_box(ax, zoom=0.8)
    ax.view_init(elev=30, azim=-60)


np.random.seed(random_seed)
rng = np.random.default_rng(random_seed)

# tilt combined with wobble about the tilted axis, square figure
fig_tw = plt.figure(figsize=[10, 10])
fig_tw.subplots_adjust(left=0, right=1, bottom=0, top=0.93)
ax_tw = fig_tw.add_subplot(projection="3d")
for i, wobble in enumerate(wobbles):
    for j, tilt in enumerate(tilts):
        epitope = np.array([j * spacing, i * spacing, 0.0])
        draw_site(ax_tw, epitope, place_probes(epitope, tilt, wobble), tilt, rng)
style_tilt_axis(ax_tw)
ax_tw.set_yticks([i * spacing for i in range(len(wobbles))])
ax_tw.set_yticklabels([f"{w}°" for w in wobbles])
ax_tw.set_ylabel("Wobble about tilted axis", labelpad=25)
ax_tw.set_title("Tilt combined with wobble")
add_legend(ax_tw)

# normals figure
fig_n = plt.figure(figsize=[20, 10])
fig_n.subplots_adjust(left=0.09, right=1, top=0.95, bottom=0, wspace=0.02, hspace=0.05)
axs_n = fig_n.subplots(
    2, 4, subplot_kw={"projection": "3d"},
    gridspec_kw={"width_ratios": [1, 1, 0.65, 1]},
)
for row, method in enumerate(["By local plane", "By scaling"]):
    axs_n[row, 0].text2D(-0.2, 0.5, method, transform=axs_n[row, 0].transAxes,
                         rotation=90, va="center", ha="center", fontsize=24)


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
show(points_cyl[body_cyl], normals_cyl_local[body_cyl], "Cylinder", axs_n[0, 2], length=10, arrow_length_ratio=0, view_init=[0,-60,0], zoom=0.9)

normals_cyl_scaling = normals_by_scaling(points_cyl)
body_cyl = np.abs(normals_cyl_scaling[:, 2]) < 80
show(points_cyl[body_cyl], normals_cyl_scaling[body_cyl], "Cylinder", axs_n[1, 2], length=3, arrow_length_ratio=0, view_init=[0,-60,0], zoom=0.9)

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
my_experiment.structure.show_target_labels(with_normals=True, with_assembly_atoms=True, assembly_fraction=0.001, view_init=[20,0,0], axesoff=True, show_axis=False, reference_point=False, axis_object=axs_n[1, 3])
axs_n[1, 3].set_title("NPC (7R5K)")
fill_box(axs_n[1, 3], zoom=1.5)
# by plane
n_neighbours_npc = 2
my_experiment.set_structure_normal_params(mode="local_plane", n_neighbours=n_neighbours_npc)
my_experiment.structure.assign_normals2targets()
my_experiment.structure.show_target_labels(with_normals=True, with_assembly_atoms=True, assembly_fraction=0.001, view_init=[20,0,0], axesoff=True, show_axis=False, reference_point=False, axis_object=axs_n[0, 3])
axs_n[0, 3].set_title("NPC (7R5K)")
fill_box(axs_n[0, 3], zoom=1.5)

filename = my_experiment.date_as_string + "figS_probe_tilt_wobble_normals"
output_dir = os.path.join(my_experiment.output_directory, filename)
os.makedirs(output_dir, exist_ok=True)
filename2 = os.path.join(output_dir, filename)
fig_tw.savefig(filename2 + "_tilt_wobble.png", dpi=300)
plt.close(fig_tw)
fig_n.savefig(filename2 + "_normals.png", dpi=300, bbox_inches="tight")
plt.close(fig_n)

legend = (
    "Tilt and wobble. Probe orientation on epitopes on the z = 0 plane, with normals along +z. "
    f"The probe axis is tilted from the normal by {', '.join(str(t) for t in tilts)}° (along x) with a "
    "random azimuth, so before wobble it lies on the tilt cone (red circle), and is then wobbled about "
    f"the tilted axis by {', '.join(str(w) for w in wobbles)}° (along y). Each site shows the emitters "
    f"(blue) of {n_probes} probes with length {probe_length}; black dots mark the epitopes.\n\n"
    "Normals. Normals estimated for target sites by fitting a local plane to neighbouring sites (top row) "
    "and by scaling about the centroid of the sites (bottom row), on a flat lattice, a sphere, "
    "a cylinder (side view) and the nuclear pore complex (PDB 7R5K, Nup96 C-terminal "
    f"sequence ELAVGSL; local plane fitted to {n_neighbours_npc} neighbouring sites). "
    "Scaling gives the surface normal only for convex shapes centred on the centroid: on the flat "
    "lattice its normals lie in the plane, and on the cylinder they bend towards the ends. "
    "Fitting a local plane gives the surface normal on each of these surfaces.\n"
)
with open(filename2 + "_legend.txt", "w") as f:
    f.write(legend)