import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from vlab4mic import experiments
from scipy import optimize
import pandas as pd
import seaborn as sns


def fit_circle_fixed_radius(points, radius):
    """
    Least-squares centre of a circle of known radius through 2-D points.

    The starting point is an algebraic (Kasa) circle fit, so the fit does not
    drift away from the ring when its radius differs from `radius`.
    """
    xy = np.asarray(points, dtype=float)[:, :2]
    if len(xy) >= 3:
        A = np.c_[2 * xy, np.ones(len(xy))]
        start = np.linalg.lstsq(A, (xy**2).sum(axis=1), rcond=None)[0][:2]
    else:
        start = xy.mean(axis=0)
    fit = least_squares(lambda c: np.linalg.norm(xy - c, axis=1) - radius, start)
    return fit.x


def count_npc_corners(points, radius, r_min=30, r_max=70, n_corners=8, n_steps=1000, min_per_section=0):
    """
    Number of detected corners of a nuclear pore (Thevathasan et al. 2019).

    1. A circle of fixed radius is fitted and coordinates converted to (phi_i, r_i).
    2. Localisations with r_i < r_min or r_i > r_max are discarded as background.
    3. phi_rot minimises the summed squared angular distance between phi_rot and
       phi_i modulo 2*pi/n_corners (wrapped, so both sides of an edge count as close).
    4. Corners = segments of width 2*pi/n_corners, starting at phi_rot - pi/n_corners,
       that contain at least one localisation.

    Parameters
    ----------
    points : numpy.ndarray
        Nx2 or Nx3 localisations of one pore, in nm.
    radius : float
        Fixed circle radius in nm (mean radius of the data set).
    r_min, r_max : float
        Radial window kept around the centre, in nm.
    n_corners : int
        Rotational symmetry (8 for the NPC).
    n_steps : int
        Candidate rotations tried within one segment.

    Returns
    -------
    dict with
        n_corners      : int, number of occupied segments
        points         : input coordinates
        centre         : fitted centre (x, y), None if no points
        radius, r_min, r_max, symmetry : parameters used
        r, phi         : polar coordinates of all localisations
        kept           : bool mask, localisations inside [r_min, r_max]
        phi_rot        : ring rotation in radians (None if nothing kept)
        aligned        : Nx2 coords of all localisations, centred and rotated by
                         -phi_rot (segment k is centred on angle k*2*pi/n_corners)
        segment        : segment index of each kept localisation
        segment_counts : localisations per segment
        segment_edges  : n_corners+1 segment edges, radians, original frame
    """
    if points is None:
        result = dict(
            n_corners=0, points=0, centre=None, radius=0, r_min=0,
            r_max=0, symmetry=0, r=np.empty(0), phi=np.empty(0),
            kept=np.zeros(1, bool), phi_rot=None, aligned=np.empty((0, 2)),
            segment=np.empty(0, int), segment_counts=np.zeros(0, int),
            segment_edges=None,
        )
        return result
    else: 
        points = np.asarray(points, dtype=float).reshape(-1, np.shape(points)[-1] if np.ndim(points) == 2 else 3)
        xy = points[:, :2]
        period = 2 * np.pi / n_corners
        result = dict(
            n_corners=0, points=points, centre=None, radius=radius, r_min=r_min,
            r_max=r_max, symmetry=n_corners, r=np.empty(0), phi=np.empty(0),
            kept=np.zeros(len(xy), bool), phi_rot=None, aligned=np.empty((0, 2)),
            segment=np.empty(0, int), segment_counts=np.zeros(n_corners, int),
            segment_edges=None,
        )
        if len(xy) == 0:
            return result

        # 1. fixed-radius circle fit, polar coordinates
        centre = fit_circle_fixed_radius(xy, radius)
        dx, dy = (xy - centre).T
        r = np.hypot(dx, dy)
        phi = np.arctan2(dy, dx)
        # 2. radial filter
        kept = (r >= r_min) & (r <= r_max)
        result.update(centre=centre, r=r, phi=phi, kept=kept)
        if not kept.any():
            return result

        # 3. rotation of the ring
        phi_k = phi[kept]
        candidates = np.linspace(0, period, n_steps, endpoint=False)
        diff = (phi_k[None, :] - candidates[:, None] + period / 2) % period - period / 2
        phi_rot = candidates[np.argmin((diff**2).sum(axis=1))]

        # 4. occupied segments
        start = phi_rot - period / 2
        segment = np.floor(((phi_k - start) % (2 * np.pi)) / period).astype(int) % n_corners
        counts = np.bincount(segment, minlength=n_corners)

        c, s = np.cos(-phi_rot), np.sin(-phi_rot)
        aligned = np.c_[dx, dy] @ np.array([[c, s], [-s, c]])
        result.update(
            #n_corners=int((counts > 0).sum()), phi_rot=float(phi_rot),
            n_corners=int((counts >= min_per_section).sum()), phi_rot=float(phi_rot),
            aligned=aligned, segment=segment, segment_counts=counts,
            segment_edges=start + period * np.arange(n_corners + 1),
        )
        return result


def corner_histogram(particles, radius, n_corners=8, **kwargs):
    """
    Histogram of detected corners over many pores.

    Returns
    -------
    histogram : numpy.ndarray
        Number of pores with 0, 1, ..., n_corners corners.
    results : list of dict
        count_npc_corners output for each pore.
    """
    results = [count_npc_corners(p, radius, n_corners=n_corners, **kwargs) for p in particles]
    counts = np.array([res["n_corners"] for res in results], dtype=int)
    return np.bincount(counts, minlength=n_corners + 1), results


def plot_npc_corners(result, axes=None):
    """
    (a) input, fixed-radius circle and radial window (excluded in grey);
    (b) kept localisations centred and rotated by -phi_rot;
    (c) kept localisations coloured by segment, segment edges, detected corners
        (filled circle = occupied, light circle = empty; numbers = counts).
    """
    if axes is None:
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    else:
        fig = axes[0].figure
    xy = result["points"][:, :2]
    kept = result["kept"]
    n = result["symmetry"]
    period = 2 * np.pi / n
    t = np.linspace(0, 2 * np.pi, 200)
    circle = np.c_[np.cos(t), np.sin(t)]

    ax = axes[0]
    ax.scatter(*xy[~kept].T, s=6, c="0.7", label="excluded")
    ax.scatter(*xy[kept].T, s=6, c="C0", label="kept")
    if result["centre"] is not None:
        c = result["centre"]
        ax.plot(*(c + result["radius"] * circle).T, "k-", lw=1, label="fixed-radius fit")
        for rr in (result["r_min"], result["r_max"]):
            ax.plot(*(c + rr * circle).T, "k:", lw=0.8)
        ax.plot(*c, "k+", ms=10)
    ax.set_title("(a) input and circle fit")
    ax.legend(fontsize=7, loc="upper right")

    al = result["aligned"]
    r_edge = result["r_max"]
    ax = axes[1]
    if len(al):
        ax.scatter(*al[kept].T, s=6, c="C0")
        ax.plot(*(result["radius"] * circle).T, "k-", lw=1)
        ax.set_title(f"(b) aligned, phi_rot = {np.degrees(result['phi_rot']):.1f} deg")
    else:
        ax.set_title("(b) aligned (no localisations kept)")

    ax = axes[2]
    if len(al):
        ax.scatter(*al[kept].T, s=8, c=result["segment"], cmap="tab10", vmin=-0.5, vmax=9.5)
        for k in range(n):
            edge = (k + 0.5) * period
            ax.plot([0, r_edge * np.cos(edge)], [0, r_edge * np.sin(edge)], "k-", lw=0.5)
            mid = k * period
            occupied = result["segment_counts"][k] > 0
            ax.plot(result["radius"] * np.cos(mid), result["radius"] * np.sin(mid),
                    marker="o", ms=12, ls="none",
                    mfc=(0, 0, 0, 0.15) if occupied else "none",
                    mec="k" if occupied else "0.6")
            ax.text(1.15 * r_edge * np.cos(mid), 1.15 * r_edge * np.sin(mid),
                    str(result["segment_counts"][k]), ha="center", va="center", fontsize=8)
    ax.set_title(f"(c) segments: {result['n_corners']} corners")

    lim = 1.3 * r_edge
    for ax in axes[1:]:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
    fig.tight_layout()
    return fig

def split_by_midpointZ(positions, split=True):
    if positions is None:
        return None
    else:
        if split:
            if positions.shape[1] == 3:
                midpointZ = positions[:,2].mean()
                row_condition = positions[:,2] > midpointZ
                locs = np.compress(row_condition, positions, axis=0)
                return locs
            else:
                print("Needs 3 dimensions")
                return None
        else:
            return positions

def fit_circle(xcoord, ycoord):
    if xcoord is None or ycoord is None:
        return None, None, None, None, None
    else:
        x_m = np.mean(xcoord)
        y_m = np.mean(ycoord)
        center_estimate = x_m, y_m
        
        def calc_R(xc, yc, xdata, ydata):
            """ calculate the distance of each 2D points from the center (xc, yc) """
            return np.sqrt((xdata-xc)**2 + (ydata-yc)**2)

        def f_2(starting_center, xdata, ydata):
            """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
            Ri = calc_R(*starting_center, xdata, ydata)
            return Ri - Ri.mean()
        center_2, ier = optimize.leastsq(func=f_2, x0=[x_m, y_m], args=(xcoord, ycoord))

        xc_2, yc_2 = center_2
        Ri_2       = calc_R(*center_2, xcoord, ycoord)
        R_2        = Ri_2.mean()
        residu_2   = sum((Ri_2 - R_2)**2)
        return xc_2, yc_2, Ri_2, R_2, residu_2


def sector_counting(replicates,
                    modality="SMLM",
                    channel="ch0",
                    split_z = False,
                    n_sites=8,
                    mean_radius_NPC_nm=52,
                    min_locs_per_segment=4,
                    min_locs_to_circlefit=10,
                    max_radius=80,
                    ):
    rows = []
    images_list = []
    for result in replicates:
        record = list(result["positions"][modality][channel].values())[0]
        locs = split_by_midpointZ(record["localisations"], split=split_z)
        image = result["images"][modality][channel]
        n_locs = 0 if locs is None else len(locs)
        npc_corner_counts_result = count_npc_corners(
            points=locs,
            radius=mean_radius_NPC_nm,
            n_corners=n_sites,
            min_per_section=min_locs_per_segment)
        if n_locs >= min_locs_to_circlefit:
            xc_2, yc_2, Ri_2, radius_fit, residu_2 = fit_circle(locs[:,0],locs[:,1])
            if radius_fit >= max_radius:
                radius_fit = None
        else:
            xc_2 = yc_2 = Ri_2 = radius_fit = residu_2 = None    
        rows.append(dict(
            replicate=result["replicate"],
            n_localisations=n_locs,
            radius_fit = radius_fit,
            #npc_corner_counting_result = npc_corner_counts_result,
            n_corners=npc_corner_counts_result["n_corners"],
            points = npc_corner_counts_result["points"],
            centre = npc_corner_counts_result["centre"],
            radius= npc_corner_counts_result["radius"],
            r_min= npc_corner_counts_result["r_min"],
            r_max= npc_corner_counts_result["r_max"],
            symmetry= npc_corner_counts_result["symmetry"],
            r= npc_corner_counts_result["r"],
            phi= npc_corner_counts_result["phi"],
            kept= npc_corner_counts_result["kept"],
            phi_rot= npc_corner_counts_result["phi_rot"],
            aligned= npc_corner_counts_result["aligned"],
            segment= npc_corner_counts_result["segment"],
            segment_counts= npc_corner_counts_result["segment_counts"],
            segment_edges= npc_corner_counts_result["segment_edges"],
        ))
        images_list.append(image)

    results = pd.DataFrame(rows)
    rows_results = rows
    return results, images_list, rows


###########

RANDOM_SEED = 7
PRECISION_NM = 10
MIN_SEPARATION_NM = 2 * PRECISION_NM  # peak separation and DBSCAN radius
DBSCAN_MIN_SAMPLES = 3
PEAK_THRESHOLD_REL = 0.2  # default of count_resolved_sites
N_SEGMENTS = 8  # maximum number of segments of the ring in 2D
nlocalisations=10
EXAMPLE_REPLICATE = 0
RENDERING_KERNEL = 3
depth_of_field_nm = 20
depth_of_field_nm = 1000

labelling_efficiencies = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
replicates_per_le = dict()
parameters_le = dict()
n_replicates_le = 100
probe_DoL = None
import copy

for labeff in labelling_efficiencies:
    primary_le = dict(
        probe_template="GFP_w_nanobody",
        probe_name="GFP_w_nanobody",
        probe_target_type="Sequence",
        probe_target_value="ELAVGSL",
        probe_distance_to_epitope = 0,
        labelling_efficiency=labeff,
        probe_DoL=probe_DoL,
    )
    # GFP_w_nanobody has no default target: Nup96 C-terminus
    _, _, my_experiment_le = experiments.image_vsample(
        structure="7R5K",
        **primary_le,
        random_rotations=True,
        multimodal=["SMLM"],
        number_of_particles=1,
        sample_dimensions=[500, 500, 1],
        sample_inital_orientation=[0, 0, 1],  # pore seen from above
        run_simulation=False,
        clear_experiment=True,
        random_seed=RANDOM_SEED,
    )
    my_experiment_le.update_modality(
        modality_name="SMLM",
        simulate_localistations=True,
        nlocalisations=nlocalisations,
        depth_of_field_nm=depth_of_field_nm,
        lateral_precision=PRECISION_NM,
        axial_precision=PRECISION_NM,
        precision_model="fixed",
        lateral_resolution_nm=RENDERING_KERNEL,
    )
    
    replicates_le = my_experiment_le.run_replicates(n_replicates=n_replicates_le, modality="SMLM")
    replicates_per_le[labeff] = copy.copy(replicates_le)
    parameters_le[labeff] = my_experiment_le.get_parameters()



####### 
counting_measures_per_le = dict()
split_z = False
min_locs_per_segment = 4
sectorcounting_params_le = dict(
    modality="SMLM",
    channel="ch0",
    split_z = split_z,
    n_sites=8,
    mean_radius_NPC_nm=52,
    min_locs_per_segment=min_locs_per_segment,
    min_locs_to_circlefit=nlocalisations
)

for i in labelling_efficiencies:
    results, images_list, rows = sector_counting(
    replicates=replicates_per_le[i], 
    **sectorcounting_params_le
    )
    counting_measures_per_le[i] = dict()
    counting_measures_per_le[i]["dataframe"] = results
    counting_measures_per_le[i]["images"] = images_list
    counting_measures_per_le[i]["rows"] = rows

dfs = dict() 
for lef in labelling_efficiencies:
    dfs[lef] = counting_measures_per_le[lef]["dataframe"]

merged = pd.concat(
    [d.assign(source=name) for name, d in dfs.items()],
    ignore_index=True,
)
######## plot 1
sns.boxplot(data=merged, x="source", y="n_corners")
plt.xlabel("Labelling efficiency")
plt.ylabel("Number of Corners detected")
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
plt.title("Corners detected on NPC simuations" + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))


#### plot 2
n_corners = 8
lab_eff_reference = 0.3
counts = counting_measures_per_le[lab_eff_reference]["dataframe"].n_corners
hist =  np.bincount(counts, minlength=n_corners + 1)
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar(np.arange(9), hist, label="detected")
ax.set_xlabel("corners per NPC")
ax.set_ylabel("NPCs")
ax.legend()
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
ax.set_title("Corners detected on NPC simuations. LE = " + str(lab_eff_reference) + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))

### plot 3
from scipy.stats import binom

n_corners = 8
lab_eff_reference = 0.4
counts = counting_measures_per_le[lab_eff_reference]["dataframe"].n_corners
hist =  np.bincount(counts, minlength=n_corners + 1)

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar(np.arange(9), hist, label="detected")
ax.set_xlabel("corners per NPC")
ax.set_ylabel("NPCs")
ax.legend()
if split_z:
    sides = "One side of NPC"
else:
    sides = "Both sides of NPC"
ax.set_title("Corners detected on NPC simuations. LE = " + str(lab_eff_reference) + "\n" + sides + "\n" + "Min locs per sector = " + str(min_locs_per_segment))

### plot 4

n_le = len(labelling_efficiencies)
nrows = 5
fig, ax = plt.subplots(nrows, n_le, figsize=(30,10))
for row in range(nrows):
    for i, lef in enumerate(labelling_efficiencies):
        ax[row,i].imshow(counting_measures_per_le[lef]["images"][row][0], cmap = "grey")
        ax[row,i].set_title(f"LE = {lef}")
        ax[row,i].set_xticks([])
        ax[row,i].set_yticks([])
filename = my_experiment_le.date_as_string + 'vlab4mic_NPCs_across_LE.png'
filename2 = os.path.join(my_experiment_le.output_directory, filename)
if probe_DoL is None:
    maintitle = "probe_DoL: max" + ". DoF: " + str(depth_of_field_nm)
else:
    maintitle = "probe_DoL: " + str(probe_DoL) + ". DoF: " + str(depth_of_field_nm)
fig.suptitle(maintitle)
plt.tight_layout()

## plot 5
# e.g. pores where only 3 or fewer corners were detected
#threshold_for_plot = 4
threshold_for_plot = [4, 6, 8]
lab_eff_reference = 0.3
_df = counting_measures_per_le[lab_eff_reference]["dataframe"]
for ncorners in threshold_for_plot:
    corners = results[_df.n_corners == ncorners]
    ids = list(_df.loc[_df["n_corners"] == ncorners].replicate)
    for example in ids[:1]:
        plot_npc_corners(counting_measures_per_le[lab_eff_reference]["rows"][example])

# plot 6

sns.boxplot(data=merged, x="source", y="radius_fit")
plt.xlabel("Labelling efficiency")
plt.ylabel("Radius fit")