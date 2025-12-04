import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axes_grid1
from ipywidgets import interact, widgets
import numpy as np
import copy
from scipy.ndimage import rotate


def add_ax_scatter(plotobj, trgt_dictionary, fraction=1):
    if trgt_dictionary["coordinates"].shape[0] == 0:
        print(f"No emitters to show.")
    else:
        if fraction == 1:
            plotobj.scatter(
                trgt_dictionary["coordinates"][:, 0],
                trgt_dictionary["coordinates"][:, 1],
                trgt_dictionary["coordinates"][:, 2],
                c=trgt_dictionary["plotcolour"],
                label=trgt_dictionary["label_name"],
                s=trgt_dictionary["plotsize"],
                alpha=trgt_dictionary["plotalpha"],
                marker=trgt_dictionary["plotmarker"],
                depthshade=True,
            )
        else:
            n = ((trgt_dictionary["coordinates"]).shape)[0]
            print(f"Showing {n*fraction} atoms for {trgt_dictionary['label_name']}")
            ids = np.random.choice(np.arange(0, n), int(n * fraction), replace=False)
            subset = trgt_dictionary["coordinates"][ids, :]
            plotobj.scatter(
                subset[:, 0],
                subset[:, 1],
                subset[:, 2],
                c=trgt_dictionary["plotcolour"],
                label=trgt_dictionary["label_name"],
                s=trgt_dictionary["plotsize"],
                alpha=trgt_dictionary["plotalpha"],
                marker=trgt_dictionary["plotmarker"],
            )


def draw1nomral_segment(points_normal, figure, lenght=100, colors=["g", "y"]):
    # points_normals is a list of 2 elements
    # first element are the normals
    # second element are ponts in space
    starts = points_normal["pivot"]
    normal = points_normal["direction"]  # this ones might not be normalized to 1
    normalized = normal / np.linalg.norm(normal)
    ends = starts + normalized * lenght
    figure.plot(
        [starts[0], ends[0]],
        [starts[1], ends[1]],
        [starts[2], ends[2]],
        color=colors[0],
    )
    figure.scatter(ends[0], ends[1], ends[2], color=colors[1], marker="o")


def draw_nomral_segments(points_normal, figure, lenght=100, colors=["g", "y"], **kwargs):
    # points_normals is a list of 2 elements
    # first element are the normals
    # second element are ponts in space
    starts = points_normal[1]
    normals = points_normal[0]  # this ones might not be normalized to 1

    for i in range(starts.shape[0]):
        normalized = normals[i]
        normalized /= np.linalg.norm(normalized, axis=0)
        ends = starts[i] + normalized * lenght
        figure.plot(
            [starts[i][0], ends[0]],
            [starts[i][1], ends[1]],
            [starts[i][2], ends[2]],
            color=colors[0],
        )
        # figure.scatter(ends[0], ends[1], ends[2], color = 'r', marker = "")
    figure.scatter(
        starts[:, 0], starts[:, 1], starts[:, 2], color=colors[1], marker="o"
    )
    return figure


def stack_projection(stack, angle=45, axes=(1,2), method = "sd", zdepth = False):
    projection = None
    zstack = copy.copy(stack)
    rotated = rotate(zstack, angle=angle, axes=axes, reshape=False, order=3)  # order=3: cubic interpolation
    if zdepth:
        nslices = rotated.shape[2]
        for i in range(1, nslices):
            #intensity_factor = np.log(1.8 + ((nslices - i)/(nslices)))
            intensity_factor = 0.6 * np.exp(((nslices - i)/(nslices)) * (1/(i+1)))
            rotated[:,:,i] = rotated[:,:,i]*intensity_factor
    if method == "max":
        projection = np.max(rotated, axis=2)
    elif method == "sd":
        projection = np.std(rotated, axis=2)
    elif method == "mean":
        projection = np.mean(rotated, axis=2)
    elif method == "sum":
        projection = np.sum(rotated, axis=2)
    return projection

def plot_projection(stack, angle=45, plane="XY", method = "max", zdepth=False):
    if plane == "XY":
        axes = (1, 2)
    elif plane == "XZ": 
        axes = (0, 2)
    elif plane == "YZ":
        axes = (0, 1)
    else:
        raise ValueError("Invalid plane. Choose from 'XY', 'XZ', or 'YZ'.")
    projection = stack_projection(stack, angle=angle, axes=axes, method=method, zdepth=zdepth)
    plt.imshow(projection, cmap='gray')
    plt.axis('off')
    plt.show()