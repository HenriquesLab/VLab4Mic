import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axes_grid1
from ipywidgets import interact, widgets
import numpy as np



def slider_normalised(stack, dimension):
    def frame_slider_norm(frame):
        '''
        stack is assumed to be of the shape MxNxZ, where Z is the axial direction
        '''
        stack_max = np.max(stack)
        fig = plt.figure()
        grid = axes_grid1.AxesGrid(
            fig, 111, nrows_ncols=(1, 1), axes_pad=0.5, cbar_location="right",
            cbar_mode="each", cbar_size="15%", cbar_pad="5%",)
        if dimension == 0:
            im0 = grid[0].imshow(stack[frame-1,:,:], cmap='gray', interpolation='none', vmin=0, vmax=stack_max)
        elif dimension == 1:
            im0 = grid[0].imshow(stack[:,frame-1,:], cmap='gray', interpolation='none', vmin=0, vmax=stack_max)
        elif dimension == 2:
            im0 = grid[0].imshow(stack[:,:,frame-1], cmap='gray', interpolation='none', vmin=0, vmax=stack_max)
        grid.cbar_axes[0].colorbar(im0)

    interact(frame_slider_norm, frame=widgets.IntSlider(min=1, max=stack.shape[dimension], step=1, value=0, continuous_update=False))


def add_ax_scatter(plotobj, trgt_dictionary, fraction=1):
    if fraction == 1:
        plotobj.scatter(trgt_dictionary["coordinates"][:, 0],
                        trgt_dictionary["coordinates"][:, 1],
                        trgt_dictionary["coordinates"][:, 2],
                        c=trgt_dictionary["plotcolour"], label=trgt_dictionary["label_name"],
                        s=trgt_dictionary["plotsize"], alpha=trgt_dictionary["plotalpha"],
                        marker=trgt_dictionary["plotmarker"],
                        depthshade=True)
    else:
        n = ((trgt_dictionary["coordinates"]).shape)[0]
        print(f"Showing {n*fraction} atoms for {trgt_dictionary['label_name']}")
        ids = np.random.choice(np.arange(0, n), int(n*fraction), replace=False)
        subset = trgt_dictionary["coordinates"][ids, :]
        plotobj.scatter(subset[:, 0],
                        subset[:, 1],
                        subset[:, 2],
                        c=trgt_dictionary["plotcolour"], label=trgt_dictionary["label_name"],
                        s=trgt_dictionary["plotsize"], alpha=trgt_dictionary["plotalpha"],
                        marker=trgt_dictionary["plotmarker"])


def draw1nomral_segment(points_normal, figure, lenght=100, colors=['g', 'y']):
    # points_normals is a list of 2 elements, first are the normals, and seconds are ponts in space
    starts = points_normal["pivot"]
    normal = points_normal["direction"]  # this ones might not be normalized to 1
    normalized = normal / np.linalg.norm(normal)
    ends = starts + normalized * lenght
    figure.plot([starts[0], ends[0]], [starts[1], ends[1]], [starts[2], ends[2]], color=colors[0])
    figure.scatter(ends[0], ends[1], ends[2], color = colors[1], marker = "o")


def draw_nomral_segments(points_normal, figure, lenght=100, colors=['g', 'y']):
    # points_normals is a list of 2 elements, first are the normals, and seconds are ponts in space
    starts = points_normal[1]
    normals = points_normal[0]  # this ones might not be normalized to 1

    for i in range(starts.shape[0]):
        normalized = normals[i]
        normalized /= np.linalg.norm(normalized, axis=0)
        ends = starts[i] + normalized * lenght
        figure.plot([starts[i][0], ends[0]], [starts[i][1], ends[1]], [starts[i][2], ends[2]], color = colors[0])
        # figure.scatter(ends[0], ends[1], ends[2], color = 'r', marker = "")
    figure.scatter(starts[:,0], starts[:,1], starts[:,2], color = colors[1], marker = "o")
    return figure