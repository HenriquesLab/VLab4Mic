import matplotlib.pyplot as plt
from vlab4mic import experiments
from mpl_toolkits import axes_grid1
import os
random_seed=24

modalities = ["Widefield", "STED", "SMLM"]
dof = [10, 50, 200, 300, 500, 1000]
images , noiseless_,  experiment = experiments.image_vsample(
    structure="7R5K",
    probe_template="NPC_Nup96_Cterminal_direct",
    clear_experiment=True,
    sample_dimensions=[1000,1000,150],
    particle_positions=[[0.2,0.2,0.1],[0.8,0.8,0.9]],
    multimodal=modalities,
    run_simulation=False,
    random_seed=random_seed
)

fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(121, projection="3d")
experiment.coordinate_field.show_field(view_init=[0,0,0], axis_object=ax, initial_pos=False, emitters_plotsize=10)
ax.set_xticklabels([])
ax.set_zticks([0,50, 100,150])
ax.set_yticks([0,1000])
ax.set_xlabel(None)
ax = fig.add_subplot(122, projection="3d")
experiment.coordinate_field.show_field(view_init=[90,0,0], axis_object=ax, initial_pos=False, emitters_plotsize=10)
ax.set_zticklabels([])
ax.set_zlabel(None)
fig.tight_layout()

filename = os.path.join(experiment.output_directory, 'figS_depth_of_field_panelA.pdf')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()

def generate_depth_examples(experiment = None, modality_name="STED", depth_of_field_nm_list=None, ):
    images_list = []
    depths_list = []
    minmax_list = []
    max_value = 0
    min_value = 1000
    #experiment.set_modality_acq(modality_name=modality_name, exp_time=0.001, noise=True)
    for depth_of_field_nm in depth_of_field_nm_list:
        experiment.update_modality(
            modality_name=modality_name,
            depth_of_field_nm=depth_of_field_nm,)
        images, noiseless = experiment.run_simulation(modality=modality_name)
        img = images[modality_name]["ch0"][0]
        depth_slices = experiment.imager.modalities[modality_name]["psf"]["depth"]
        pixelsize = experiment.imager.modalities[modality_name]["psf"]["voxelsize"][2]
        depths_list.append(depth_slices*pixelsize)
        minmax_list.append((img.min(), img.max()))
        if img.max() > max_value:
            max_value = img.max()
        if img.min() < min_value:
            min_value = img.min()
        images_list.append(img)
    return {"images_list": images_list, "depths_list": depths_list, "minmax_list":minmax_list, "modality_name":modality_name}


def plot_examples(images_list, modality_name, depths_list, minmax_list, with_cbar=True):   
    fig = plt.figure(figsize=(20, 5))
    cbar_mode = "each"
    if with_cbar:
        grid = axes_grid1.AxesGrid(
                fig,
                111,
                nrows_ncols=(1, int(len(depths_list))),
                axes_pad=1,
                cbar_location="bottom",
                cbar_mode=cbar_mode,
                cbar_size="15%",
                cbar_pad="5%",
            )
    else:
        grid = axes_grid1.AxesGrid(
                fig,
                111,
                nrows_ncols=(1, int(len(depths_list))),
                axes_pad=1,
                cbar_mode=None,
            )

    normalise = True
    for i, images in enumerate(images_list):
        if normalise:
            im0 = grid[i].imshow(images, cmap="gray")
        else:
            im0 = grid[i].imshow(images, cmap="gray", vmin=minmax_list["min"], vmax=minmax_list["max"])
        grid[i].set_title(f"Field of Depth: {depths_list[i]} nm")
        grid[i].axis('off')
        if with_cbar:
            cbar = grid.cbar_axes[i].colorbar(im0)

    fig.suptitle("Chaning Field of Depth with fixed focused plane for " + modality_name + " modality")
    plt.tight_layout()
    
    filename = os.path.join(experiment.output_directory, modality_name + '_figS_depth_of_field_panelB.pdf')
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

sted_examples = generate_depth_examples(experiment=experiment, modality_name="STED", depth_of_field_nm_list=dof)
plot_examples(**sted_examples, with_cbar=False)
smlm_examples = generate_depth_examples(experiment=experiment, modality_name="SMLM", depth_of_field_nm_list=dof)
plot_examples(**smlm_examples, with_cbar=False)
widefield_examples = generate_depth_examples(experiment=experiment, modality_name="Widefield", depth_of_field_nm_list=dof)
plot_examples(**widefield_examples, with_cbar=False)