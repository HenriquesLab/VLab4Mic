import numpy as np
import os
import matplotlib.pyplot as plt
from vlab4mic import experiments
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
np.random.seed(44)

structure = "PATH/TO/YOUR/LOCAL/MODEL/AF-P00520-F1-model_v6.cif"

# Check if the structure file exists
if not os.path.isfile(structure):
    raise FileNotFoundError(
        f"Structure file not found at: {structure}\n\n"
        "## Instructions to obtain a custom structure for this example:\n"
        "# 1. Go to AlphaFold Protein Structure Database:\n"
        "#    https://alphafold.ebi.ac.uk/\n"
        "# 2. Search for your protein of interest using its UniProt ID\n"
        "#    for instance: https://alphafold.ebi.ac.uk/entry/AF-P00520-F1\n"
        "#    (For this example, we used AF-P00520-F1-v6)\n"
        "# 3. Download the CIF (mmCIF) file of the predicted structure,\n"
        "# 4. Set 'structure' to the path of the downloaded CIF file.\n"
        "# Example:\n"
        "#    structure = 'PATH/TO/YOUR/LOCAL/MODEL/AF-P00520-F1-model_v6.cif'\n"
        "# 5. Run the python script with the new structure path."
    )

modalities = ["STED", "SMLM",]
# Run simulation
images , noiseless_,  experiment = experiments.image_vsample(
    structure=structure,
    structure_is_path=True,
    clear_experiment=True,
    multimodal=modalities,
    run_simulation=True
)

# expand sample for better visualization
experiment.coordinate_field.expand_isotropically(factor=10.0)
experiment.build(modules=["imager",])
images , noiseless =  experiment.run_simulation()

# Plotting figure with 3 panels
# panel 1
fig = plt.figure(figsize=[20,20])
ax = fig.add_subplot(131, projection="3d")
experiment.structure.show_target_labels(
    axis_object=ax,
    with_assembly_atoms=True,
    assembly_fraction=0.1,
    atoms_alpha=0.1,
    show_axis=False,
    axesoff=False,
    reference_point=False
)
ax.set_title('Target sites on structure \n' + f' ({structure_name})')
ax.set_ylabel('Angstroms')
ax.set_xticklabels([])
# panel 2
ax = fig.add_subplot(132, projection="3d")
experiment.coordinate_field.show_field(axis_object=ax, view_init=[90,0,0], emitters_plotsize=4)
ax.set_title('Virtual Sample \n (expanded 10x)')
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_zlabel(None)
ax.set_ylabel(None)
ax.set_xlabel('Nanometers')
#panel 3

ax = fig.add_subplot(133)
ax.imshow(images["SMLM"]["ch0"][0], cmap="grey")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title('Simulated SMLM Image')
scalebar = AnchoredSizeBar(ax.transData,
                           100, '500 nm', 'lower right', 
                           pad=2,
                           color='white',
                           frameon=False,
                           size_vertical=5)
ax.add_artist(scalebar)
filename = os.path.join(experiment.output_directory, 'vlab4mic_alphafold_from_local_file.pdf')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()