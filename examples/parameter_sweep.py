from vlab4mic.sweep_generator import run_parameter_sweep
import matplotlib.pyplot as plt

sweep_gen = run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],  # Probe template tailored to 7R5K
    sweep_repetitions=3,
    # parameters for sweep
    labelling_efficiency=(0, 1, 0.5),  # values between 0 and 1 with step of 0.5
    structural_integrity=(0, 1, 0.5),  # values between 0 and 1 with step of 0.5
    structural_integrity_small_cluster=[300,],  # 1 single value 
    structural_integrity_large_cluster=[600,],  # 1 single value 
    exp_time=[0.001, 0.01,],  # 2 values
    # output and analysis
    output_name="vlab_example_sweep",
    return_generator=True,
    save_sweep_images=True,  # By default, the saving directory is set to the home path of the user
    save_analysis_results=True,
    run_analysis=True
)
## Once our sweep has been generated, we can preview images corresponding to specific parameter combinations.
## In order to retrieve images based on their parameters used, we use the method `preview_image_output_by_ID`,
## which takes as input the IDs corresponding to each parameter (except for the structure) in the following order:

## - probe_model
## - probe_parameters
## - structural_integrity_parameters
## - virtual_sample_parameters
## - Modality_template
## - Modality_parameters
## - Acquisition_parameters

## You can use sweep_gen.acquisition_outputs_parameters get a dictionary where the key is
## composed by each parameter ID separated by underscores.
## For instance, the first combination is always "0_0_0_0_0_0_0". When a parameter is not used in the sweep,
## its ID is always 0.

## For this example, since we used 3 values for labelling efficiency and structural_integrity,
## the possible IDs for labelling efficiency are 0, 1 or 2 (for 0, 0.5 and 1, respectively). 
## Similarly, the IDs for structural_integrity are also 0, 1, 2 (for 0, 0.5 and 1, respectively).

## We will show the simulations for SMLM modality (modality_template ID = 3) in defaults.
fig, axs = plt.subplots(1, 2)
# First image: labelling efficiency of 1 and no structural_integrity
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=2, structural_integrity_parameters=0, modality_template=3, return_image=True)
axs[0].imshow(image, cmap="grey")
axs[0].set_xticks([])
axs[0].set_yticks([])
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[0].set_title(title)
# Second image: labelling efficiency of 1 and structural_integrity of 0.5
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=2, structural_integrity_parameters=1, modality_template=3, return_image=True)
axs[1].imshow(image, cmap="grey")
axs[1].set_xticks([])
axs[1].set_yticks([])
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[1].set_title(title)

plt.show()