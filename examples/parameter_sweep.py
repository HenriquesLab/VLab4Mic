from vlab4mic.sweep_generator import run_parameter_sweep
import matplotlib.pyplot as plt

sweep_gen = run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],  # Probe template tailored to 7R5K
    sweep_repetitions=3,
    # parameters for sweep
    labelling_efficiency=(0, 1, 3),  # 3 linearly spaced values between 0 and 1
    defect=(0, 1, 3),  # 3 linearly spaced values between 0 and 1
    defect_small_cluster=[300,],  # 1 single value 
    defect_large_cluster=[600,],  # 1 single value 
    exp_time=[0.001, 0.01,],  # 2 values
    # output and analysis
    output_name="vlab_example_sweep",
    return_generator=True,
    save_sweep_images=True,  # By default, the saving directory is set to the home path of the user
    save_analysis_results=True,
    run_analysis=True
)

fig, axs = plt.subplots(1, 2)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=2, defect_parameters=0, modality_template=3, return_image=True)
axs[0].imshow(image, cmap="grey")
axs[0].set_xticks([])
axs[0].set_yticks([])
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Incomplete Labelling: " + str(parameters[3]["defect"])
axs[0].set_title(title)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=2, defect_parameters=1, modality_template=3, return_image=True)
axs[1].imshow(image, cmap="grey")
axs[1].set_xticks([])
axs[1].set_yticks([])
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Incomplete Labelling: " + str(parameters[3]["defect"])
axs[1].set_title(title)

plt.show()