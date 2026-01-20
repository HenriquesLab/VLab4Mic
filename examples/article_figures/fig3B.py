from vlab4mic import sweep_generator
import matplotlib.pyplot as plt
import os
import numpy as np
np.random.seed(44)

sweep_gen = sweep_generator.run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],
    sweep_repetitions=20,
    # parameters for sweep
    labelling_efficiency=(0,1,0.25),
    structural_integrity=(0,1,0.25),
    structural_integrity_small_cluster=[300,],
    structural_integrity_large_cluster=[600,],
    exp_time=[0.001, 0.01,],
    # output and analysis
    output_name="vlab_script",
    return_generator=True,
    save_sweep_images=True,
    save_analysis_results=True,
    run_analysis=True
    )


plt.rcParams['figure.figsize'] = [20, 10]

fig, axs = plt.subplots(1, 4)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=4, structural_integrity_parameters=0, modality_template=3, return_image=True)
axs[0].imshow(image, cmap="grey")
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[0].set_title(title)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=1, structural_integrity_parameters=0, modality_template=3, return_image=True)
axs[1].imshow(image, cmap="grey")
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[1].set_title(title)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=4, structural_integrity_parameters=2, modality_template=3, return_image=True)
axs[2].imshow(image, cmap="grey")
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[2].set_title(title)
image, parameters = sweep_gen.preview_image_output_by_ID(probe_parameters=1, structural_integrity_parameters=2, modality_template=3, return_image=True)
axs[3].imshow(image, cmap="grey")
title = "Labelling efficiency: " + str(parameters[2]["labelling_efficiency"]) + "\n" + "Structural Integrity: " + str(parameters[3]["structural_integrity"])
axs[3].set_title(title)


filename = os.path.join(sweep_gen.output_directory, 'vlab4mic_fig3B.png')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()