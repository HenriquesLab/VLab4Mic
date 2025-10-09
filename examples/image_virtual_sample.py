from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt

modalities = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

images, noiseless, experiment = image_vsample(
    structure="7R5K",  # PDB ID code for a Nuclear Pore complex
    probe_template="Antibody",  # Probe template for an antibody
    probe_target_type="Sequence",  
    probe_target_value="ELAVGSL",  # epitope sequence
    number_of_particles=10,
    random_rotations=True,
    rotation_angles=None,
    multimodal=modalities,
    STED={"exp_time": 0.01},  # modality-specific parameters
    run_simulation=True,
    clear_experiment=True,
)
# visualize results
nmods = len(modalities)
fig, axs = plt.subplots(1, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][nframe], cmap="magma")
    axs[i].set_title(mod)
plt.show()