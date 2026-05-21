# Methods

VLab4Mic exposes two high-level functions that cover the main simulation workflows.

---

## `image_vsample` — Create and Image a Virtual Sample

`image_vsample` creates a virtual sample from a structural model, applies fluorescent labelling, and generates imaging simulations across one or more microscopy modalities.

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample()
```

The function returns:

- **`images`** — dictionary of simulated images with noise, keyed by modality name
- **`noiseless`** — dictionary of noiseless simulated images, keyed by modality name  
- **`experiment`** — object containing all modules used to generate the simulation

### Basic Usage

Running with default parameters creates a virtual sample with a default structure and probe:

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(run_simulation=True)
```

### Parameterized Usage

Specify each aspect of the virtual sample and imaging:

```python
from vlab4mic.experiments import image_vsample

modalities = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

images, noiseless, experiment = image_vsample(
    structure="7R5K",              # PDB ID for a Nuclear Pore Complex
    probe_template="Antibody",     # Probe template
    probe_target_type="Sequence",
    probe_target_value="ELAVGSL",  # Epitope sequence
    number_of_particles=10,
    random_rotations=True,
    rotation_angles=None,
    multimodal=modalities,
    STED={"exp_time": 0.01},       # Modality-specific parameters
    expansion_factor=1,            # Scale factor for structure coordinates
    random_seed=42,                # Set for reproducibility
    run_simulation=True,
    clear_experiment=True,
)
```

### Using Structure-Specific Probe Templates

For common structures, pre-configured probe templates are available:

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(
    structure="7R5K",
    probe_templates=["NPC_Nup96_Cterminal_direct"],
    run_simulation=True,
    clear_experiment=True,
)
```

See the [Templates](templates.md) page for the full list of available probes.

### Indirect Labelling (Primary + Secondary Probe)

For indirect labelling, pass `primary_probe` and `secondary_probe` instead of a single `probe_template`:

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(
    structure="3J3Y",
    primary_probe="HIV_capsid_p24_direct",
    secondary_probe="anti-p24_primary_antibody_HIV",
    multimodal=["SMLM"],
    run_simulation=True,
)
```

### Displaying Simulation Results

```python
import matplotlib.pyplot as plt

nmods = len(modalities)
fig, axs = plt.subplots(1, nmods)
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][0], cmap="magma")
    axs[i].set_title(mod)
plt.show()
```

---

## `run_parameter_sweep` — Sweep Over Parameter Combinations

`run_parameter_sweep` systematically generates and analyses simulations over combinations of parameter values.

```python
from vlab4mic.sweep_generator import run_parameter_sweep

sweep_gen = run_parameter_sweep(
    output_name="my_sweep",
    return_generator=True,
    save_sweep_images=True,
    save_analysis_results=True,
    run_analysis=True,
)
```

### Specifying Sweep Parameters

Pass parameter values as a **list** of explicit values, or as a **tuple `(min, max, nsteps)`** to generate linearly spaced values:

```python
from vlab4mic.sweep_generator import run_parameter_sweep

sweep_gen = run_parameter_sweep(
    structures=["7R5K"],
    probe_templates=["NPC_Nup96_Cterminal_direct"],
    sweep_repetitions=20,
    # Parameters to sweep
    labelling_efficiency=(0, 1, 5),              # 5 values from 0 to 1
    structural_integrity=(0, 1, 5),              # 5 values from 0 to 1
    structural_integrity_small_cluster=[300],    # Single value
    structural_integrity_large_cluster=[600],    # Single value
    exp_time=[0.001, 0.01],                      # Two values
    # Output options
    output_name="my_sweep",
    return_generator=True,
    save_sweep_images=True,
    save_analysis_results=True,
    run_analysis=True,
    random_seed=42,                              # Reproducible sweep
)
```

!!! warning "Combinatorial growth"
    All parameter combinations are tested. Setting 10 values for `labelling_efficiency` and 10 values for `structural_integrity` generates **100 combinations**. Plan accordingly.

### Available Sweep Parameters

The following parameters can be swept (pass `None` to use the default value):

| Category | Parameter |
|----------|-----------|
| **Probe** | `probe_target_type`, `probe_target_value`, `probe_target_option`, `probe_model`, `probe_fluorophore`, `probe_paratope`, `probe_conjugation_target_info`, `probe_seconday_epitope`, `peptide_motif` |
| **Probe geometry** | `probe_distance_to_epitope`, `probe_steric_hindrance`, `probe_DoL`, `probe_wobble_theta` |
| **Labelling** | `labelling_efficiency` |
| **Structural integrity** | `structural_integrity`, `structural_integrity_small_cluster`, `structural_integrity_large_cluster` |
| **Virtual sample** | `sample_dimensions`, `particle_orientations`, `rotation_angles`, `minimal_distance` |
| **Imaging** | `pixelsize_nm`, `lateral_resolution_nm`, `axial_resolution_nm`, `psf_voxel_nm`, `depth_of_field_nm` |
| **Acquisition** | `exp_time` |
