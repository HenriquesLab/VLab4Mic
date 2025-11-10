<img src="src/logo/logo.png" align="right" width="200" style="margin-left: 20px;"/>

# VLab4Mic: A Virtual Laboratory for Microscopy

VLab4Mic is a library for validation of experimental design, from choosing a macromolecule to study, to its imaging through several imaging modalities.

Current features include:
- **Creation** of structural models from PDB/CIF data
- **Direct and indirect labeling** of PDB/CIF structures with various probes
- **Modeling of structural defects** and variations
- **Image generation** for common microscopy modalities (widefield, confocal, STED, etc.)
- **Parameter sweeps** for experimental design optimization
- **Quantitative analysis** tools for validation

# Installation

VLab4Mic is compatible with **Python 3.9, 3.10, 3.11, and 3.12** on macOS, Windows, and Linux.

**We strongly recommend** creating a new Conda environment to use VLab4Mic. To create a new Python environment using Conda, refer to the official [Miniconda documentation](https://docs.anaconda.com/miniconda/).

## Step 1: Create a New Environment

Run the following command (replace `MYENV` with your desired environment name):

```bash
conda create --name MYENV python=3.11
```

Then activate it:
```bash
conda activate MYENV
```

## Step 2: Install VLab4Mic

Currently, VLab4Mic is available through Test PyPI (PyPI release coming soon!):
Run the following command to install vlab4mic with the necesary libraries to support our Jupyter notebooks:

```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "vlab4mic[jupyter]"
```

# How to Use VLab4Mic

VLab4Mic can be used in three ways, depending on your coding experience:
1. **Jupyter Notebooks with GUI** - No coding required (recommended for beginners)
2. **Jupyter Notebooks with code** - Moderate coding (recommended for intermediate users)
3. **Python scripts** - Full coding control (recommended for advanced users)

Refer to the [manual](https://github.com/HenriquesLab/VLab4Mic/blob/main/manual.md) for detailed instructions to use.


## Python Script Example
VLab4Mic is a Python library that can be used directly through our workflows and analysis modules.
For detailed usage examples, see our [example scripts](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/examples).
Here's a simple example to get you started with VLab4Mic:

```python
from vlab4mic.experiments import image_vsample
import matplotlib.pyplot as plt

modalities = ["Widefield", "Confocal", "STED", "SMLM"]

## use a high-level function to parameterise an virtual sample and its image simulation
images, noiseless, experiment = image_vsample(
    structure="7R5K",  # PDB ID code for a Nuclear Pore complex
    probe_template="Antibody",  # Probe template for an antibody
    probe_target_type="Sequence",  # Epitope type
    probe_target_value="ELAVGSL",  # Epitope sequence
    number_of_particles=10, # Number of indpependent copies of a labelled structure in the sample
    random_rotations=True, # Rotation in XY plane
    rotation_angles=None,
    multimodal=modalities, # Imaging modalities
    STED={"exp_time": 0.01},  # modality-specific parameters
    run_simulation=True, # Simulate image acquisition after creating virtual sample
    clear_experiment=True, # Clear default values of the experiment
)
# visualize results
nmods = len(modalities)
fig, axs = plt.subplots(1, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][nframe], cmap="magma")
    axs[i].set_title(mod)
plt.show()
```

## Jupyter Notebook Usage

It is also possible to use Jupyter Notebooks for running custom code such as the example above in scripts. For GUI-based usage without coding in a Jupyter Lab, run the following code boxes in separate code cells:

```python
# Initialisaiton cell
from vlab4mic.jupyter_widgets import experiment_widgets
from vlab4mic import experiments

## Initialize experiment
my_experiment = experiments.ExperimentParametrisation()
```
VLab4Mic provides custom jupyter widgets to interact with the experiment object.
```python
# Display interactive widgets for specific modules, for instance to select structure
experiment_widgets.select_structure_widget(my_experiment).show()
```
```python
# Run experiment
experiment_widgets.run_experiment_widget(my_experiment).show()
```

## Interactive Jupyter Notebooks

| Category | Description | Notebook | Colab Link |
| --- | --- | --- | --- |
| **Main Interface** | Create virtual samples and simulate image acquisition with multiple imaging modalities | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) |
| **Parameter Sweeps** | Generate and analyze simulations over parameter ranges for optimization | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) |


# Documentation & Tutorials

VLab4Mic is designed as a collection of independent modules that work together through our workflows. Comprehensive documentation and tutorials are available on our [Wiki Tutorial page](https://github.com/HenriquesLab/SupraMolecularSimulator/wiki).


## Pre-configured Examples

We provide ready-to-use models and configurations for:

### Structures:
- **Nuclear Pore Complex** - Multi-protein assemblies with 8-fold symmetry
- **Clathrin Coated Pit** - Membrane trafficking structures
- **Matured HIV-capsid core** - Viral capsid assemblies
- **Bacteriophage T4 capsid shell** - Capsid shell of Bacteriophage T4

### Labelling Probes:
- **Primary/Secondary Antibodies** - Traditional immunofluorescence
- **Nanobodies** - Single-domain antibodies for improved access
- **Fluorescent Proteins** (GFP, mCherry, etc.) - Genetically encoded tags
- **Chemical Linkers** - NHS-ester, click chemistry probes
- **And many more** specialized labeling strategies

### Imaging Modalities:
- **Widefield Microscopy** - Basic fluorescence imaging
- **Confocal Microscopy** - Optical sectioning capability
- **AiryScan** - Enhanced resolution confocal
- **STED** - Super-resolution nanoscopy
- **Single Molecule Localization Microscopy** (SMLM) - Nanometer precision

Refer to the [manual](https://github.com/HenriquesLab/VLab4Mic/blob/main/manual.md) for detailed tables of templates and parameters.

# Contributing

Contributions are very welcome! Please read our [Contribution Guidelines](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/CONTRIBUTING.md) to learn how to get involved.

# License

Distributed under the terms of the [MIT license](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/LICENSE.txt), VLab4Mic is free and open source software.

# Issues & Support

Should you encounter any problems, please don't hesitate to [let us know](https://github.com/HenriquesLab/SupraMolecularSimulator/issues). When reporting issues, please include:
- A detailed description of the problem
- Steps to reproduce the issue
- Your operating system and Python version
- Any error messages or screenshots

For general questions, please use our [Discussions](https://github.com/HenriquesLab/SupraMolecularSimulator/discussions) section.

