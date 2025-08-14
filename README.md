# VLab4Mic: A virtual laboratory for Microscopy

VLab4Mic is a library for validation of experimental design from choosing a macromolecule to study, to its imaging through several imaging modalities.

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

Currently, VLab4Mic is available through GitHub (PyPI release coming soon!):

```bash
pip install git+https://github.com/HenriquesLab/SupraMolecularSimulator.git
```

**Coming soon:** Google Colab integration for browser-based usage!

# Documentation & Tutorials

VLab4Mic is designed as a collection of independent modules that work together through our workflows. Comprehensive documentation and tutorials are available on our [Wiki Tutorial page](https://github.com/HenriquesLab/SupraMolecularSimulator/wiki).

## Pre-configured Examples

We provide ready-to-use models and configurations for:

### Structures:
- **Nuclear Pore Complex** - Multi-protein assemblies with 8-fold symmetry
- **Clathrin Coated Pit** - Membrane trafficking structures
- **Matured HIV-capsid core** - Viral capsid assemblies

### Labeling Probes:
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

## Usage Options

VLab4Mic can be used in three ways, depending on your coding experience:
1. **Jupyter Notebooks with GUI** - No coding required (recommended for beginners)
2. **Jupyter Notebooks with code** - Moderate coding (recommended for intermediate users)
3. **Python scripts** - Full coding control (recommended for advanced users)


# Quick Start

## Python Script Example

Here's a simple example to get you started with VLab4Mic:

```python
import supramolsim as vlm

# Load a structure (e.g., Nuclear Pore Complex)
structure, params = vlm.workflows.load_structure(
    structure_id="5A9Q",  # NPC structure
    config_dir="path/to/your/config"
)

# Create a labeled instance
labeled_structure = vlm.workflows.create_labeled_structure(
    structure=structure,
    probe_config="antibody_primary",
    labeling_strategy="epitope_based"
)

# Generate a virtual sample
virtual_sample = vlm.workflows.create_virtual_sample(
    labeled_structure=labeled_structure,
    n_particles=50,
    random_orientations=True,
    sample_size=[10, 10, 2]  # micrometers
)

# Simulate imaging
imager = vlm.workflows.simulate_imaging(
    virtual_sample=virtual_sample,
    modality="confocal",
    pixel_size=20,  # nanometers
    psf_config="standard_confocal"
)

# Generate and save images
images = imager.generate_images()
imager.save_results("output_directory/")
```

## Jupyter Notebook Usage

For GUI-based usage without coding:

```python
from supramolsim.jupyter_widgets import experiment_widgets
from supramolsim import experiments

# Initialize experiment
my_experiment = experiments.ExperimentParametrisation()

# Display interactive widgets
experiment_widgets.select_structure_widget(my_experiment).show()
experiment_widgets.select_probe_widget(my_experiment).show()
experiment_widgets.select_sample_parameters_widget(my_experiment).show()
```

# How to Use VLab4Mic

VLab4Mic is a Python library that can be used directly through our workflows and analysis modules.
For detailed usage examples, see our [example scripts](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/examples).

## Interactive Jupyter Notebooks

| Category | Description | Notebook | Colab Link |
| --- | --- | --- | --- |
| **Main Interface** | Create virtual samples and simulate image acquisition with multiple imaging modalities | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) |
| **Parameter Sweeps** | Generate and analyze simulations over parameter ranges for optimization | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) |

## Getting Started with Notebooks

### Step 1: Activate Your Environment
Following the installation instructions above, activate your Python environment:
```bash
conda activate MYENV  # Replace MYENV with your environment name
```

### Step 2: Launch Jupyter Lab
```bash
jupyter lab
```
This will open Jupyter Lab in your web browser.

### Step 3: Download and Open Notebooks
Download the notebooks from our [repository](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/notebooks) and open them in Jupyter Lab.

### Step 4: Start Experimenting!
Follow the interactive widgets and instructions in each notebook. For questions, refer to our [documentation](https://github.com/HenriquesLab/SupraMolecularSimulator/wiki) or [create an issue](https://github.com/HenriquesLab/SupraMolecularSimulator/issues).

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

