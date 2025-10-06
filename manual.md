# VLab4Mic: A virtual laboratory for Microscopy

# 1. Installation

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


# 2. How to Use VLab4Mic

VLab4Mic is a Python library that can be used directly through our workflows and analysis modules or with interactive jupyter notebooks.
For detailed usage examples, see our [example scripts](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/examples).

## 2.1 Interactive Jupyter Notebooks

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

## 2.2 Using VLab4Mic as python library

A full documentation of the code can be found [here] <br>
vLab4Mic main functionality can be accessed with our high-level methods. For instance, the method "image_vsample" creates a virtual sample model and generates multimodal imaging acquisitions as shown in this example:

```python
from supramolsim.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()

```

Get a description of the available parameters of each method by using the built-in python method "help()":

```python
help(image_vsample)
```

# 3. VLab4Mic Methods


## 3.1 Image virtual sample


## 3.2 Run a parameter sweep analysis

