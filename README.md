# SupraMolSim

A Universal Validation Tool for Macromolecular Complexes in Fluorescence Microscopy

SupraMolSim is a library that focused on modelling Macromolecule Labelling and simulates their detection across imaging modalities to validate feature recovery under variable conditions.

Current features include:
- Creattion of a structural model from PDB data
- Direct and indirect labelling of PDB structures
- Modelling structural defects
- Structure detection as Image Generation for common microscopy modalities
- Parameter sweeps 
- And more!

# Codeless Jupyter Notebooks

[General workflow notebooks](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/basic_demo.ipynb)

Analysis notebooks: 

# Installation

SupraMolSim is compatible with Python 3.9, 3.10, 3.11 and 3.12 in MacOS, Windows and Linux.
It is advisable to create a new Conda environment to use SupraMolSim. To create a new python environment using Conda, refer to the official [website] (https://docs.anaconda.com/miniconda/)

1.- Create new environment
Run the following command (replace MYENV with your desired name for the environment):

```bash
conda create --name MYENV python=3.11
```

Then activate it by running:
```bash
conda activate MYENV
```


2.- Install SupraMolSim:

(Option 1)
Install SupraMolSim by simply running:
```shell
pip install SupraMolecularSimulator[jupyter]
```
You can remove "[jupyter]" from the command if your environment already has jupyter installed


(Option 2) To install the lastest release, run the following command:
```shell
pip install git+https://github.com/HenriquesLab/SupraMolecularSimulator.git
```

(Option 3)
To install from source distribution
-   Download the distribution from [here](https://github.com/HenriquesLab/SupraMolecularSimulator/dist/)
-   In your terminal, move to the directory containing the wheel and run 
```shell
pip install DIST_FILE_NAME.whl
```


# Wiki Tutorial

SupraMolSim, at its core, is designed as a collection of independent modules, each focused on a specific task within a typicall imaging experiment,
from choosing a macromolecule of interest, deciding between strategies to label it, to simulating the detection of a virtual sample across imaging modalities.

Notheless, integrative modules connect these functionalities into a complex workflow that allow us to branch the simulation and analysis.

We provide example configuration data for:
-   Models:
    -   Nuclear Pore Complex
    -   Immunoglobulin M
-   Labelling strategies
    -   Antibody labelling of Nup96
    -   NHS-ester direct labelling
-   Imaging modalities
    -   Widefield
    -   Confocal
    -   AiryScan
    -   STED
    -   Single Molecule Localisation Microscopy
-   Parameters
    -   Labelling efficiency
    -   Particle defects
    -   Particle orientation
-   Image metrics
    - SSIM

You can use SupraMolSim library through a jupyter notebook, with a graphical user interphase, or with a python stript according to your coding exprience.


### SupraMolSim Modules
- Atomic Structure parser
- Label design
- Labelled particle generator
- Virtual sample generator (Particle field)
- Virtual Microscope
- Experiment designer
- Analysis

### Integrated Workflows
- Load and parse an atomic model
- Create a labelled particle (multilabelling supported)
- Introduce particle defects
- Design a virtual sample
- Build your virtual microscope
- Generate multi-modal images of your virtual sample
- Use parameter sweeps acrross the workflow


# How to use SupraMolSim

SupraMolSim is a python library that can be used directly through our workflows and analysis module.
For a detail use see see our [script](https://github.com/HenriquesLab/SupraMolecularSimulator/)

## Codeless Jupyter Notebooks
However, part of experimental design requires a close examination of the models and parameters.
We provide codeless jupyter notebooks that requires no coding experience and allows to preview each step of the 
workflow.

### 1.- Following the installation instructions, activate your python environment where you installed SupraMolSim.
### 2.- Lauch Jupyter Lab by running this command:
```shell
jupyter lab
```
This command will open Jupyter Lab interfase in your web browser. In here, locate the folder containing SupraMolSim notebooks. You can download them from [here](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/)

### 3.- Select the notebook and start using SupraMolSim!
If you have any doubt during the use of the notebook refer to the [code documentation](https://github.com/HenriquesLab/SupraMolecularSimulator/) or rise an [issue](https://github.com/HenriquesLab/SupraMolecularSimulator/issues/).

# Contributing

# Licence

# Issues
Shoud you encounter any problem, do not hesitate on [letting us know](https://github.com/HenriquesLab/SupraMolecularSimulator/issues/). Don't forget to incude a detail description.

