# VLab4Mic: A virtual laboratory for Microscopy

VLab4Mic is a library for validation of experimental design from choosing a macromolecule to study, to its imaging through several imaging modalities.

Current features include:
- Creattion of a structural model from PDB/CIF data
- Direct and indirect labelling of PDB/CIF structures
- Modelling structural defects
- Structure detection as Image Generation for common microscopy modalities
- Parameter sweeps 
- And more!

# Installation

VLab4Mic is compatible with Python 3.9, 3.10, 3.11 and 3.12 in MacOS, Windows and Linux.
It is advisable to create a new Conda environment to use VLab4Mic. To create a new python environment using Conda, refer to the official [website](https://docs.anaconda.com/miniconda/)

1.- Create new environment
Run the following command (replace MYENV with your desired name for the environment):

```bash
conda create --name MYENV python=3.11
```

Then activate it by running:
```bash
conda activate MYENV
```


2.- Install VLab4Mic:

Currently, you can only access VLab4Mic if
you have access to the repository. If this is your case, run the following command:
```shell
pip install git+https://github.com/HenriquesLab/SupraMolecularSimulator.git
```
Soon you'll be able to acces it in a goolge colab! Stay tuned!

# Wiki Tutorial

VLab4Mic, at its core, is designed as a collection of independent modules, each focused on a specific task within a typical imaging experiment,
from choosing a macromolecule of interest, deciding between strategies to label it, to simulating the detection of a virtual sample across imaging modalities.

Notheless, integrative modules connect these functionalities into a complex workflow that allow us to branch the simulation and analysis.

We also provide example models:
-   Structure:
    -   Nuclear Pore Complex
    -   Clathrin Coated Pit
    -   Matured HIV-capsid core
-   Probes
    -   Antibody
    -   Nanobody
    -   GFP
    -   Linker
    -   NHS-ester
    -   ... and more!
-   Imaging modalities
    -   Widefield
    -   Confocal
    -   AiryScan
    -   STED
    -   Single Molecule Localisation Microscopy

You can use VLab4Mic library through a jupyter notebook, with a graphical user interphase, or with a python stript according to your coding exprience.


# How to use VLab4Mic

VLab4Mic is a python library that can be used directly through our workflows and analysis module.
For a detail use see our [script](https://github.com/HenriquesLab/SupraMolecularSimulator/)

## Codeless Jupyter Notebooks:
| Category | Description | Last test | Notebook | Colab Link |
| --- | --- | --- | --- | --- |
| Main interphase | Create a virtual sample and simulate its image acquisition with multiple imaging modalities |  | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb)
| Parameter sweep | Generate and analyse simulations over a range of parameter combinatinos  |  | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb)


### 1.- Following the installation instructions, activate your python environment where you installed VLab4Mic.
### 2.- Lauch Jupyter Lab by running this command:
```shell
jupyter lab
```
This command will open a Jupyter Lab interphase in your web browser. Once here, locate the folder containing VLab4Mic notebooks. You can download them from [here](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/)


### 3.- Select the notebook and start using VLab4Mic!
If you have any doubt during the use of the notebook refer to the [code documentation](https://github.com/HenriquesLab/SupraMolecularSimulator/) or rise an [issue](https://github.com/HenriquesLab/SupraMolecularSimulator/issues/).

# Contributing
Contributions are very welcome. Please read our [Contribution Guidelines](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/CONTRIBUTING.md) to know how to proceed.

# Licence
Distributed under the terms of the [MIT license](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/LICENSE.txt), VLab4Mic is free and open source software

# Issues
Shoud you encounter any problem, do not hesitate on [letting us know](https://github.com/HenriquesLab/SupraMolecularSimulator/issues/). Don't forget to incude a detail description.

