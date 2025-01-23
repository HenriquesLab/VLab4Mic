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
It is advisable to create a separate environment to use SupraMolSim


To install latest development version:

```shell
pip install git+https://github.com/HenriquesLab/SupraMolecularSimulator.git
```

To install from source distribution
-   Download the distribution from [here](https://github.com/HenriquesLab/SupraMolecularSimulator/dist/)
-   Locate the folder containing the wheel and run 
```shell
pip install DIST_FILE_NAME.whl
```


# Tutorial

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

## Codeless SupraMolSim functionalities through graphical user interphase

Getting Started


## Codemuch SupraMolSim scripts


# Contributing

# Licence

# Issues
Shoud you encounter any problem, do not hesitate on [letting us know](https://github.com/HenriquesLab/SupraMolecularSimulator/issues/). Don't forget to incude a detail description.

