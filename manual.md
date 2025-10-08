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


# 2. How to use VLab4Mic

VLab4Mic can be used as a Python library or through codeless jupyter notebooks.
For detailed usage examples, see our [example scripts](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/examples).

## 2.1 Codeless jupyter notebooks

VLab4Mic codeless notebooks allows you to use its methods without the need to code or modify scripts. 

## Where to find and use VLab4Mic codeless notebooks

You can use a notebook on a local installation or on Google Colab. 
For using VLab4Mic in colab, click on the "Open in Colab" and follow the instructions on the Notebook.

| Category | Description | Notebook | Colab Link |
| --- | --- | --- | --- |
| **Main Interface** | Create virtual samples and simulate image acquisition with multiple imaging modalities | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_main.ipynb) |
| **Parameter Sweeps** | Generate and analyze simulations over parameter ranges for optimization | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/SupraMolecularSimulator/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) |

## Use codeless notebooks on a local installation

To use a notebook locally, you need to install VLab4Mic as detailed in the Installation section.

### Step 1: Activate your environment
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
Download the notebooks from our [repository here](https://github.com/HenriquesLab/SupraMolecularSimulator/tree/main/notebooks), or refer to the table above.
Once downloaded, use the Jupyter Lab interphase to find the notebook and open it.

### Step 4: Start Experimenting!
Follow the interactive widgets and instructions in each notebook. For questions, refer to our [documentation](https://github.com/HenriquesLab/SupraMolecularSimulator/wiki) or [create an issue](https://github.com/HenriquesLab/SupraMolecularSimulator/issues).

## 2.2 Using VLab4Mic as Python library

When using VLab4Mic as a python library, you can paramtererise your experiment or parameter sweep through its high level functions. These functions are covered in dettail in the Methods section.


Following the installation instructions above, activate your Python environment.


Access VLab4Mic main methods such as image_vsample as follows:

```python
from supramolsim.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()

```

Get a description of the available parameters of each method by using the built-in python method "help()":

```python
help(image_vsample)
```

# 3. VLab4Mic Methods

vLab4Mic main functionallity can be accessed with 2 main methods.
- Create and simulate imaging of a virtual sample
- Set up and run parameter sweep


## 3.1 Create and image a virtual sample

VLab4Mic allows you to create a virtual sample and generate imaging simulations of it on one or many modalities. <br> To do so, use "image_vsample" method as follows:

```python
from supramolsim.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()

```
This example will first create a default virtual sample and return imaging simulations of this sample with and without noise. Besides the imaging simulations,  "image_vsample" method returns an experiment object containing all the modules used to generate the simulations.

## Parameterise your virtual sample and imaging

VLab4Mic method "image_vsample" allows you to specify each step of the virtual sample such as:
-   Structure of interest
-   Type of probe (e.g. Antibody, Linker, GFP,...)
-   Target site on the structure (Epitope sequence, or residues)
-   Labelling efficiency of probe
-   Distance of probe to epitope
-   Number of labelled particles
-   Sample dimensions
-   ...

Find a complete and detailed list of parameters in the documentation of the method, or refer to the [table of parameters]()

```python
from supramolsim.experiments import image_vsample
# parameters
structure = "7R5K" # PDB ID code for a Nuclear Pore complex
probe_template = "Antibody" # example probe template name
modalities = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

# Paramterise and run simulation
images, noiseless, experiment = image_vsample(
    structure=structure,
    probe_template=probe_template,
    number_of_particles = 10,
    random_rotations=True,
    rotation_angles=None,
    multimodal=modalities,
    run_simulation=True
    clear_experiment=True,
)
```


## Parameters table:


### **<ins>Structure parameters:</ins>** 
- **Structure**: The 4 Letter code that identifies the atomic structure, for instance, "7R5K"

### **<ins>Probes and probe parameters:</ins>**
- **Probe template**: Name for a template probe such as “Antibody”, “NHS ester”, etc. This parameter defines the model to be used as probe. Other optional parameters control specific behaviour of the probe
- **Labelling efficiency**: Efficiency of probe to targeting an epitope in the structure
- ... 

### **<ins>Particle field parameters:</ins>**
- **Number of particles**: Number of labelled particles to be distributed in the sample volume.
- **Sample dimensions**: Dimensions of the volume for the virtual sample
- ... 

The full list of parameters and their in-depth usage can be found [here](https://github.com/jdmartinez24/vlabwikis/wiki/Parameter-tables).



## 3.2 Run a parameter sweep analysis

vLab4Mic allows you to test parameter combinations to use for sample and image simulations. <br>
The main method to run a parameter sweep is "run_parameter_sweep" and can be used as in the following example:


```python
from vlab4mic import sweep_generator

sweep_gen = sweep_generator.run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],
    sweep_repetitions=20,
    # parameters for sweep
    labelling_efficiency=(0,1,5),
    defect=(0,1,5),
    defect_small_cluster=[300,],
    defect_large_cluster=[600,],
    exp_time=[0.001, 0.01,],
    # output and analysis
    output_name="vlab_script",
    return_generator=True,
    save_sweep_images=True,
    save_analysis_results=True,
    run_analysis=True
    )
```

From "run_parameter_sweep" it is possible to parameterise the while parameter sweep. All parameters available for this method can be found in the Parameter sweep section.


### **<ins>Setting up parameter values</ins>**


In order to specify the parameters to include and the values to iterate in two ways. 
For instance, you can set to use 10 equally distributed values for Labelling Efficiencies between 0 and 1; or you can set them explicitly: 0.4, 0.8 and 1.0, for instance. Each parameter has a default value.

When running a parameter sweep, all possible parameter combinations will be used. For instance, setting 10 values for Labelling Efficiencies and 10 values for particle defects will generate 100 combinations. 
