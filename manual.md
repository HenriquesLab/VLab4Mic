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

When using VLab4Mic as a python library, you can paramtererise your experiment or parameter sweep through its high level functions. These functions are covered in detail in the Methods section.


Following the installation instructions above, activate your Python environment.


Access VLab4Mic main methods such as image_vsample as follows:

```python
from vlab4mic.experiments import image_vsample

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
from vlab4mic.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()

```

The former example will first create a default virtual sample and return imaging simulations of this sample with and without noise. Besides the imaging simulations,  "image_vsample" method returns an experiment object containing all the modules used to generate the simulations.


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

The following code example shows how to use some of these parameters.

```python
from vlab4mic.experiments import image_vsample

# Paramterise and run simulation
modalities = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

images, noiseless, experiment = image_vsample(
    structure="7R5K", # PDB ID code for a Nuclear Pore complex
    probe_template="Antibody", # Probe template for an antibody
    probe_target_type="Sequence",  
    probe_target_value="ELAVGSL",  # epitope sequence
    number_of_particles = 10,
    random_rotations=True,
    rotation_angles=None,
    multimodal=modalities,
    run_simulation=True,
    clear_experiment=True,
)
```
Find a complete and detailed list of parameters in the documentation of the method, or refer to the table of parameters section.

We also provide probe templates tailored for specific structures, such as the following example:

```python
from vlab4mic.experiments import image_vsample

# Paramterise and run simulation
images, noiseless, experiment = image_vsample(
    structure="7R5K", # PDB ID code for a Nuclear Pore complex
    probe_templates=["NPC_Nup96_Cterminal_direct",], # Pre-set probe for 7R5K
    run_simulation=True,
    clear_experiment=True,
)
```

A table of these and other templates can be found in the Templates tables.


Imaging simulations output is a dictionary where each key is the imaging modality. To display the results you can use the following example code:

```python
import matplotlib.pyplot as plt
nmods = len(modalities)
fig, axs = plt.subplots(1, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][nframe], cmap="magma")
    axs[i].set_title(mod)
```


## 3.2 Run a parameter sweep analysis

vLab4Mic allows you to test parameter combinations to use for sample and image simulations. <br>
The main method to run a parameter sweep is "run_parameter_sweep" and can be used as in the following example:

```python
from vlab4mic.sweep_generator import run_parameter_sweep

sweep_gen = run_parameter_sweep(
    # output and analysis
    output_name="vlab_script",
    return_generator=True,
    save_sweep_images=True,
    save_analysis_results=True,
    run_analysis=True
    )
```

This example will generate and save a default parameter sweep. The "sweep_gen" object contains all parameters used to set and run the sweep. 

## Specify values for each parameter 
Similar to "image_vsample" method, you can parameterise the sweep by specifying the values or ranges to use on each parameter. To do so, pass a list of the values or a tuple with (min, max, nsteps) to generate linearly spaced values.

For instance, you can set to use 10 equally distributed values for Labelling Efficiencies between 0 and 1; or you can set them explicitly: 0.4, 0.8 and 1.0. Each parameter has a default value.


```python
from vlab4mic.sweep_generator import run_parameter_sweep

sweep_gen = run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],  # Probe template tailored to 7R5K
    sweep_repetitions=20,
    # parameters for sweep
    labelling_efficiency=(0,1,5),  # 5 linearly spaced values between 0 and 1
    defect=(0,1,5), # 5 linearly spaced values between 0 and 1
    defect_small_cluster=[300,],    # 1 single value 
    defect_large_cluster=[600,],    # 1 single value 
    exp_time=[0.001, 0.01,],    # 2 values
    # output and analysis
    output_name="vlab_script",
    return_generator=True,
    save_sweep_images=True,  # By default, the saving directory is set to the home path of the user
    save_analysis_results=True,
    run_analysis=True
    )
```
**Note: when running a parameter sweep, all possible parameter combinations will be used. For instance, setting 10 values for Labelling Efficiencies and 10 values for particle defects will generate 100 combinations.**

## Available parameters for sweeping

Each parameter in the following list is available for sweep (one or more values passed).
Note: If a parameter is None, it will not be swept and will use default values.
- probe_target_type
- probe_target_value
- probe_target_option
- probe_model
- probe_fluorophore
- probe_paratope
- probe_conjugation_target_info
- probe_seconday_epitope
- peptide_motif
- probe_distance_to_epitope
- probe_steric_hindrance
- probe_conjugation_efficiency
- probe_wobble_theta
- labelling_efficiency
- defect
- defect_small_cluster
- defect_large_cluster
- sample_dimensions
- particle_orientations
- rotation_angles
- minimal_distance
- pixelsize_nm
- lateral_resolution_nm
- axial_resolution_nm
- psf_voxel_nm
- depth_of_field_nm
- exp_time




# 4. Parameter tables

# Structure parameters
| Parameter name | Description | Notes |
| --- | --- | --- |
| **structure** | 4-letter ID for a structural model from PDB database|  |  

# Probe parameters
| Parameter name | Description |
| --- | --- | 
|**probe_template**| Name of a pre-set probe |
|**probe_name**| Optional. Alternative name to identify this probe |  
|**labelling_efficiency**| Probability of probe to bind a given epitope | 
|**probe_target_type** | Type of target: "Sequence", "Residue" or "Primary"|
|**probe_target_value**  | Specific value depending on the target type |  
|**probe_distance_to_epitope** | Fix distance offset between epitope and probe anchor site |    
|**probe_model** | 4-letter code for an atomic model on which to base this probe |   
|**probe_fluorophore** | Fluorophore name for emitters |    
|**probe_paratope** | If probe_model is defined, the paratope defines the anchor point of the probe |    
|**probe_conjugation_target_info** | If probe_model is defined, this dictionary specifies the sites to use as probe emitters  |    
|**probe_seconday_epitope** | If probe is secondary, this aminoacid sequence defines the epitope on the primary antibody model|    
|**probe_wobbling** | Flag to specify is proble wobbling is allowed |


# Defects parameters
| Parameter name | Description |
| --- | --- | 
| **defect** | Average fraction of the particle rendered inaccessible to probes |  
| **defect_small_cluster** | Maximum distance between epitopes for first groupping |
| **defect_large_cluster** | Minimum distance between epitopes for second groupping |   



# Virtual sample parameters
| Parameter name | Description |
| --- | --- | 
|**virtual_sample_template** | |  
|**sample_dimensions** | List of XYZ dimensions for the virtual sample  |  
|**number_of_particles** | Number of particles to place in the sample space |    
|**particle_positions** | List of relative positions per particle |   
|**random_orientations** |  |
|**random_placing** |  |  


# Modalities parameters
| Parameter name | Description | 
| --- | --- | 
| **pixelsize** | Pixelsize of image |   
| **lateral_resolution_nm** | lateral resolution |   
| **axial_resolution_nm** | axial resolution nm |   
| **psf_voxel_nm** | Voxel size use to render the PSF. Defualts to 10 nm. Smaller voxel size will increase computation time |   
| **depth_of_field_nm** | |  

# Acquisition parameters
| Parameter name | Description | 
| --- | --- | 
| **exp_time** | Exposure time in seconds | 
| **noise** | Whether to use noise at detection or not | 
| **nframes** | Number of frames | 




# 5. Templates


# Structures


Any atomic structure in PDB/CIF format can be use in vLab4Mic, nonetheless, this table list example structures readlily available for codeless notebooks and for which example probes exist.

| Structure id | Description | Database | Format used |
| --- | --- | --- | --- |
| 1XI5 | Clathrin D6 coat with auxilin J-domain | RCSB | CIF format |
| 7R5K | Human nuclear pore complex (constricted) | RCSB | CIF format |
| 3J3Y | Atomic-level structure of the entire HIV-1 capsid (186 hexamers + 12 pentamers) | RCSB | CIF format |
| 8GMO | Bacteriophage T4 capsid shell  | RCSB | CIF format |

# Probes


Each probe is a configuration file whose parameterisation models a particular type of labelling. For instance, NHS_ester models a direct label that targets lysine residues on any atomic model. In the case of indirect probes, they include an atomic model for the probe itself, as is the case for GFP or an Antibody.

## Structure-independent probes

Structure-independent probes are model probes whose target epitope is not pre-defined or is a residue.

| Probe name | Description | Targets | Probe Model | Notes |
| --- | --- | --- | --- | --- |
| **NHS_ester** | Direct label for Lysine resiudes | Residues: LYS | NA | The location of each lysine is used as position for emitters |
| **Linker** | Indirect label modelled as a rigid spacer that binds to the epitope | NA | NA | Requires specification of target type and values|
| **Antibody** | Indirect label based on the crystallographc structure of IGG antibody against HIV-1 isolates | NA | 1HZH | Requires specification of target type and values|
| **GFP** | Indirect label based on the crystallographc structure of  GFP | NA | 6YLQ | Requires specification of target type and values|

## Structure-specific probes

This probes are special cases of structure-independent probe models which aim to target a specific structure.


| Probe name | Description | Targets | Probe Model | Notes |
| --- | --- | --- | --- | --- |
| **CCP_heavy_chain_Cterminal** | Direct probe for C terminal site of Heavy Chain of Clathin Coated Pit | Sequence: EQATETQ | NA | Based on model 1XI5 |
| **NPC_Nup96_Cterminal_direct** | Direct probe C terminal site of Nup96 in Nuclear Pore Complex | Sequence: ELAVGSL | NA | Based on model 7R5K |
| **HIV_capsid_p24_direct** | Direct probe for HIV-1 Capsid monomer (p24) | Sequence: SPRTLNA | NA | Based on model 3J3Y |
| **anti-p24_primary_antibody_HIV** | Antibody for HIV-1 Capsid monomer (p24) | Sequence: SPRTLNA | 1HZH | Based on model 3J3Y |






# Virtual sample

| Virtual sample name | Description | Notes |
| --- | --- | --- |
| square1x1um_randomised | Virtual sample of square dimensions and a single particle placed in the center |  |

# Modalities

| Modality name | Description | PSF shape X,Y,Z (nm) | Image pixelsize (nm) | Notes |
| --- | --- | --- | --- | --- |
| **Widefield** | Widefield microscope | 94,94,331 | 100 | |
| **Confocal** | Confocal microscope | 94,94,331 | 70 | |
| **STED** | STED microscope | 20,20,20 | 15 | |
| **SMLM** | SMLM image model | 8,8,8 | 5 | Emulates effective image |
| **Referece** | Idealised microscope | 5,5,5 | 5 | |
