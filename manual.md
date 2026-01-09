# VLab4Mic: A Virtual Laboratory for Microscopy

## What to expect from this manual

This manual covers minimal steps required for running imaging simulations and parameteres sweeps with VLab4Mic, either by using our codeless notebooks or using it as a python package.

### Usage Options

VLab4Mic main usage can be done in three ways, depending on your coding experience:
- **Jupyter Notebooks with GUI** - No coding required (recommended for beginners)
- **Jupyter Notebooks with code** - Moderate coding (recommended for intermediate users)
- **Python scripts** - Full coding control (recommended for advanced users)

This tutorial will cover **Jupyter Notebooks with GUI** and **Python Scripts**.


# VLab4Mic with Jupyter Notebooks with GUI

VLab4Mic codeless notebooks allow you to use its methods **without writing code**.
The following table lists our notebooks tailored to use VLab4Mic in Google Colab or in a local installation (your own computer, for instance).
 

## Table of available Notebooks 

| Category | Description | Notebook | Colab Link |
| --- | --- | --- | --- |
| **Main Interface** | Create virtual samples and simulate image acquisition with multiple imaging modalities | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_main.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/VLab4Mic/blob/main/notebooks/vLab4Mic_main.ipynb) |
| **Parameter Sweeps** | Generate and analyze simulations over parameter ranges for optimization | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) |

## Option 1: 📚 Use Jupyter Notebooks in Google Colab
Google Colab is an online service that provides computing resources to run Jupyter Notebooks without prior set up. This means that you don't need to install VLab4Mic in your computer, nor create virtual environments for it.  

To use VLab4Mic in Colab, click on the "Open in Colab" badge. This will open a new tab in your browser with the selected notebook. 

When the notebooks opens, follow the instructions in the notebook to get started!

### Troubleshooting:

If the Colab Link does not work, you can still use the codeless notebooks from your personal Google Colab by following these steps:
- Click on the "jupyter" badge of the notebook you want to use
- Download the notebook file by clicking on "Download raw file"
- Log in to Google Colab with your own Google account
- Upload a copy of the notebook from File > Upload Notebook 

The notebook should open once it finish uploading.

## Option 2: 📚 Use Jupyter Notebooks in a local installation

To use a notebook locally, you need to install VLab4Mic. Here's summary of the steps to follow:

- Install VLab4Mic in your local enviornment
- Download the python notebook
- Open Jupyter Lab in your local enviornment
- Locate and run your local copy of the notebook

Follow the nexts steps to set up VLab4Mic and use jupyter notebooks locally.

### Step 1️⃣: Create and activate a New Environment

VLab4Mic is compatible with **Python 3.9, 3.10, 3.11, and 3.12** on macOS, Windows, and Linux.

> 💡 **We strongly recommend** creating a new virtual environment to use VLab4Mic.  
> You can use either Conda ([Miniconda installation](https://docs.conda.io/en/latest/miniconda.html)) or Python's built-in `venv` module.

Click in the arrows to expand the instructions for using either conda or venv.


#### Create the virtual environment
<details> 
<summary> 🐍 Using Conda </summary> 

Run the following command (replace `myenv` with your desired environment name):

```bash
conda create --name myenv python=3.11
```

Then activate it:
```bash
conda activate myenv
```
</details>

<details>
<summary> 🐍 Using venv (built-in Python module) </summary> 

Run the following commands (replace `myenv` with your desired environment name):

```bash
python -m venv myenv
```

Activate the environment:

- On macOS/Linux:
    ```bash
    source myenv/bin/activate
    ```
- On Windows:
    ```bash
    myenv\Scripts\activate
    ```
</details>


> Make sure to activate your environment before continuing to the next step.


### Step 2️⃣: Install VLab4Mic

Run the following command to install vlab4mic with the necesary libraries to support our Jupyter notebooks:

```bash
pip install vlab4mic vlab4micjupyter
```


### Step 3️⃣: Launch Jupyter Lab
Run the next command to launch jupyter lab.

```bash
jupyter lab
```
This will open Jupyter Lab in your web browser.

### Step 4️⃣ Download and Open Notebooks

Download the notebooks from our [repository here](https://github.com/HenriquesLab/VLab4Mic/tree/main/notebooks), or refer to the table above.  
Once downloaded, use the Jupyter Lab interface to find and open the notebook.

### Step 5️⃣: Start Experimenting!

Follow the instructions in each notebook.  
For questions, refer to our [documentation](https://github.com/HenriquesLab/VLab4Mic/wiki) or [create an issue](https://github.com/HenriquesLab/VLab4Mic/issues).


# VLab4Mic as Python Library (scripts) 

When using VLab4Mic as a Python library, you can parameterize your experiment or parameter sweep through its high-level functions. These functions are covered in detail in the Methods section.

Here's summary of the steps to follow for using VLab4Mic as a python library:

- Create a virtual environment
- Install VLab4Mic in your virtual environment
- Start using VLab4Mic!

Follow the following steps to set up VLab4Mic.

### Step 1️⃣: Create and activate a New Environment

VLab4Mic is compatible with **Python 3.9, 3.10, 3.11, and 3.12** on macOS, Windows, and Linux.

> 💡 **We strongly recommend** creating a new virtual environment to use VLab4Mic.  
> You can use either Conda ([Miniconda installation](https://docs.conda.io/en/latest/miniconda.html)) or Python's built-in `venv` module.


Click in the arrows to expand the instructions for using either conda or venv.

#### Create the virtual environment
<details> 
<summary> 🐍 Using Conda </summary> 

Run the following command (replace `myenv` with your desired environment name):

```bash
conda create --name myenv python=3.11
```

Then activate it:
```bash
conda activate myenv
```
</details>

<details>
<summary> 🐍 Using venv (built-in Python module) </summary> 

Run the following commands (replace `myenv` with your desired environment name):

```bash
python -m venv myenv
```

Activate the environment:

- On macOS/Linux:
    ```bash
    source myenv/bin/activate
    ```
- On Windows:
    ```bash
    myenv\Scripts\activate
    ```
</details>

> Make sure your environment is active before continuing to the next step.


### Step 2️⃣: Install VLab4Mic

Run the following command to install vlab4mic:

```bash
pip install vlab4mic
```

Alternatively, run the next command to include necesary dependencies to support jupyter notebooks as well:

```bash
pip install vlab4mic vlab4micjupyter
```

### Step 3️⃣: Start using VLab4Mic!


### 🖥️ Option 1: Run python scripts from the command line

You can run python scripts from the command line in your virtual envornment.

For detailed usage examples, see our [example scripts](https://github.com/HenriquesLab/VLab4Mic/tree/main/examples).

Here, we exemplify the case of running a script named "vlab_script.py" located in the current directory.

```python 
python vlab_script.py
```



The contents of the script can be tailored for your needs. Here's a quick example of such a script that runs an imaging simulation with the function "image_vsample":

```python
# contents of vlab_script.py
from vlab4mic.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()
```

### 🖥️ Option 2: Use VLab4Mic in an interactive python interpreter
Start a python interpreter by runing "python" from the command line with your virtual environment.

Start using VLab4Mic through the interpreter, for instance:

```bash
>>> from vlab4mic.experiments import image_vsample
>>> image_outputs, image_noiseless_outputs, experiment = image_vsample()
```

Get a description of the available parameters of each method by using Python's built-in `help()`:

```python
help(image_vsample)
```

# VLab4Mic Methods 
This section of the tutorial describes VLab4Mic high-level functions to:
- 🧬 **Create and simulate imaging of a virtual sample**
- 🔬 **Set up and run parameter sweep**

The Table of Available Notebooks contains a notebook focused on each high-level method. The examples that follow describe VLab4Mic usage as a python package.

## 🧬 Create and image a virtual sample

VLab4Mic allows you to create a virtual sample and generate imaging simulations of it on one or many modalities through the method "image_vsample".

```python
from vlab4mic.experiments import image_vsample

image_outputs, image_noiseless_outputs, experiment = image_vsample()
```

The example above creates a default virtual sample and returns imaging simulations with and without noise.  
The `image_vsample` method also returns an experiment object containing all modules used to generate the simulations.

### ⚙️ Parameterize your virtual sample and imaging

The `image_vsample` method allows you to specify each step of the virtual sample, such as:
- Structure of interest
- Type of probe (e.g. Antibody, Linker, GFP, ...)
- Target site on the structure (epitope sequence or residues)
- Labelling efficiency of probe
- Distance of probe to epitope
- Number of labelled particles
- Sample dimensions
- ...and more

Example:

```python
from vlab4mic.experiments import image_vsample

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
```

> 📖 Find a complete and detailed list of parameters in the documentation of the method, or refer to the Parameter Tables section.

We also provide probe templates tailored for specific structures, such as direct labelling of Nup96 protein of the Nuclear Pore Complex (NPC):

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(
    structure="7R5K",  # PDB ID code for a Nuclear Pore complex
    probe_templates=["NPC_Nup96_Cterminal_direct",],  # Pre-set probe for 7R5K
    run_simulation=True,
    clear_experiment=True,
)
```

A table of these and other templates can be found in the Templates section.

Run the following example to display the results of an imaging simulation:

```python
import matplotlib.pyplot as plt
nmods = len(modalities)
fig, axs = plt.subplots(1, nmods)
nframe = 0
for i, mod in enumerate(modalities):
    axs[i].imshow(images[mod][nframe], cmap="magma")
    axs[i].set_title(mod)
```


## 🔬 Run a parameter sweep analysis

VLab4Mic allows you to test parameter combinations for sample and image simulations.

The main method to run a parameter sweep is `run_parameter_sweep`. This is an example of its usage with default parameters:

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

This generates and saves a default parameter sweep. The `sweep_gen` object contains all parameters used to set and run the sweep.

##  Specify values for each parameter

You can parameterize the sweep by specifying the values or ranges to use for each parameter.  
You can parameterise by passing the values:
- as a list,
- as a tuple with the format (min, max, nsteps) to generate linearly spaced values.

Example:

```python
from vlab4mic.sweep_generator import run_parameter_sweep

sweep_gen = run_parameter_sweep(
    structures=["7R5K",],
    probe_templates=["NPC_Nup96_Cterminal_direct",],  # Probe template tailored to 7R5K
    sweep_repetitions=20,
    # parameters for sweep
    labelling_efficiency=(0, 1, 5),  # 5 linearly spaced values between 0 and 1
    defect=(0, 1, 5),  # 5 linearly spaced values between 0 and 1
    defect_small_cluster=[300,],  # 1 single value 
    defect_large_cluster=[600,],  # 1 single value 
    exp_time=[0.001, 0.01,],  # 2 values
    # output and analysis
    output_name="vlab_script",
    return_generator=True,
    save_sweep_images=True,  # By default, the saving directory is set to the home path of the user
    save_analysis_results=True,
    run_analysis=True
)
```

> ⚠️ **Note:** When running a parameter sweep, all possible parameter combinations will be used.  
> For example, setting 10 values for labelling efficiencies and 10 values for particle defects will generate 100 combinations.

### 📝 Available parameters for sweeping

Each parameter in the following list is available for sweep (one or more values passed).  
If a parameter is `None`, it will not be swept and will use default values.

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



# 📋 Parameter tables

## 🏗️ Structure parameters

| Parameter name | Description | Notes |
| --- | --- | --- |
| **structure** | 4-letter ID for a structural model from PDB database, or path to a structure file |  |  

## 🧪 Probe parameters

| Parameter name | Description |
| --- | --- | 
| **probe_template** | Name of a pre-set probe configuration file (filename) |
| **probe_name** | Optional. Alternative name to identify this probe |
| **labelling_efficiency** | Efficiency of probe labelling (probability of probe to bind a given epitope). Float, default 1.0 |
| **probe_target_type** | Type of target: "Sequence", "Residue", or "Primary" |
| **probe_target_value** | Specific value depending on the target type (e.g., sequence string, residue info, etc.) |
| **probe_target_option** | Additional option for the probe target, used for secondary epitopes |
| **probe_distance_to_epitope** | Distance from the probe to the epitope (float, minimal distance set from epitope and probe paratope) |
| **probe_model** | 4-letter code for an atomic model on which to base this probe |
| **probe_fluorophore** | Fluorophore name for emitters (e.g., "AF647") |
| **probe_paratope** | If probe_model is defined, the paratope defines the anchor point of the probe |
| **probe_conjugation_target_info** | If probe_model is defined, this dictionary specifies the sites to use as probe emitters |
| **probe_conjugation_efficiency** | Efficiency of conjugation of emitters (float) |
| **probe_seconday_epitope** | If probe is secondary, this amino acid sequence defines the epitope on the primary antibody model |
| **probe_wobble_theta** | Enable probe wobbling (float or None) |
| **probe_steric_hindrance** | Steric hindrance value or configuration (distance between epitopes) |
| **peptide_motif** | Dictionary specifying motif extraction for probe target sequence |
| **as_primary** | Whether to treat the probe as a primary linker (bool) |

## 🧩 Defects parameters

| Parameter name | Description |
| --- | --- | 
| **defect** | Fraction of the particle rendered inaccessible to probes (float) |
| **defect_small_cluster** | Maximum distance between epitopes for first grouping (in Å) |
| **defect_large_cluster** | Minimum distance between epitopes for second grouping (in Å) |

## 🧱 Virtual sample parameters

| Parameter name | Description |
| --- | --- | 
| **virtual_sample_template** | Name of the configuration file for the virtual sample template |
| **sample_dimensions** | In nanometers, define the X, Y, and Z sizes of the field (list of float) |
| **number_of_particles** | Number of independent copies of a particle to create and distribute (int) |
| **particle_positions** | List of relative positions of particles in the field (list of np.array or list) |
| **random_orientations** | If True, each particle will be randomly assigned a new orientation (bool, default False) |
| **xy_orientations** | List of orientation directions on the xy plane to choose from |
| **xz_orientations** | List of orientation directions on the xz plane to choose from |
| **yz_orientations** | List of orientation directions on the yz plane to choose from |
| **axial_offset** | List of Z positions to choose from as offsets from the focus plane |
| **random_placing** | Define if position in field is random or the center of field (bool, default False) |
| **random_rotations** | If True, apply random rotations to particles (bool, default False) |
| **rotation_angles** | List of rotation angles to choose from |
| **minimal_distance** | Minimal allowed distance between particles (float) |

## 🔬 Modalities parameters

| Parameter name | Description | 
| --- | --- | 
| **pixelsize_nm** | Pixel size of image in nanometers (float or int) |
| **lateral_resolution_nm** | Lateral resolution in nanometers (float or int) |
| **axial_resolution_nm** | Axial resolution in nanometers (float or int) |
| **psf_voxel_nm** | Voxel size used to render the PSF in nanometers (float or int, defaults to 10 nm) |
| **depth_of_field_nm** | Depth of field in nanometers (float or int) |

## 📸 Acquisition parameters

| Parameter name | Description | 
| --- | --- | 
| **exp_time** | Exposure time in seconds (float) |
| **noise** | Whether to use noise at detection or not (bool) |
| **nframes** | Number of frames (int) |

---

# 5. 🗂️ Templates

## 🏗️ Structures

Any atomic structure in PDB/CIF format can be used in VLab4Mic.  
Below are example structures readily available for codeless notebooks and for which example probes exist.

| Structure id | Description | Database | Format used |
| --- | --- | --- | --- |
| 1XI5 | Clathrin D6 coat with auxilin J-domain | RCSB | CIF format |
| 7R5K | Human nuclear pore complex (constricted) | RCSB | CIF format |
| 3J3Y | Atomic-level structure of the entire HIV-1 capsid (186 hexamers + 12 pentamers) | RCSB | CIF format |
| 8GMO | Bacteriophage T4 capsid shell  | RCSB | CIF format |

## 🧪 Probes

Each probe is a configuration file whose parameterization models a particular type of labelling.  
For instance, NHS_ester models a direct label that targets lysine residues on any atomic model.  
Indirect probes include an atomic model for the probe itself, as is the case for GFP or an Antibody.

### Structure-independent probes

| Probe name | Description | Targets | Probe Model | Notes |
| --- | --- | --- | --- | --- |
| **NHS_ester** | Direct label for Lysine residues | Residues: LYS | NA | The location of each lysine is used as position for emitters |
| **Linker** | Indirect label modelled as a rigid spacer that binds to the epitope | NA | NA | Requires specification of target type and values |
| **Antibody** | Indirect label based on the crystallographic structure of IGG antibody against HIV-1 isolates | NA | 1HZH | Requires specification of target type and values |
| **Nanobody** | Indirect label based on the crystallographic structure of a synthetic nanobody | NA | 7MFV | Requires specification of target type and values |
| **GFP** | Indirect label based on the crystallographic structure of GFP | NA | 6YLQ | Requires specification of target type and values |
| **GFP_w_nanobody** | Indirect label based on the crystallographic structure of a nanobody in complex with eGFP | NA | 6XZF | Requires specification of target type and values |
| **mMaple** | Indirect label based on the crystallographic structure of the fluorescent protein mTFP1 | NA | 2HQK | Requires specification of target type and values |
| **SNAP-tag** | Indirect label based on the crystallographic structure of SNAP-tag with labelled with a fluorophore | NA | 6Y8P | Requires specification of target type and values |

### Structure-specific probes

| Probe name | Description | Targets | Probe Model | Notes |
| --- | --- | --- | --- | --- |
| **CCP_heavy_chain_Cterminal** | Direct probe for C terminal site of Heavy Chain of Clathrin Coated Pit | Sequence: EQATETQ | NA | Based on model 1XI5 |
| **NPC_Nup96_Cterminal_direct** | Direct probe C terminal site of Nup96 in Nuclear Pore Complex | Sequence: ELAVGSL | NA | Based on model 7R5K |
| **HIV_capsid_p24_direct** | Direct probe for HIV-1 Capsid monomer (p24) | Sequence: SPRTLNA | NA | Based on model 3J3Y |
| **anti-p24_primary_antibody_HIV** | Antibody for HIV-1 Capsid monomer (p24) | Sequence: SPRTLNA | 1HZH | Based on model 3J3Y |

## 🧱 Virtual sample

| Virtual sample name | Description | Notes |
| --- | --- | --- |
| **square1x1um_randomised** | Virtual sample of square dimensions and a single particle placed in the center |  |

## 🔬 Modalities

| Modality name | Description | PSF shape X,Y,Z (nm) | Image pixelsize (nm) | Notes |
| --- | --- | --- | --- | --- |
| **Widefield** | Widefield microscope | 94,94,331 | 100 | |
| **Confocal** | Confocal microscope | 94,94,331 | 70 | |
| **STED** | STED microscope | 20,20,20 | 15 | |
| **SMLM** | SMLM image model | 8,8,8 | 5 | Emulates effective image |
| **Reference** | Idealised microscope | 5,5,5 | 5 | |
