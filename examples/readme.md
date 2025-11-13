# VLab4Mic example scripts

To run the examples use this command:

```bash
python script_name.py
```

## Table of available scripts 

| Category | Description | File |
| --- | --- | --- |
| **Image virutal sample** | Create a virtual sample and simulate its imaging acquisition| [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/image_virtual_sample.py)|
| **Parameter sweep** | Set up and run a parameter sweep analysis | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/parameter_sweep.py)|
| **Real vs Simulation** | Create a simulated image based on a real experiment image | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/real_vs_simulatione.py)|
| **Custom structure** | Image a virutal sample based on a custom PDB/CIF file. | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/custom_structure.py)|

## Real vs Simulation

In this example, VLab4Mic creates a virtual sample where the position of the labelled particles are estimated from a real STED microscopy image.  The, we simulate the acquisition of this virtual sample with a model for STED microscopy. The script generates an image with a side to side comparison of the real image and the simulated image of the same modality.