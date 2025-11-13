# VLab4Mic example scripts

To run the examples use this command:

```bash
python script_name.py
```

## Table of available scripts 

| Category | Description | File |
| --- | --- | --- |
| **Real vs Simulation** | Create a simulated image based on a real experiment image | [![script](https://img.shields.io/badge/script-grey)]()|

## Real vs Simulation

In this example, VLab4Mic creates a virtual sample where the position of the labelled particles are estimated from a real STED microscopy image.  The, we simulate the acquisition of this virtual sample with a model for STED microscopy. The script generates an image with a side to side comparison of the real image and the simulated image of the same modality.