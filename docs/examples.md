# Examples

Ready-to-run example scripts demonstrating common VLab4Mic workflows. Before running, make sure you have [installed VLab4Mic](getting-started.md).

Run any script with:

```bash
python path/to/script.py
```

---

## Available Scripts

| Category | Description | Script |
| --- | --- | --- |
| **Image virtual sample** | Create a virtual sample and simulate its imaging acquisition | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/image_virtual_sample.py) |
| **Parameter sweep** | Set up and run a parameter sweep analysis | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/parameter_sweep.py) |
| **Real vs Simulation** | Create a simulated image based on a real experimental image | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/real_vs_simulation.py) |
| **Custom structure** | Image a virtual sample using a custom local PDB/CIF file | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/custom_structure.py) |

---

## Image Virtual Sample

Demonstrates the core `image_vsample` function. In this example, 10 Nuclear Pore Complexes (NPC, model 7R5K) are labelled at the C-terminal of Nup96 and imaged across multiple modalities.

Key things shown:
- Selecting a structural model from PDB
- Using a pre-configured probe template
- Simulating acquisition across Widefield, Confocal, STED, and SMLM

[View script on GitHub →](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/image_virtual_sample.py)

---

## Parameter Sweep

Demonstrates `run_parameter_sweep`. In this example, labelling efficiency and structural integrity are swept across 3 values each (9 combinations total), and each output is compared against a reference image via Pearson correlation.

Key things shown:
- Configuring sweep parameters as lists or ranges
- Running repeated simulations
- Generating heatmap analysis outputs

[View script on GitHub →](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/parameter_sweep.py)

---

## Real vs Simulation

Shows how to base a simulated image on the parameters estimated from a real experimental acquisition — useful for validation and benchmarking.

[View script on GitHub →](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/real_vs_simulation.py)

---

## Custom Structure

Shows how to load and use a local PDB or CIF file as the structural model, instead of downloading from the PDB database.

[View script on GitHub →](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/custom_structure.py)
