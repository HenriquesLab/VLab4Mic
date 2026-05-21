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

## Article Figure Scripts

Scripts used to generate the figures in the VLab4Mic publication. Each script is self-contained and demonstrates a specific feature or validation scenario.

| Figure | Description | Script |
| --- | --- | --- |
| **Fig. 1** | HIV capsid (3J3Y) imaged across all modalities | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/fig1_hiv_examples.py) |
| **Fig. 2** | NPC (7R5K) panel with SMLM and quantification | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/fig2_npc_panel.py) |
| **Fig. 3 (A–C)** | Parameter sweep panels for structural comparison | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/fig3_abc_panels.py) |
| **Fig. 3 (D–E)** | Parameter sweep with labelling efficiency | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/fig3_de_panels.py) |
| **Fig. 4** | Clathrin-coated structure (CCS) imaging | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/fig4_ccs.py) |
| **Fig. S — AlphaFold** | Using a locally downloaded AlphaFold model as input structure | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_alphafold_from_localfile.py) |
| **Fig. S — Depth of field** | Effect of depth-of-field parameter on simulated images | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_depth_of_field.py) |
| **Fig. S — Image positioning** | Using a reference image to position epitopes in the virtual sample | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_image_for_positioning.py) |
| **Fig. S — Labelling efficiency & hindrance** | Labelling efficiency and steric hindrance effects | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_LabEff_hindrance.py) |
| **Fig. S — SMLM localisation (CCS)** | Localisation table simulation for clathrin-coated structures | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_localisation_simulation_ccs.py) |
| **Fig. S — Multicolour** | Two-colour imaging of the same structure with different probes | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_multicolour.py) |
| **Fig. S — Multiple structures** | Imaging multiple different structures in the same field of view | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_multiple_structures.py) |
| **Fig. S — NPC quantification** | NPC ring quantification simulation | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_NPC_quantification.py) |
| **Fig. S — Primary & secondary** | Indirect labelling with primary and secondary antibody models | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_primaryandsecondary.py) |
| **Fig. S — Probe examples** | Gallery of all probe templates across different structures | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_probe_examples.py) |
| **Fig. S — Site-specific labelling** | Site-specific probe conjugation on a defined residue | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_site_specific_labelling.py) |
| **Fig. S — Structural integrity** | Structural integrity parameter effects on NPC and capsid | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_structural_integrity.py) |
| **Fig. S — Structure normals** | Structure normal directions and probe positioning | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_structure_normals_and_probe_positioning.py) |
| **Fig. S — Structure parsing** | Structure parsing, asymmetric units, and multiple targets | [![script](https://img.shields.io/badge/script-grey)](https://github.com/HenriquesLab/VLab4Mic/blob/main/examples/article_figures/figS_structure_parsing_and_targets.py) |

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
