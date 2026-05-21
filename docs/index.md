---
hide:
  - toc
---

<div class="hero" markdown>

# VLab4Mic

**A Virtual Laboratory for Microscopy**

Simulate fluorescence microscopy experiments — build virtual samples, apply fluorescent labelling, and image across modalities — before stepping into the microscope room.

[Get Started](getting-started.md){ .md-button .md-button--primary }
[View on GitHub](https://github.com/HenriquesLab/VLab4Mic){ .md-button }

</div>

---

## Why VLab4Mic?

<div class="grid cards" markdown>

-   :material-dna:{ .lg .middle } **Virtual Samples from Atomic Structures**

    ---

    Build virtual samples directly from PDB/CIF structural models and populate fields of view with particles.

-   :material-tag-multiple:{ .lg .middle } **Realistic Fluorescent Labelling**

    ---

    Model direct and indirect labelling with antibodies, nanobodies, GFP, SNAP-tags and more — including labelling efficiency and steric hindrance.

-   :material-microscope:{ .lg .middle } **Multi-Modality Imaging Simulation**

    ---

    Simulate image acquisition across Widefield, Confocal, STED, AiryScan, SMLM, and custom modalities with physically-grounded PSFs and noise models.

-   :material-chart-scatter-plot:{ .lg .middle } **Parameter Sweeps**

    ---

    Systematically test parameter combinations and compare outputs — quantify how labelling efficiency, probe distance, or exposure time affect your images.

</div>

---

## Quick Start

=== "pip"

    ```bash
    pip install vlab4mic
    ```

=== "pip (with Jupyter)"

    ```bash
    pip install vlab4mic vlab4micjupyter
    ```

Then run your first imaging simulation:

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(run_simulation=True)
```

---

## Available Notebooks

| Category | Description | Notebook | Colab |
| --- | --- | --- | --- |
| **Main Interface** | Create virtual samples and simulate image acquisition with multiple imaging modalities | [![Jupyter](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_main.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/VLab4Mic/blob/main/notebooks/vLab4Mic_main.ipynb) |
| **Parameter Sweeps** | Generate and analyze simulations over parameter ranges for optimization | [![Jupyter](https://img.shields.io/badge/jupyter-blue.svg?style=flat&logo=jupyter&logoColor=white)](https://github.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://githubtocolab.com/HenriquesLab/VLab4Mic/blob/main/notebooks/VLab4Mic_parameter_sweeps.ipynb) |

---

## Usage Options

| Option | Skills Needed | Best For |
|--------|---------------|----------|
| **Google Colab** | None | Beginners, teaching, quick experiments |
| **Local Jupyter Notebooks** | Basic installation | Smooth widget-based use |
| **Python Scripts** | Comfortable with code | Full flexibility & automation |
