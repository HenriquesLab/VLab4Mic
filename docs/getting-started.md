# Getting Started

VLab4Mic is compatible with **Python 3.10, 3.11, and 3.12** on macOS, Windows, and Linux.

---

## Step 1 — Create and Activate a Virtual Environment

!!! tip
    We strongly recommend creating a dedicated virtual environment to avoid dependency conflicts.

=== "Conda"

    Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you haven't already, then run:

    ```bash
    conda create --name myenv python=3.11
    conda activate myenv
    ```

=== "venv (built-in)"

    ```bash
    python -m venv myenv
    ```

    Activate the environment:

    === "macOS / Linux"

        ```bash
        source myenv/bin/activate
        ```

    === "Windows"

        ```bash
        myenv\Scripts\activate
        ```

Make sure your environment is active before continuing.

---

## Step 2 — Install VLab4Mic

**For Python scripts only:**

```bash
pip install vlab4mic
```

**Including Jupyter notebook support:**

```bash
pip install vlab4mic vlab4micjupyter
```

---

## Verify Installation

Check that VLab4Mic is installed correctly:

```python
import vlab4mic
print(vlab4mic.__version__)
```

Or run a minimal simulation:

```python
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(run_simulation=True)
print("Installation successful!")
```

---

## Next Steps

<div class="grid cards" markdown>

-   :material-notebook:{ .lg .middle } **Use Jupyter Notebooks**

    ---

    Run VLab4Mic without writing code in Google Colab or a local Jupyter Lab.

    [Notebooks guide →](notebooks.md)

-   :material-code-braces:{ .lg .middle } **Use Python Scripts**

    ---

    Run VLab4Mic from the command line or a Python interpreter.

    [Python usage guide →](python-usage.md)

</div>
