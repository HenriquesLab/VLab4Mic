# Python Usage

When using VLab4Mic as a Python library, you get full programmatic control over every step of the simulation.

First make sure you have [installed VLab4Mic](getting-started.md) in your virtual environment.

---

## Option 1 — Run Scripts from the Command Line

Write your simulation in a Python script and run it from the command line:

```bash
python vlab_script.py
```

Here is a minimal script that runs an imaging simulation:

```python title="vlab_script.py"
from vlab4mic.experiments import image_vsample

images, noiseless, experiment = image_vsample(run_simulation=True)
```

For more complete examples, see the [Examples](examples.md) page.

---

## Option 2 — Interactive Python Interpreter

Start an interactive Python session:

```bash
python
```

Then use VLab4Mic interactively:

```python
>>> from vlab4mic.experiments import image_vsample
>>> images, noiseless, experiment = image_vsample(run_simulation=True)
```

---

## Getting Help

Use Python's built-in `help()` to explore the available parameters for any function:

```python
from vlab4mic.experiments import image_vsample
help(image_vsample)
```

```python
from vlab4mic.sweep_generator import run_parameter_sweep
help(run_parameter_sweep)
```

---

## Next Steps

<div class="grid cards" markdown>

-   :material-function:{ .lg .middle } **Methods Reference**

    ---

    Detailed usage guide for `image_vsample` and `run_parameter_sweep`.

    [Methods →](methods.md)

-   :material-table:{ .lg .middle } **Parameter Reference**

    ---

    Complete list of all parameters available for simulation and sweeps.

    [Parameter reference →](parameters.md)

-   :material-code-tags:{ .lg .middle } **Example Scripts**

    ---

    Browse ready-to-run example scripts covering common use cases.

    [Examples →](examples.md)

</div>
