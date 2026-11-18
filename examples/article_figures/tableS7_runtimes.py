"""Supp. Table S7: run times of the Fig. 2 analysis.

Times the clathrin distinguishability analysis of Fig. 2 with the same
models, probe, geometry, sweep machinery and scoring
(fig2_clathrin_distinguishability.py), on the current machine:

1. one realisation of one lattice imaged in five imaging methods;
2. one Fig. 2 condition (reference: 100% labelling, structural integrity 1):
   N realisations of each lattice in five imaging methods;
3. scoring that condition (classifier, bootstrap and permutation test);
4. the whole Fig. 2 run (every sweep block for both lattices), simulation
   and scoring timed separately. Set VLAB4MIC_TIME_FULL_RUN=0 to skip it.

The hardware, software versions, number of worker processes, number of
simulated images and peak memory are recorded with each time.

Outputs, in the VLab4Mic output directory: <date>tableS7_runtimes.csv
(one row per task; source data for Supp. Table S7).

Run:   python tableS7_runtimes.py
Environment variables: as fig2_clathrin_distinguishability.py
(VLAB4MIC_N_REALISATIONS, VLAB4MIC_N_PERMUTATIONS, VLAB4MIC_N_WORKERS), and
VLAB4MIC_TIME_FULL_RUN (1).
"""

import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from importlib.metadata import version

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig2_clathrin_distinguishability as fig2  # noqa: E402
from vlab4mic.analysis.distinguishability import (  # noqa: E402
    distinguishability_from_features,
    feature_metrics,
)
from vlab4mic.experiments import ExperimentParametrisation, output_path  # noqa: E402
from vlab4mic.sweep_generator import run_parameter_sweep  # noqa: E402

TIME_FULL_RUN = os.environ.get("VLAB4MIC_TIME_FULL_RUN", "1") not in ("0", "", "false", "False")
RANDOM_SEED = fig2.RANDOM_SEED
REFERENCE = dict(labelling_efficiency=[fig2.LABELLING_EFFICIENCY],
                 structural_integrity=[fig2.STRUCTURAL_INTEGRITY], **fig2.SI_CLUSTERS)


def hardware():
    """CPU model, logical cores, RAM (GB) and operating system."""
    cpu = platform.processor() or platform.machine()
    try:
        if sys.platform == "darwin":
            cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                 capture_output=True, text=True, check=True).stdout.strip()
        elif sys.platform == "win32":
            import winreg

            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                                 r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            cpu = winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
        elif os.path.exists("/proc/cpuinfo"):
            with open("/proc/cpuinfo") as handle:
                cpu = next(line.split(":", 1)[1].strip() for line in handle
                           if line.startswith("model name"))
    except Exception:
        pass
    try:
        if sys.platform == "win32":
            import ctypes

            class MemoryStatus(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

            status = MemoryStatus()
            status.dwLength = ctypes.sizeof(MemoryStatus)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
            ram_gb = status.ullTotalPhys / 1e9
        else:
            ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
    except Exception:
        ram_gb = float("nan")
    return dict(cpu=cpu, logical_cores=os.cpu_count(), ram_gb=round(ram_gb, 1),
                os=platform.platform())


def peak_memory_gb():
    """
    Peak resident memory (GB): of this process and its finished worker
    processes on macOS and Linux; of this process only on Windows (needs
    psutil); NaN if it cannot be measured.
    """
    try:
        import resource  # not available on Windows

        scale = 1e9 if sys.platform == "darwin" else 1e6  # ru_maxrss: bytes on macOS, kB on Linux
        own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / scale
        workers = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / scale
        return round(max(own, workers), 2)
    except ImportError:
        pass
    try:
        import psutil

        return round(psutil.Process().memory_info().peak_wset / 1e9, 2)
    except Exception:
        return float("nan")


def simulate_condition(path, n_realisations, modalities, seed):
    """Feature table of one lattice at the reference condition, as in Fig. 2."""
    sweep = run_parameter_sweep(
        structures=[path], structure_is_path=True, modalities=modalities,
        custom_metrics=feature_metrics(fig2.FEATURES, fov_nm=fig2.FOV_NM),
        random_seed=seed, save_sweep_images=False, output_directory=None,
        **{**fig2.COMMON, "sweep_repetitions": n_realisations}, **fig2.VENTRAL, **REFERENCE,
    )
    return sweep.analysis["dataframes"]


def score(features, by):
    return distinguishability_from_features(
        features, "model", ("dome", "flat"), fig2.FEATURES, by=by, controls=fig2.CONTROLS,
        replicate_column="Replica", n_permutations=fig2.N_PERMUTATIONS,
        random_state=RANDOM_SEED, n_workers=fig2.N_WORKERS, progress=False,
    )


def main():
    rows = []

    def record(task, settings, seconds, n_images, n_conditions, n_realisations):
        rows.append(dict(task=task, settings=settings, n_conditions=n_conditions,
                         realisations_per_lattice=n_realisations, n_images=n_images,
                         seconds=round(seconds, 1), minutes=round(seconds / 60, 2),
                         seconds_per_image=round(seconds / n_images, 3) if n_images else np.nan,
                         peak_memory_gb=peak_memory_gb()))
        print(f"{task}: {seconds:.1f} s", flush=True)

    paths = {name: fig2.download(url) for name, url in fig2.MODELS.items()}  # not timed
    seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(4)
    n_mod = len(fig2.MODALITIES)

    # 1. one realisation, one lattice, five imaging methods
    start = time.perf_counter()
    simulate_condition(paths["dome"], 1, fig2.MODALITIES, int(seeds[0]))
    record("One realisation, five imaging methods",
           "dome lattice, one particle, ventral geometry, 100% labelling",
           time.perf_counter() - start, n_mod, 1, 1)

    # 2. one Fig. 2 condition: N realisations of each lattice, five methods
    start = time.perf_counter()
    tables = []
    for k, (model, path) in enumerate(paths.items()):
        table = simulate_condition(path, fig2.N_REALISATIONS, fig2.MODALITIES, int(seeds[1 + k]))
        table["model"] = model
        tables.append(table)
    features = pd.concat(tables, ignore_index=True)
    n_images = 2 * fig2.N_REALISATIONS * n_mod
    record("One Fig. 2 condition, simulation",
           f"dome and flat lattices, {fig2.N_REALISATIONS} realisations each, five imaging methods",
           time.perf_counter() - start, n_images, 1, fig2.N_REALISATIONS)

    # 3. scoring that condition
    start = time.perf_counter()
    score(features, by=["modality_name"])
    record("One Fig. 2 condition, distinguishability score",
           f"five imaging methods; 5-fold cross-validation, 1,000 bootstrap resamples, "
           f"{fig2.N_PERMUTATIONS} permutations",
           time.perf_counter() - start, n_images, 1, fig2.N_REALISATIONS)

    # 4. the whole Fig. 2 run
    if TIME_FULL_RUN:
        blocks = [(model, path, block) for block in fig2.SWEEPS for model, path in paths.items()]
        block_seeds = np.random.SeedSequence(RANDOM_SEED).generate_state(len(blocks))
        jobs = [block + (int(seed), None) for block, seed in zip(blocks, block_seeds)]
        start = time.perf_counter()
        if fig2.N_WORKERS > 1:
            with ProcessPoolExecutor(fig2.N_WORKERS) as pool:
                tables = list(pool.map(fig2.simulate, jobs))
        else:
            tables = [fig2.simulate(job) for job in jobs]
        seconds = time.perf_counter() - start
        features = pd.concat(tables, ignore_index=True)
        n_images = len(features)
        # a condition is one parameter set, imaged in every method of its block
        n_conditions = len(features.drop_duplicates(
            [c for c in fig2.CONDITIONS if c != "modality_name"]))
        record("Whole Fig. 2 run, simulation",
               f"{len(fig2.SWEEPS)} sweep blocks x 2 lattices, {fig2.N_REALISATIONS} realisations "
               "per lattice and condition",
               seconds, n_images, n_conditions, fig2.N_REALISATIONS)
        start = time.perf_counter()
        score(features, by=fig2.CONDITIONS)
        record("Whole Fig. 2 run, distinguishability scores",
               f"every condition and imaging method; {fig2.N_PERMUTATIONS} permutations",
               time.perf_counter() - start, n_images, n_conditions, fig2.N_REALISATIONS)

    results = pd.DataFrame(rows)
    for key, value in hardware().items():
        results[key] = value
    results["workers"] = fig2.N_WORKERS
    results["python"] = platform.python_version()
    results["vlab4mic"] = version("vlab4mic")  # installed metadata; may lag an editable install
    try:
        results["git_commit"] = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__)),
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        results["git_commit"] = ""
    results["date"] = datetime.now().strftime("%Y-%m-%d")
    prefix = os.path.join(str(output_path), ExperimentParametrisation().date_as_string)
    path = prefix + "tableS7_runtimes.csv"
    results.to_csv(path, index=False)
    print(results.to_string(index=False))
    print(f"saved {path}")


if __name__ == "__main__":
    main()
