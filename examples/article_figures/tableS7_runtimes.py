"""Supp. Table S7: run times.

Times, on the current machine:
- one virtual sample (nuclear pore, 7R5K, GFP-nanobody, 10 particles in
  1 x 1 um) imaged in five modalities;
- N realisations of one particle for a distinguishability score (SMLM);
- a parameter sweep of 3 labelling efficiencies x 3 repetitions (STED).
The hardware is recorded with the results.
"""

import os
import platform
import time

import pandas as pd

from vlab4mic import experiments, sweep_generator

N_REALISATIONS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 100))
RANDOM_SEED = 7
MODALITIES = ["Widefield", "Confocal", "AiryScan", "STED", "SMLM"]

rows = []

start = time.perf_counter()
_, _, experiment = experiments.image_vsample(
    structure="7R5K", probe_template="GFP_w_nanobody", probe_target_type="Sequence",
    probe_target_value="ELAVGSL", number_of_particles=10, multimodal=MODALITIES,
    clear_experiment=True, random_seed=RANDOM_SEED,
)
rows.append(dict(task="One virtual sample, five imaging methods",
                 settings="7R5K, GFP-nanobody, 10 particles, 1 x 1 um",
                 seconds=time.perf_counter() - start))

start = time.perf_counter()
experiments.run_replicates(
    N_REALISATIONS, structure="7R5K", probe_template="NPC_Nup96_Cterminal_direct",
    number_of_particles=1, multimodal=["SMLM"], random_seed=RANDOM_SEED,
)
rows.append(dict(task=f"{N_REALISATIONS} realisations for one distinguishability score (one class)",
                 settings="7R5K, one particle, SMLM",
                 seconds=time.perf_counter() - start))

start = time.perf_counter()
sweep_generator.run_parameter_sweep(
    structures=["7R5K"], probe_templates=["NPC_Nup96_Cterminal_direct"],
    modalities=["STED"], labelling_efficiency=[0.3, 0.6, 1.0], sweep_repetitions=3,
    save_sweep_images=False, save_analysis_results=False, analysis_plots=False,
    run_analysis=True, random_seed=RANDOM_SEED,
)
rows.append(dict(task="Parameter sweep", settings="3 conditions x 3 repeats, STED",
                 seconds=time.perf_counter() - start))

results = pd.DataFrame(rows)
results["hardware"] = f"{platform.processor() or platform.machine()}, {platform.platform()}"
results["python"] = platform.python_version()
out = experiment.output_directory
results.to_csv(os.path.join(out, experiment.date_as_string + "tableS7_runtimes.csv"), index=False)
print(results)
