"""Supp. Table S1: linkage error and steric size of the probe templates.

Each probe template labels the Nup96 C-terminal motif (ELAVGSL) of the
nuclear pore (7R5K) at full labelling efficiency, with the wobble and
degree of labelling of the template (secondary antibodies: 15 degrees).
For every fluorophore, the linkage error is the distance to the nearest
target site; it is pooled over N_BUILDS labelled particles. The steric size
is the largest dimension of the probe model (estimate_probe_size), which is
the exclusion distance between bound probes.
"""

import contextlib
import io
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from vlab4mic import experiments
from vlab4mic.generate import labels
from vlab4mic.utils.io.yaml_functions import load_yaml
from vlab4mic.workflows import estimate_probe_size, probe_model

N_BUILDS = int(os.environ.get("VLAB4MIC_N_REALISATIONS", 20))
RANDOM_SEED = 3
TARGET = dict(probe_target_type="Sequence", probe_target_value="ELAVGSL")
PRIMARY = dict(probe_template="Antibody", probe_name="Primary", **TARGET)

PROBES = {
    "Primary antibody (Antibody, 1HZH)": [dict(probe_template="Antibody", **TARGET)],
    "Primary + secondary antibody": [
        dict(PRIMARY, as_primary=True),
        dict(probe_template="Antibody", probe_name="Secondary",
             probe_target_type="Primary", probe_target_value="Primary",
             probe_wobble_theta=15),
    ],
    "Primary antibody + secondary nanobody": [
        dict(PRIMARY, as_primary=True),
        dict(probe_template="Nanobody", probe_name="Secondary",
             probe_target_type="Primary", probe_target_value="Primary"),
    ],
    "Nanobody (7MFV)": [dict(probe_template="Nanobody", **TARGET)],
    "GFP-nanobody (6XZF)": [dict(probe_template="GFP_w_nanobody", **TARGET)],
    "GFP (6YLQ)": [dict(probe_template="GFP", **TARGET)],
    "SNAP-tag (6Y8P)": [dict(probe_template="SNAP-tag", **TARGET)],
    "Halo-tag (5UXZ)": [dict(probe_template="Halo-tag", **TARGET)],
    "mMaple (2HQK)": [dict(probe_template="mMaple", **TARGET)],
    "Linker": [dict(probe_template="Linker", **TARGET)],
}


def steric_size_nm(template, configuration_path):
    params = load_yaml(os.path.join(configuration_path, "probes", template + ".yaml"))
    label_object, label_params = labels.construct_label(label_config_dictionary=params)
    if not label_object.model or not label_object.model.get("ID"):
        return np.nan
    probe, emitters, *_ = probe_model(
        model=label_object.model, binding=label_object.binding,
        conjugation_sites=label_object.conjugation, epitope=label_object.epitope,
        config_dir=configuration_path,
    )
    return estimate_probe_size({"structural_atoms": probe.assembly_atoms, "coordinates": emitters}) / 10


with contextlib.redirect_stdout(io.StringIO()):
    _, experiment = experiments.generate_virtual_sample(
        structure="7R5K", probe_template="Antibody", clear_experiment=True,
        random_seed=RANDOM_SEED, **TARGET,
    )

rows = []
for name, probe_list in PROBES.items():
    distances = []
    with contextlib.redirect_stdout(io.StringIO()):
        experiment.remove_probes()
        for probe in probe_list:
            experiment.add_probe(labelling_efficiency=1.0, **probe)
        for _ in range(N_BUILDS):
            experiment.build(modules=["particle"])
            particle = experiment.particle
            targets = np.concatenate(
                [np.asarray(t["coordinates"]) for t in particle.source["targets"].values()]
            )
            for emitters in particle.emitters.values():
                emitters = np.asarray(emitters).reshape(-1, 3)
                if len(emitters):
                    distances.append(cKDTree(targets).query(emitters)[0])
    d = np.concatenate(distances) / 10 if distances else np.array([np.nan])  # Angstrom -> nm
    steric_template = probe_list[0]["probe_template"]
    rows.append(dict(
        probe=name,
        linkage_min_nm=round(float(np.nanmin(d)), 1),
        linkage_mean_nm=round(float(np.nanmean(d)), 1),
        linkage_median_nm=round(float(np.nanmedian(d)), 1),
        linkage_max_nm=round(float(np.nanmax(d)), 1),
        steric_size_nm=round(steric_size_nm(steric_template, experiment.configuration_path), 1),
        n_fluorophores=int(np.isfinite(d).sum()),
    ))

table = pd.DataFrame(rows)
path = os.path.join(experiment.output_directory, experiment.date_as_string + "tableS1_probe_geometry.csv")
table.to_csv(path, index=False)
experiment.save_parameters(experiment.output_directory, name="tableS1_probe_geometry")
print(table.to_string(index=False))
