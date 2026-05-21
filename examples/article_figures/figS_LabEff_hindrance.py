import matplotlib.pyplot as plt
from vlab4mic import experiments
import numpy as np
from IPython.utils import io


import os
random_seed=24

probe_parameters = dict(
    probe_template="Antibody",
    probe_target_type="Sequence",
    probe_target_value="EQATETQ",
    probe_DoL=4,
    probe_distance_to_epitope=10,
)
strcture = "1XI5"

images , noiseless_,  experiment = experiments.image_vsample(
    structure=strcture,
    **probe_parameters,
    clear_experiment=True,
    run_simulation=False,
    random_seed=random_seed
)

fig = plt.figure(figsize=[20,10])
total = 5
for n, efficiency in enumerate(np.linspace(start=0.2, stop=1, num=total)):
    ax = fig.add_subplot(1, total, n+1, projection="3d")
    with io.capture_output() as captured:
        experiment.remove_probes()
        experiment.add_probe(
            **probe_parameters,
            labelling_efficiency=efficiency)
        experiment.build(modules = ["particle",])
    experiment.particle.gen_axis_plot(axis_object=ax, emitter_plotsize=10, with_sources=True)
    title = strcture + "\n" + probe_parameters["probe_template"]
    title += " \n  Labelling Efficiency: " + str(efficiency)[0:4]
    ax.set_title(title)

filename = os.path.join(experiment.output_directory, 'figS_LabellingEfficiency.pdf')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=[20,10])
total = 5
emitter_plotsize = 10
labelling_efficiency = 1
for n, hindrance in enumerate(np.linspace(start=0, stop=200, num=total-1)):
    ax = fig.add_subplot(1, total, n+1, projection="3d")
    with io.capture_output() as captured:
        experiment.remove_probes()
        experiment.add_probe(
            **probe_parameters,
            labelling_efficiency=labelling_efficiency,
            probe_steric_hindrance=hindrance)
        experiment.build(modules = ["particle",])
    experiment.particle.gen_axis_plot(axis_object=ax, emitter_plotsize=emitter_plotsize, with_sources=True)
    title = strcture + "\n" + probe_parameters["probe_template"]
    title += " \n  Labelling Efficiency: " + str(labelling_efficiency)[0:4]
    title += " \n  Probe Steric Hindrance (A): " + str(experiment.particle.labels["Antibody"]["minimal_distance"])[0:4]
    ax.set_title(title)
### Estimate from size
ax = fig.add_subplot(1, total, total, projection="3d")
with io.capture_output() as captured:
    experiment.remove_probes()
    experiment.add_probe(
       **probe_parameters,
        labelling_efficiency=labelling_efficiency,
        probe_steric_hindrance="estimate")
    experiment.build(modules = ["particle",])
experiment.particle.gen_axis_plot(axis_object=ax, emitter_plotsize=emitter_plotsize, with_sources=True)
title = title = strcture + "\n" + probe_parameters["probe_template"]
title += " \n  Labelling Efficiency: " + str(labelling_efficiency)[0:4]
title += " \n  Estimated Steric Hindrance (A): " + str(experiment.particle.labels["Antibody"]["minimal_distance"])[0:6]
ax.set_title(title)

filename = os.path.join(experiment.output_directory, 'figS_ProbeHindrance.pdf')
fig.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()