from vlab4mic import experiments
import matplotlib.pyplot as plt
import os
random_seed=24

probe1 = {
    "probe_template" : "Nanobody",
    "probe_name" : "Nanobody",
    "probe_DoL": 2 
}
probe2 = {
    "probe_template":"Antibody",
    "probe_name" : "Antibody",
    "probe_DoL": 3
}
probe3 = {
    "probe_template" : "GFP",
    "probe_name" : "GFP",
    "probe_DoL": None
}
probe4 = {
    "probe_template" : "GFP_w_nanobody",
    "probe_name" : "GFP_w_nanobody",
    "probe_DoL": 2 
}

vsample, my_experiment = experiments.generate_virtual_sample(
    probe_list=[probe1, probe2, probe3, probe4],
    clear_experiment=True,
    random_seed=random_seed
)
target_colour="#4163c8"
target_marker= "o"
plotting_parameters = dict(
    Nanobody = {
        "emitter_plotsize": 200,
        "emitter_plotcolour": target_colour,
        "emitter_plotmarker": target_marker,
        "atoms_plotsize":50,
        "atoms_fraction":1,
        "atoms_plotalpha" : 0.1,
        "central_axis_length": 100,
        "central_axis": False,
        "view_init":[0,0,0],
        "use_dol":False,
        "xlims": [-20, 20],
        "ylims":[-20, 20],
        "zlims":[-20, 20],
    },
    Antibody = {
        "emitter_plotsize": 100,
        "emitter_plotcolour": target_colour,
        "emitter_plotmarker": target_marker,
        "atoms_plotsize" : 50,
        "atoms_fraction":1,
        "atoms_plotalpha" : 0.01,
        "central_axis_length" :100,
        "central_axis" : False,
        "view_init":[-20,40,0],
        "use_dol":False,
        "xlims": [-100, 100],
        "ylims":[-100, 100],
        "zlims":[-50, 50],
    },
    GFP = {
        "emitter_plotsize" : 200,
        "emitter_plotcolour": target_colour,
        "emitter_plotmarker": target_marker,
        "atoms_plotsize" : 50,
        "atoms_fraction":1,
        "atoms_plotalpha" : 0.1,
        "central_axis_length" : 100,
        "central_axis" : False,
        "view_init":[20,60,0],
        "use_dol":False,
         "xlims": [-20, 20],
        "ylims":[-20, 20],
        "zlims":[-20, 20],
    },
    GFP_w_nanobody = {
        "emitter_plotsize" : 200,
        "emitter_plotcolour": target_colour,
        "emitter_plotmarker": target_marker,
        "atoms_plotsize" : 50,
        "atoms_fraction":1,
        "atoms_plotalpha" : 0.03,
        "central_axis_length" : 100,
        "central_axis" : False,
        "view_init":[20,90,30],
        "use_dol":False,
         "xlims": [-20, 20],
        "ylims":[-20, 20],
        "zlims":[-20, 20],
    }
)

fig = plt.figure(figsize=[30,20])
nprobes = len(list(plotting_parameters.keys()))
for i, probename in enumerate(plotting_parameters):
    ax = fig.add_subplot(1, nprobes, i+1, projection="3d")
    print(i, probename)
    my_experiment.particle.show_probe(
        axis_object=ax,
        probe_name=probename,
        with_structural_atoms=True,
        **plotting_parameters[probename],
        axesoff=True,)
plt.tight_layout()

filename = my_experiment.date_as_string + 'figS_probe_examples.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')