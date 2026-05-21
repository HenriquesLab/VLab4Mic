from vlab4mic.utils.transform.points_transforms import rotate_point
import numpy as np
from vlab4mic import experiments
import matplotlib.pyplot as plt
import os

def create_spaced_unitary_vectors(
        unitary_point = np.array([0, 0, 1]),
        axis_of_gyration = np.array([1, 0, 0]),
        number_of_vectors=1,
        start_angle_rad = 0,
        stop_angle_rad = 2*np.pi
        ):
    angles =np.linspace(start=start_angle_rad, stop=stop_angle_rad, num=number_of_vectors)
    vectors_list = []
    for i in angles:
        new = rotate_point(point=unitary_point, axis=axis_of_gyration, angle=i)
        vectors_list.append(new)
    return vectors_list

def create_xy_relative_positions(h_positions, v_positions, margin=0.1, spacing = 0.8, **kwargs):
    parameter_grid2 = np.mgrid[0:v_positions-1, 0:h_positions-1]
    parameter_grid_normalised2 = np.array(parameter_grid2)
    rotations_nomralised =  parameter_grid_normalised2[0]/(v_positions-1)
    axialpositions_normalised =  parameter_grid_normalised2[1]/(h_positions-1)
    #
    parameter_grid_normalised3 = [rotations_nomralised, axialpositions_normalised]
    #
    x_relative_positions2 = np.array(parameter_grid_normalised3[0]).flatten()
    y_relative_positions2 = np.array(parameter_grid_normalised3[1]).flatten()
    relative_positions2 = [np.array([margin + x*spacing, margin + y*spacing,0]) for x, y in zip(x_relative_positions2, y_relative_positions2)]
    # 
    return relative_positions2
    
    
def create_parameter_grid_as_list(h_positions, v_positions, margin=0.1, min_z = 0, max_z=100, spacing=1, start_angle_rad=0, stop_angle_rad=2*np.pi):
    # 
    vector_index = np.arange(stop=v_positions-1)
    zpositions = np.linspace(min_z, max_z, h_positions-1)
    # create 2D positions
    relative_positions_xyz = create_xy_relative_positions(h_positions, v_positions, margin, spacing)
    # create list of vectors
    vectors_list = create_spaced_unitary_vectors(number_of_vectors = v_positions, start_angle_rad=start_angle_rad, stop_angle_rad=stop_angle_rad)
    # Prepare grid
    X, Y = np.meshgrid(zpositions, vector_index)
    orientation_per_particle = []
    for i, rotation_idx_zposition in enumerate(zip(X.flatten(), Y.flatten())):
        orientation_per_particle.append(vectors_list[rotation_idx_zposition[1]])
        relative_positions_xyz[i][2] = (rotation_idx_zposition[0]/max_z)
    return relative_positions_xyz, orientation_per_particle

v_positions = 6
h_positions = 6
sample_dimensions=[1000,1000,100]
#max_z = 100
relative_positions_xyz, orientation_per_particle =create_parameter_grid_as_list(
    h_positions=h_positions,
    v_positions=v_positions,
    margin=0.1,
    min_z=0,
    max_z=sample_dimensions[2],
    spacing=1,
    stop_angle_rad=np.pi
    )

random_seed = 24
vsample, myexperiment2 = experiments.generate_virtual_sample(
    structure = "3J3Y",
    probe_template="Antibody",
    probe_DoL=4,
    probe_target_type="Sequence", 
    probe_target_value="SPRTLNA",
    particle_positions=relative_positions_xyz,
    orientation_per_particle=orientation_per_particle,
    sample_dimensions=sample_dimensions,
    clear_experiment=True,
    random_seed=random_seed
)
myexperiment2.add_modality(modality_name="SMLM")
myexperiment2.build(modules=["imager",])
myexperiment2.update_modality(modality_name="SMLM", depth_of_field_nm=100, lateral_precision=2, axial_precision=2, nlocalisations=30)
myexperiment2.set_modality_acq(modality_name="SMLM", exp_time=0.0001)
images2, noiseless2 = myexperiment2.run_simulation()


fig = plt.figure(figsize=(25, 15))
ax = fig.add_subplot(121, projection="3d")
myexperiment2.coordinate_field.show_field(view_init=[90,0,0], axis_object=ax, emitters_plotsize=1)
ax.set_xlabel("Orientation")
ax.set_ylabel("Zposition")
ax.set_zlabel("")
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_zticklabels([])
ax.xaxis


ax = fig.add_subplot(122)
ax.imshow(images2["SMLM"]["ch0"][0], cmap="Greys_r")
ax.set_title("SMLM")
ax.set_xticklabels([])
ax.set_yticklabels([])

filename = myexperiment2.date_as_string + 'vlab4mic_fig3_panelA.png'
filename2 = os.path.join(myexperiment2.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

#######################################################

fig2 = plt.figure(figsize=[20,30])
ax1 = fig2.add_subplot(1, 4, 1, projection='3d')
myexperiment2.remove_probes()
myexperiment2.add_probe(
    probe_template="Linker",
    probe_distance_to_epitope=0,
    probe_target_type="Sequence", 
    probe_target_value="SPRTLNA",
)
myexperiment2.build(modules=["particle",])
myexperiment2.structure.show_target_labels(
    axis_object=ax1,
    assembly_fraction=0.01,
    atoms_alpha=0.02,
    view_init=[90,0,0], 
    with_assembly_atoms=True,
    target_size=0, atoms_size=5, show_axis=False, reference_point=False)

ax2 = fig2.add_subplot(1, 4, 2, projection='3d')
myexperiment2.structure.show_target_labels(
    axis_object=ax2,
    assembly_fraction=0.01,
    atoms_alpha=0.01,
    view_init=[90,0,0], 
    with_assembly_atoms=False,
    target_size=5, atoms_size=5, show_axis=False, reference_point=False)

ax3= fig2.add_subplot(1,4, 3, projection='3d')
myexperiment2.remove_probes()
myexperiment2.add_probe(
    probe_template="anti-p24_primary_antibody_HIV",
    probe_DoL=6,
)
myexperiment2.build(modules=["particle",])
myexperiment2.particle.gen_axis_plot(axis_object=ax3, with_sources=True, source_plotsize=5, source_plotmarker="o", view_init=[90,0,0],
                                     xlim=[0,1000], ylim=[0,1000], zlim=[0,600], axesoff=True, emitter_plotsize=10)

ax4 = fig2.add_subplot(1,4, 4, projection='3d')
myexperiment2.remove_probes()
myexperiment2.add_probe()
myexperiment2.build(modules=["particle",])

myexperiment2.particle.gen_axis_plot(axis_object=ax4, with_sources=False, source_plotsize=0, source_plotmarker="o", view_init=[90,0,0],
                                     xlim=[0,1000], ylim=[0,1000], zlim=[0,600], emitter_plotsize=1)

filename = myexperiment2.date_as_string + 'vlab4mic_fig3_panelB.png'
filename2 = os.path.join(myexperiment2.output_directory, filename)
fig2.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

#############################
myexperiment2.remove_probes()
myexperiment2.add_probe(
    probe_template="HIV_capsid_p24_direct",
    labelling_efficiency=0.1
)
myexperiment2.clear_virtual_sample()
myexperiment2.set_virtualsample_params(
    sample_dimensions=[500,500,10],
    particle_positions=[[0.5,0.5,0]],
)
myexperiment2.update_modality(modality_name="SMLM", lateral_precision=2, lateral_resolution_nm=5)
myexperiment2.add_modality(modality_name="AiryScan")
myexperiment2.add_modality(modality_name="STED")
myexperiment2.build(modules=["particle","coordinate_field", "imager"])

images, noiseless = myexperiment2.run_simulation()


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(131)
ax.imshow(images["STED"]["ch0"][0], cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(132)
ax.imshow(images["SMLM"]["ch0"][0], cmap="grey")
ax.set_axis_off()
ax = fig.add_subplot(133)
ax.imshow(images["AiryScan"]["ch0"][0], cmap="grey")
ax.set_axis_off()

filename = myexperiment2.date_as_string + 'vlab4mic_fig3_panelC_mods.png'
filename2 = os.path.join(myexperiment2.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()