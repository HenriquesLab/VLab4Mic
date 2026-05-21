# Parameter Reference

Complete reference for all parameters accepted by `image_vsample` and `run_parameter_sweep`.

---

## Structure Parameters

| Parameter | Description | Notes |
| --- | --- | --- |
| **structure** | 4-letter PDB ID code, or path to a local structure file | Any PDB/CIF format file |
| **structure_is_path** | Treat `structure` as a local file path rather than a PDB ID (bool, default `False`) | — |
| **structure_axis_euler** | Euler angles `[α, β, γ]` (degrees) to rotate the structure's principal axis before positioning | Default `[0, 0, 0]` |
| **structure_global_normal_orientation** | Assign a fixed global normal direction for all structure epitopes (list or `None`) | — |
| **structure_format** | File format of the structure (`"CIF"` or `"PDB"`). Inferred from file extension when `structure_is_path` is `True` | Default `"CIF"` |

---

## Probe Parameters

| Parameter | Description |
| --- | --- |
| **probe_template** | Name of a pre-set probe configuration file (filename) |
| **probe_name** | Optional. Alternative name to identify this probe |
| **labelling_efficiency** | Probability that the probe binds a given epitope (float, default `1.0`) |
| **probe_target_type** | Type of target: `"Sequence"`, `"Residue"`, or `"Primary"` |
| **probe_target_value** | Specific value for the target type (e.g., amino acid sequence, residue info) |
| **probe_target_option** | Additional option for the probe target; used for secondary epitopes |
| **probe_distance_to_epitope** | Minimum distance from probe to epitope (float) |
| **probe_model** | 4-letter PDB ID for the atomic model of the probe itself |
| **probe_fluorophore** | Fluorophore name for emitters (e.g., `"AF647"`) |
| **probe_paratope** | If `probe_model` is set, the residue anchor point of the probe |
| **probe_conjugation_target_info** | If `probe_model` is set, dictionary specifying emitter conjugation sites |
| **probe_DoL** | Degree of labelling — efficiency of fluorophore conjugation (float) |
| **probe_seconday_epitope** | Amino acid sequence defining the epitope on a primary antibody model |
| **probe_wobble_theta** | Enable probe angular wobble (float or `None`) |
| **probe_steric_hindrance** | Minimum distance between adjacent epitopes (float) |
| **peptide_motif** | Dictionary specifying motif extraction for probe target sequence |
| **as_primary** | Treat the probe as a primary linker (bool) |

---

## Structural Integrity Parameters

| Parameter | Description |
| --- | --- |
| **structural_integrity** | Fraction of the particle rendered inaccessible to probes (float, 0–1) |
| **structural_integrity_small_cluster** | Maximum distance between epitopes for the first (small) inaccessible cluster grouping (Å) |
| **structural_integrity_large_cluster** | Minimum distance between epitopes for the second (large) inaccessible cluster grouping (Å) |

---

## Virtual Sample Parameters

| Parameter | Description |
| --- | --- |
| **virtual_sample_template** | Name of the configuration file for a virtual sample template |
| **sample_dimensions** | X, Y, Z field dimensions in nanometers (list of float) |
| **sample_inital_orientation** | Initial orientation applied to all labelled structures in the virtual sample (list or `None`) |
| **number_of_particles** | Number of independent particle copies distributed in the field (int) |
| **particle_positions** | Explicit list of relative particle positions (list of arrays) |
| **random_orientations** | Randomly assign orientations to each particle (bool, default `False`) |
| **orientation_per_particle** | Explicit list of orientations, one per particle |
| **xy_orientations** | List of orientation directions on the XY plane to sample from |
| **xz_orientations** | List of orientation directions on the XZ plane to sample from |
| **yz_orientations** | List of orientation directions on the YZ plane to sample from |
| **axial_offset** | List of Z offsets from the focus plane to sample from |
| **random_placing** | Place particles at random positions (bool, default `False`) |
| **random_rotations** | Apply random 3-D rotations to particles (bool, default `False`) |
| **rotation_per_particle** | Explicit list of rotations, one per particle |
| **rotation_angles** | List of rotation angles to sample from |
| **expansion_factor** | Scale factor applied to the structure coordinates (float, default `1`) |
| **minimal_distance** | Minimum allowed distance between particles (float) |
| **clear_probes** | Clear any probes already added to the experiment before adding new ones (bool, default `False`) |
| **clear_experiment** | Fully reset the experiment object after the simulation (bool, default `False`) |

---

## Modalities Parameters

| Parameter | Description |
| --- | --- |
| **pixelsize_nm** | Image pixel size in nanometers (float or int) |
| **lateral_resolution_nm** | Lateral (XY) resolution in nanometers (float or int) |
| **axial_resolution_nm** | Axial (Z) resolution in nanometers (float or int) |
| **psf_voxel_nm** | Voxel size used to render the PSF in nanometers (float or int, default `10`) |
| **depth_of_field_nm** | Depth of field in nanometers (float or int) |

---

## Multi-Probe Parameters

| Parameter | Description |
| --- | --- |
| **primary_probe** | Probe template name for the primary label in indirect labelling |
| **secondary_probe** | Probe template name for the secondary label in indirect labelling |
| **probe_list** | List of probe configuration dictionaries for adding multiple probes in a single call |
| **as_primary** | Treat this probe as a primary linker target for a secondary probe (bool) |

---

## Reproducibility

| Parameter | Description |
| --- | --- |
| **random_seed** | Integer seed passed to `numpy.random.seed` before any stochastic step (int or `None`) |

---

## Acquisition Parameters

| Parameter | Description |
| --- | --- |
| **exp_time** | Exposure time in seconds (float) |
| **noise** | Include photon and detector noise in the simulation (bool) |
| **nframes** | Number of frames to simulate (int) |
