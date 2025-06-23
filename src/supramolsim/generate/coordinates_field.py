import numpy as np
import copy
import matplotlib.pyplot as plt
import yaml
from ..analysis import metrics

from .labelled_instance import LabeledInstance
from .molecular_structure import MolecularReplicates
from .labels import Label

from ..utils.transform.points_transforms import transform_displace_set
from ..utils.sample import arrays as sampl
from ..utils.data_format.visualisation import format_coordinates
from ..utils.visualisation.matplotlib_plots import add_ax_scatter


class Field:
    def __init__(self):
        self.params = {}
        self.params["dimension_sizes"] = [1000.0, 1000.0, 10.0]
        self.params["scale"] = 1e-9
        self.params["relative_reference_point"] = [0.5, 0.5, 0.01]
        self.params["absolute_reference_point"] = None

        self.molecules = None

        self.molecules_params = {}
        self.molecules_params["nMolecules"] = 1
        self.molecules_params["relative_positions"] = np.array(
            [
                [0.5, 0.5, 0.01],
            ]
        ).reshape((1, 3))
        self.molecules_params["absolute_positions"] = None
        self.molecules_params["orientations"] = None
        self.molecules_params["minimal_distance"] = None

        self.emitters_per_fluorophore = {}  # this can be exported

        #
        self.ximage_dimensions = None  # fov_info["imagesize"][0]  # in pixels
        self.yimage_dimensions = None  # fov_info["imagesize"][1]  # in pixels
        self.zimage_dimensions = None  # fov_info["imagesize"][2]  # in pixels
        self.z_offset = 0  # in scale units
        self.pixelsize = None  # fov_info["pixelsize"]  # in scale units
        self.relative_pos = None  # fov_info["relative_positions"]
        self.orientations = []
        self.absolute_pos = []
        self.fov_relative_ref_pt_xy = np.array([0.5, 0.5, 0.01])
        self.fov_absolute_ref_pt_xyself = None
        self.scale = 1e-9  # in meters
        self.particles = []
        self.label_names = []
        self.colours = {}
        self.fluorophre_emitters = {}
        self.fluo2labels = []
        self.plotting_params = dict()
        self.fluoparams = dict()
        self.random_placing = False
        self.random_orientations = False
        # print(f'Working scale of the Field of View is {self.scale} meters')

    # methods to initialise field parameteres
    def _set_dimension_sizes(self, dims):
        # limit x y and z in units of scale
        self.params["dimension_sizes"] = dims

    def _set_relative_reference(self, rel_ref: list):
        self.params["relative_reference_point"] = np.array(rel_ref)

    def _set_absolute_reference(self, abs_ref: np.ndarray):
        self.params["absolute_reference_point"] = abs_ref

    def _set_scale(self, scale: float):
        self.params["scale"] = scale

    def set_params(self, **kwargs):
        """
        :param kwargs: same as self.params.keys()
        """
        for key, value in kwargs.items():
            self.params[key] = value

    def init_from_file(self, field_yaml):
        with open(field_yaml, "r") as f:
            field_params = yaml.safe_load(f)
        params = dict(field_params["params"])
        molecules = dict(field_params["molecules_params"])
        export = dict(field_params["export_params"])
        # print("Field parameters for placing and orientation of labelled particles:")
        for key, value in molecules.items():
            print(key, ": ", value)
        # set general parameters
        self.set_params(**params)
        # set molecules params
        self.set_molecules_params(**molecules)

    def create_minimal_field(self, nmolecules=1, random_placing=False, random_orientations=False, **kwargs):
        fluo_name = "AF647"
        self.calculate_absolute_reference()
        if "minimal_distance" in kwargs.keys():
            self.molecules_params["minimal_distance"] = kwargs["minimal_distance"]
        if "sample_dimensions" in kwargs.keys():
            self._set_dimension_sizes(kwargs["sample_dimensions"])
            self.calculate_absolute_reference()
        if "relative_positions" in kwargs.keys():
            print("Using relative positons to initialise field")
            nmolecules = len(kwargs["relative_positions"])
            self.set_molecule_param("nMolecules", nmolecules)
            self.molecules_params["relative_positions"] = kwargs["relative_positions"]
            self._gen_abs_from_rel_positions()
            self.fluorophre_emitters = {
                fluo_name: self.get_molecule_param("absolute_positions")
            }
            # ignore nmolecules and randomplacing
            # use thise values to set particles
        elif nmolecules > 1:
            print("More than 1 position needed. Randomising.")
            self.set_molecule_param("nMolecules", nmolecules)
            # if no list was passed, then we will use the default nmolecules parameter
            # having more than one particle necesarily use random placing
            self.random_placing = True
            self.generate_random_positions()
            self._gen_abs_from_rel_positions()
            self.fluorophre_emitters = {
                fluo_name: self.get_molecule_param("absolute_positions")
            }
        else:
            # if none of the above, this is the case of a single particle in the center
            #point = self.get_field_param("absolute_reference_point")
            self.set_molecule_param("nMolecules", 1)
            if random_placing:
                self.random_placing = True
                self.generate_random_positions()
            self._gen_abs_from_rel_positions()
            point = self.get_molecule_param("absolute_positions")
            self.fluorophre_emitters = {fluo_name: point.reshape(1, 3)}
        self._set_fluo_plotting_params(fluo_name)
        self.random_orientations = random_orientations



    def calculate_absolute_reference(self):
        relative_pt = self.get_field_param("relative_reference_point")
        dimensions = self.get_field_param("dimension_sizes")
        abs_posx = relative_pt[0] * dimensions[0]
        abs_posy = relative_pt[1] * dimensions[1]
        abs_posz = relative_pt[2] * dimensions[2]
        self._set_absolute_reference(np.array([abs_posx, abs_posy, abs_posz]))

    def get_field_param(self, parameter: str):
        return self.params[parameter]

    def change_number_of_molecules(self, n: int):
        self.set_molecule_param("nMolecules", n)
        self.generate_random_positions()

    # molecule parameters
    def get_molecule_param(self, parameter: str):
        return self.molecules_params[parameter]

    def set_molecule_param(self, parameter: str, value):
        self.molecules_params[parameter] = value

    def set_molecules_params(
        self,
        nMolecules: int,
        random_positions,
        random_orientations,
        random_rotations,
        **kwargs,
    ):
        # self.change_number_of_molecules(nMolecules)
        self.set_molecule_param("nMolecules", nMolecules)
        if random_positions:
            self.generate_random_positions()
        if random_orientations:
            self.generate_random_orientations()
        if random_rotations:
            pass
        for key, value in kwargs.items():
            self.molecules_params[key] = value

    def show_params(self):
        print("Field parameters")
        for key, val in self.params.items():
            print(key, ": ", val)
        print("Molecule objects parameters")
        for key, val in self.molecules_params.items():
            print(key, ": ", val)

    # methods to prime molecule positions, orientations...
    def generate_random_positions(self):
        """
        Generate positions randomly constrained by the dimensions of the
        field and the minimal molecule distance if exists

        Sets the molecule parameter "absolute positions"
        """
        npositions = self.get_molecule_param("nMolecules")
        print(f"nmolecules: {npositions}")
        # print(f"Total positions: {npositions}")
        # this can only work after the size of field has been established
        if self.molecules_params["minimal_distance"] is not None:
            print(f"distributing with minimal distance: {self.molecules_params['minimal_distance']}")
            self._random_pos_minimal_dist(npositions)
        else:
            print("Generating unconstrained random positions")
            xrel = np.random.uniform(size=npositions)
            yrel = np.random.uniform(size=npositions)
            zrel = np.random.uniform(size=npositions)
            rand = list(zip(xrel, yrel, zrel))
            self.molecules_params["relative_positions"] = rand
            self._gen_abs_from_rel_positions()

    def generate_random_orientations(self):
        # give new orientation
        norientations = self.get_molecule_param("nMolecules")
        orientations = []
        for i in range(norientations):
            orientations.append(np.array(sampl.sample_spherical_normalised(1, ndim=3)))
        self.set_molecule_param("orientations", orientations)

    def generate_global_orientation(self, global_orientation = None):
        # give new orientation
        if global_orientation is not None:
            norientations = self.get_molecule_param("nMolecules")
            orientations = []
            for i in range(norientations):
                orientations.append(global_orientation)
            self.set_molecule_param("orientations", orientations)

    def _set_molecule_minimal_distance(self, dist):
        if dist > 0:
            print(
                f"minimal distance from particles set to {dist}. "
                f"FOV scale unit is {self.scale} meters"
            )
            self.molecules_params["minimal_distance"] = dist

    def _random_pos_minimal_dist(self, n):
        # convert minimal distance in relative units
        selected_positions = []
        selected_relative = []
        minimal_distance = self.get_molecule_param("minimal_distance")
        dimension_sizes = self.get_field_param("dimension_sizes")
        maxabs_posx = dimension_sizes[0]
        maxabs_posy = dimension_sizes[1]
        maxabs_posz = dimension_sizes[2]
        # print(f"Max values: {maxabs_posx},{maxabs_posy}, {maxabs_posz}")
        # sample the first point
        x = np.random.uniform(0, maxabs_posx, size=1)
        y = np.random.uniform(0, maxabs_posy, size=1)
        z = np.random.uniform(0, maxabs_posz, size=1)
        v0 = np.array([x, y, z]).reshape(-1)
        selected_positions.append(v0)
        flag = 0
        while len(selected_positions) < n and flag < 10000:
            flag = flag + 1
            x = np.random.uniform(0, maxabs_posx, size=1)
            y = np.random.uniform(0, maxabs_posy, size=1)
            z = np.random.uniform(0, maxabs_posz, size=1)
            new = np.array([x, y, z]).reshape(-1)
            is_available = 1
            for pos in range(len(selected_positions)):
                # verify is the next point is within the minimum distance
                # print(np.linalg.norm(new - selected_positions[pos]))
                if np.linalg.norm(new - selected_positions[pos]) < minimal_distance:
                    is_available = 0
                    #break
            if is_available:
                selected_positions.append(new)
        for abs_pos in selected_positions:
            rel_pos = abs_pos
            rel_pos[0] = abs_pos[0]/maxabs_posx
            rel_pos[1] = abs_pos[1]/maxabs_posy
            rel_pos[2] = abs_pos[2]/maxabs_posz
            selected_relative.append(rel_pos)
        # print(f"end flag : {flag}")
        self.molecules_params["relative_positions"] = selected_relative
        #self.set_molecule_param("absolute_positions", selected_positions)
        #self.absolute_pos = selected_positions

    def _calculate_absolute_position(self, relative_pos):
        fieldsizes = self.get_field_param("dimension_sizes")
        abs_posx = relative_pos[0] * fieldsizes[0]
        abs_posy = relative_pos[1] * fieldsizes[1]
        abs_posz = relative_pos[2] * fieldsizes[2]
        return np.array([abs_posx, abs_posy, abs_posz])

    def _gen_abs_from_rel_positions(self):
        abs_pos = []
        for pos in self.molecules_params["relative_positions"]:
            abs_pos.append(self._calculate_absolute_position(pos))
        self.set_molecule_param("absolute_positions", np.array(abs_pos))

    # methods for fluorophores
    def _load_fluorophore_params(self, **fluodictionary):
        fname = fluodictionary["fluorophore_name"]
        fluo_plotting_params = dict(
            plotsize=20,
            plotalpha=1,
            plotmarker="o",
            plotcolour=fluodictionary["plotcolour"],
        )
        self.plotting_params[fname] = fluo_plotting_params
        self.fluoparams[fname] = dict(fluodictionary)
        # here add all other relevant parameters of a fluorophore in a separate atribute

    def _set_fluo_plotting_params(self, fluoname, plotcolour="m"):
        fluo_plotting_params = dict(
            plotsize=20, plotalpha=1, plotmarker="o", plotcolour=plotcolour
        )
        self.plotting_params[fluoname] = fluo_plotting_params

    def get_fluorophore_params(self):
        return self.fluoparams

    # working with macromolecules
    def create_molecules_from_InstanceObject(self, InstancePrototype: LabeledInstance):
        reps = self.get_molecule_param("nMolecules")
        particle_copy = copy.deepcopy(InstancePrototype)
        self.molecules_default_orientation = particle_copy._get_source_parameter("axis")
        particle_copy.scale_coordinates_system(self.get_field_param("scale"))
        if self.molecules_params["minimal_distance"] is None:
            self._set_molecule_minimal_distance(
                dist=particle_copy.radial_hindance
            )
        if self.random_placing:
            self.generate_random_positions()
        self._gen_abs_from_rel_positions()
        # due to constraints, the actual number of particles might not be reps
        nmolecules = len(self.molecules_params["absolute_positions"])
        if reps > nmolecules:
            reps = nmolecules
            print(f"Updating Number of molecules: {reps}")
            self.set_molecule_param("nMolecules", reps)
        # prepare plotting params
        self.labels_plotting_params = dict(particle_copy.plotting_params)
        self.fluo2labels = dict(particle_copy.fluo2labels)
        self._set_plotting_params()
        molecules = []
        for r in range(reps):
            molecules.append(copy.deepcopy(particle_copy))
        self.molecules = molecules
        if self.random_orientations:
            self.generate_random_orientations()
            # self.relabel_molecules()
        self.relabel_molecules()
        

    def _create_instances_from_pdb(self, cif_f, structure_id, label_file):
        # this is exactly the same as the high level function
        structure_dictionary = {
            "file": cif_f,
            "identifier": structure_id,
            "format": "CIF",
        }
        with open(label_file, "r") as f:
            label_params = yaml.safe_load(f)
        # create the molecular structure
        Molecularstructure = MolecularReplicates(structure_dictionary)
        nlabels = len(label_params["labels"])
        # create labels
        label_handler = [Label() for i in range(nlabels)]
        for la in range(nlabels):
            label_handler[la].set_params(**label_params["labels"][la])
        # add labels to molecular structure
        for labelobject in label_handler:
            Molecularstructure.add_label(labelobject)
        inst_builder = Molecularstructure.create_instance_builder()
        particle = LabeledInstance()
        particle.load_source(**inst_builder)  # define the source
        particle.load_label_file(label_file)  # load the labelling entities
        particle.generate_instance()
        return particle

    def show_molecules(self):
        if self.molecules is not None:
            for mol in self.molecules:
                mol.show_instance()
                mol.get_ref_point()

    def relabel_molecules(self):
        if self.molecules:
            for mol in self.molecules:
                mol.generate_instance()
                # mol.show_instance()
                mol.scale_coordinates_system(self.get_field_param("scale"))

    def relocate_molecules(self):
        if self.molecules:
            if self.get_molecule_param("absolute_positions") is not None:
                # move each molecule to their absolute position
                for mol, pos in zip(
                    self.molecules, self.get_molecule_param("absolute_positions")
                ):
                    mol.transform_translate(pos)

    def reorient_molecules(self):
        if self.get_molecule_param("orientations") is not None:
            for mol, ori in zip(
                self.molecules, self.get_molecule_param("orientations")
            ):
                mol.transform_reorient(ori)
        else:
            pass
            # print("molecule orientations has not been set. No reorientation done.")

    def _set_plotting_params(self):
        for fluoname, labelname in self.fluo2labels.items():
            print(fluoname, labelname[0])
            fluo_colour = self.labels_plotting_params[labelname[0]]["plotcolour"]
            self._set_fluo_plotting_params(fluoname, fluo_colour)

    # CONSTRUCTION OF STATIC FIELD
    # methods for creating the coordinates emitters per fluorophore species
    # the rationale is, there will be emitters that correspond to one or other emitters
    # and all of them will be able to be visualised if the optics define it

    def construct_static_field(self, relocate=True, reorient=True):
        if self.molecules:
            # verify absolute positions exist
            if self.get_molecule_param("absolute_positions") is None:
                self._gen_abs_from_rel_positions()
            if self.params["absolute_reference_point"] is None:
                self.calculate_absolute_reference()
            # the only thing that will be needed every time is the placing of the particles
            # in the field
            if relocate:
                self.relocate_molecules()
            if reorient:
                self.reorient_molecules()
            # pull emitters by fluorophores
            self._construct_channels_by_fluorophores()
        else:
            default_fluorophore = "AF647"
            #self.get_molecule_param("absolute_positions")
            self.fluorophre_emitters = {}
            self.fluorophre_emitters[default_fluorophore] = self.get_molecule_param("absolute_positions")

    def _construct_channels_by_fluorophores(self):
        """
        group all emitter positions that correspond to a fluorophore species
        """
        nmolecules = self.get_molecule_param("nMolecules")
        if nmolecules == 0:
            print("Field contains no macromolecules")
        else:
            chs = {}
            if self.fluo2labels is not None:
                for fluoname in list(self.fluo2labels.keys()):
                    chs[fluoname] = self._pull_coordinates_per_fluorophore(fluoname)
                self.fluorophre_emitters = chs
            else:
                print("No relation from fluorophore to labels")

    def _pull_coordinates_per_fluorophore(self, fluophorename):
        """
        fluorophorename is a string that identifies a specific species of fluorophore
        Iterates over the particles of InstanceLabelled and gets its emitters
        Returns all emitters corresponding to that emitter
        """
        # first, verify that there are particles (Labelled instance)
        nmolecules = self.get_molecule_param("nMolecules")
        if nmolecules == 0:
            print("Field contains no particles")
        else:
            pulled = np.empty((0, 3))
            for mol in self.molecules:
                pulled = np.vstack(
                    (pulled, mol.get_emitters_by_fluorophore(fluophorename))
                )  # each particle is an object of LabeledInstance
            return pulled

    def get_emitters_by_fluorophores(self):
        if len(self.fluorophre_emitters) == 0:
            print("Constructing channels by fluorophore type")
            self._construct_channels_by_fluorophores()
        return dict(self.fluorophre_emitters)

    # methods for visualisation
    def show_field(self, fluo_type="all", view_init=[30, 0, 0], initial_pos=True, return_fig=False):
        dimension_sizes = self.get_field_param("dimension_sizes")
        xx, yy = np.meshgrid(
            range(int(dimension_sizes[0])), range(int(dimension_sizes[1]))
        )
        zz = (yy) * 0
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(xx, yy, zz, alpha=0.2)

        if initial_pos:
            if self.molecules_params["absolute_positions"] is not None:
                add_ax_scatter(
                    ax,
                    format_coordinates(
                        np.array(self.molecules_params["absolute_positions"])
                    ),
                )
        if len(self.fluorophre_emitters) == 0:
            print("Constructing static field")
            self.construct_static_field()
        if fluo_type == "all":
            print("Showing all fluorophores")
            for fname in self.fluorophre_emitters.keys():
                print(fname)
                add_ax_scatter(
                    ax,
                    format_coordinates(
                        self.fluorophre_emitters[fname], **self.plotting_params[fname]
                    ),
                )
        else:
            add_ax_scatter(
                ax,
                format_coordinates(
                    self.fluorophre_emitters[fluo_type],
                    **self.plotting_params[fluo_type],
                ),
            )
        ax.set_box_aspect(
            [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
        )
        ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        if return_fig:
            return fig
        else:
            fig.show()

    def expand_isotropically(self, factor):
        if factor <= 1:
            print("Input factor has to be greater than 1")
        else:
            self.fov_absolute_ref_pt_xyself = self._calculate_absolute_position(
                self.fov_relative_ref_pt_xy
            )
            print(
                "Warining: this method changes the final emitter sizes"
                "without altering the scales. "
                "This is intended to be used as a last step"
            )
            for fname in self.fluorophre_emitters.keys():
                print(fname)
                print("expanding")
                # expanded[fname] = self.fluorophre_emitters[fname] * factor
                self.fluorophre_emitters[fname] = (
                    self.fluorophre_emitters[fname] * factor
                )
                print(self.fluorophre_emitters[fname])
            new_abs_ref = self.fov_absolute_ref_pt_xyself * factor
            print(f"newabs {new_abs_ref}")
            for fname in self.fluorophre_emitters.keys():
                print("adjusting")
                # expanded[fname] = self.fluorophre_emitters[fname] * factor
                self.fluorophre_emitters[fname], _ = transform_displace_set(
                    self.fluorophre_emitters[fname],
                    new_abs_ref,
                    self.fov_absolute_ref_pt_xyself,
                )
                print(self.fluorophre_emitters[fname])

    def export_field(self):
        if self.fluorophre_emitters is None:
            self.construct_static_field()
        field_emitters = self.fluorophre_emitters
        field_scale = self.get_field_param("scale")
        field_sizes = self.get_field_param("dimension_sizes")
        plotting_params = self.plotting_params
        reference_point = self.get_field_param("absolute_reference_point")
        export_field = dict(
            field_emitters=field_emitters,
            field_scale=field_scale,
            plotting_params=plotting_params,
            reference_point=reference_point,
            field_sizes=field_sizes,
        )
        return export_field


def create_min_field(number_of_particles=1, random_placing=False, random_orientations=False, prints=False, **kwargs):
    if prints:
        print("Initialising default field")
    coordinates_field = Field()
    #if molecule_pars:
    #    for key, value in molecule_pars.items():
    #        coordinates_field.molecules_params[key] = value
    coordinates_field.create_minimal_field(
        nmolecules=number_of_particles,
        random_placing=random_placing, 
        random_orientations=random_orientations,
        **kwargs
    )
    return coordinates_field


def gen_positions_from_image(img, mode="mask", pixelsize = None, **kwargs):
    npixels = list(img.shape)
    image_physical_size = np.zeros(shape=(2))
    image_physical_size[0] = pixelsize * npixels[0]
    image_physical_size[1] = pixelsize * npixels[1]
    if mode == "mask":
        if "npositions" not in kwargs.keys():
            npositions = 1
        else:
            npositions = kwargs["npositions"]
        if "min_distance" not in kwargs.keys():
            min_distance = 1
        else:
            min_distance = kwargs["min_distance"]
        #
        pixel_positions = sampl.get_random_pixels(
            img, 
            num_pixels=npositions, 
            min_distance=min_distance
        )
    elif mode == "localmaxima":
        if "background" not in kwargs.keys():
            background = None
        else:
            background = kwargs["background"]
        if "sigma" not in kwargs.keys():
            sigma = None
        else:
            sigma = kwargs["sigma"]
        if "threshold" not in kwargs.keys():
            threshold = None
        else:
            threshold = kwargs["threshold"]
        if "min_distance" not in kwargs.keys():
            min_distance = 1
        else:
            min_distance = kwargs["min_distance"]
        pixel_positions, img_processed = metrics.local_maxima_positions(
            img, 
            min_distance=min_distance, 
            threshold=threshold, 
            sigma=sigma, 
            background=background)
    xyz_relative = metrics.pixel_positions_to_relative(
            pixel_positions,
            image_sizes=image_physical_size,
            pixelsize=pixelsize
        )
    return xyz_relative, image_physical_size