import copy
import numpy as np
import warnings
import matplotlib.pyplot as plt
import yaml

from ..utils.transform.points_transforms import (
    rotate_pts_by_vector,
    transform_displace_set,
)
from ..utils.data_format.visualisation import format_coordinates, set_colorplot
from ..utils.visualisation.matplotlib_plots import add_ax_scatter, draw1nomral_segment
from ..utils.transform.cif_builder import create_instance_label
from ..utils.transform.defects import xmersubset_byclustering
from ..utils.data_format.structural_format import builder_format

class LabeledInstance:
    def __init__(self):
        self.params = {}
        self.params["ref_point"] = None
        self.params["scale"] = None
        self.source = {}
        self.source["targets"] = None
        self.source["reference_pt"] = None
        self.source["scale"] = None
        self.source["axis"] = None
        self.source["info"] = None
        self.labels = {}
        self.labelnames = list()
        # instance related
        self.emitters = []
        self.axis = dict()
        self.plotting_params = dict()
        self.radial_hindance = None
        self.defects = False
        self.defects_params = dict()
        # Contains the subset of epitopes after defect calculation
        self.defects_target_normals = None
        self.fluo2labels = []
        self.status = dict(source=False, labels=False)
        # attributes for sequencial labelling
        self.sequential_labelling = False
        self.primary = {}
        self.primary["targets"] = {}
        self.primary["reference_pt"] = None
        self.primary["scale"] = None
        self.primary["axis"] = None
        self.primary["info"] = None
        #
        self.secondary = {}

    def set_params(self, **kwargs):
        """
        Method used to set the parameters of LabeledInstance using keyword arguments.
        :param kwargs: same as self.params.keys()
        """
        for key, value in kwargs.items():
            self.params[key] = value

    # methods to define source entity
    def _set_source_targets(self, targets: dict):
        # keys are the label agent used to define such targets
        self.source["targets"] = targets

    def _set_source_reference(self, ref):
        self.source["reference_pt"] = ref

    def _set_source_scale(self, scale):
        self.source["scale"] = scale

    def _set_source_axis(self, axis: dict):
        self.source["axis"] = axis

    def _set_source_info(self, info: str):
        self.source["info"] = info

    def load_source(
        self,
        targets: dict,
        reference_point: np.ndarray,  # all points share the same reference point
        scale: float,
        axis=dict(pivot=np.array([0, 0, 0]), direction=np.array([0, 0, 1])),
        info="NA",
        **kwargs,
    ):
        """
        Method to set the ideal label locations and information
        for label targeting.

        targets: dictionary with target locations per label entity
        targets should also contain a list of surface normals for each target
        """
        self._set_source_targets(dict(targets))  # making an explicit copy
        self._set_source_reference(reference_point)
        self._set_source_scale(scale)
        self._set_source_axis(axis)
        self._set_source_info(info)
        self.status["source"] = True

    def _get_source_target_normals(self, target_name):
        if self.source["targets"][target_name]["normals"] is None:
            return None
        else:
            return self.source["targets"][target_name]["normals"]

    def _get_source_coords_normals(self, target_name):
        # output is a dictionary with keys "coordinates" and "normals"
        return self.source["targets"][target_name]

    def _get_source_target_names(self):
        return self.source["targets"].keys()

    def _get_source_target_label(self, target_name):
        if target_name in self.labels.keys():
            return self.labels[target_name]
        else:
            return None

    def _get_source_parameter(self, key):
        return self.source[key]

    def load_label_file(self, label_yaml):
        """
        Intended to be used with a predefined label structure
        that contains emitters expplicitly stated if exists
        and with those emitters the first two points
        correspond to the pivot and axis definition
        """
        with open(label_yaml, "r") as f:
            label_params = yaml.safe_load(f)
        targets = dict()
        for lab in label_params["labels"]:
            # print(lab)
            fluorophore = lab["fluorophore"]
            coordinates = None
            if "coordinates" in lab.keys():
                # print(lab["coordinates"])
                coordinates = np.array(lab["coordinates"])
                # print(coordinates, type(coordinates))
            targets[fluorophore] = coordinates
            self.load_label(
                targets=targets,
                label_name=lab["label_name"],
                labellig_efficiency=lab["labeling_efficiency"],
            )

    # methods to specify labels
    def load_label(
        self,
        targets: dict,
        scale: float = 1e-9,
        axis=dict(pivot=None, direction=None),
        label_name="NA",
        labellig_efficiency: float = 1.0,
        minimal_distance: float = 0.0,
        secondary = False,
        info="NA",
        **kwargs,
    ):
        """
        Targets is a dictionary with the fluorophore name as key
        Each of these keys have a set of emitters in the case
        of antibodies or linkers. For direct labeling the key
        value is None
        """
        # should check if the name is different
        label_name = label_name
        self.labelnames.append(label_name)
        for (
            keys,
            values,
        ) in targets.items():  # assumption: theres only one key and value
            fluorophore = keys
            emitters = values
        if emitters is None:
            label_type = "direct"
        else:
            label_type = "indirect"
        colour = set_colorplot(self.plotting_params)
        # print(f"colour assigned to label: {colour}")
        plt_params = dict(plotsize=20, plotalpha=1, plotmarker="o", plotcolour=colour)
        self.plotting_params[label_name] = plt_params
        label_params = dict(
            emitters=emitters,
            scale=scale,
            axis=axis,
            labelling_efficiency=labellig_efficiency,
            minimal_distance=minimal_distance,
            label_type=label_type,
            fluorophore=fluorophore,
        )
        if secondary:
            self.secondary[label_name] = label_params
        else:
            self.labels[label_name] = label_params
            self.status["labels"] = True

    def _get_labels(self):
        return self.labels

    def get_label_names(self):
        return self.labelnames

    def _get_label_fluorophore(self, label_name):
        return self.labels[label_name]["fluorophore"]

    def _get_label_emitters(self, label_name):
        # by default we consider each label has
        # a single type of fluorophore conjugated to it
        return self.labels[label_name]["emitters"]

    def _get_label_efficiency(self, label_name):
        return self.labels[label_name]["labelling_efficiency"]

    def _get_label_axis(self, label_name):
        return self.labels[label_name]["axis"]

    def _get_label_minimal_dist(self, label_name):
        return self.labels[label_name]["minimal_distance"]

    def _get_label_type(self, label_name):
        return self.labels[label_name]["label_type"]

    def _get_label_plotting_params(self, label_name):
        return self.plotting_params[label_name]

    # methods to change label parameters
    def _change_label_fluorophore(self, label_name, newfluorophore):
        pass

    def _change_label_efficiency(self, label_name, newefficiency):
        self.labels[label_name]["labelling_efficiency"] = newefficiency

    # methods to generate a labelled instance from a target
    def _label_source_target(self, target_name):
        """
        Given a label name, it will consider the type of labeling
        and will create an instance of emitter locations following the source
        and the label characteristics

        Method inteded to be called for each of the targets in the Source
        """
        # The target name should correspond to the label name
        # print(f"Generating instance for label {target_name}")

        if (
            self._get_source_target_normals(target_name) is None
        ):  # if no normals exist this is true
            target_type = "direct"
            # print(f"Label has no normals defined. Using as direct labelling")
        else:
            target_type = "indirect"
            # print(f"Label has normals defined. Using as indirect labelling")
        target_normals = self._get_source_coords_normals(target_name)
        ###### print(f"target_normals before breaking: {target_normals}")
        label4target = self._get_source_target_label(target_name)
        # at this point we can sample pairs of target_normals to model defects
        if self.defects:
            defect_target_normals = self._model_defects(target_normals)
            target_normals = copy.copy(defect_target_normals)
        else:
            defect_target_normals = None
        ##### print(f"target_normals after breaking: {target_normals}")
        # print(target_normals, label4target)
        labelling_realisation = None
        plotting_params = None
        fluorophore_name = None
        labelling_realisation_vectors = None
        if label4target is not None:
            labelling_realisation, labelling_realisation_vectors = create_instance_label(
                target_normals, target_type, label4target
            )
            plotting_params = self._get_label_plotting_params(target_name)
            fluorophore_name = self._get_label_fluorophore(target_name)
        else:
            warnings.warn(
                f'No Label correspond to target "{target_name}". '
                "Target will no be considered"
            )
        # generate the emitters positions
        # create_instance_label is capable of cosidering direct and indirect labelling
        # Once created, define the instance
        return (
            labelling_realisation,
            labelling_realisation_vectors,
            fluorophore_name,
            plotting_params,
            defect_target_normals,
        )

    def _generate_instance_constructor(self):
        targets_labeled_instance = {}
        targets_labeled_instance_vectors = {}
        label_fluorophore = {}
        plotting_params = {}
        # for each target create the emitters
        for target_name in self._get_source_target_names():
            # print(target_name)
            self.defects_target_normals = None
            emitters, labelling_realisation_vectors, fluorophore_name, plotting_par, defects_target_normals = (
                self._label_source_target(target_name)
            )
            if defects_target_normals is not None:
                self.defects_target_normals = defects_target_normals
            # print(f"Emitters for target {target_name}")
            # print(emitters)
            targets_labeled_instance[target_name] = emitters
            targets_labeled_instance_vectors[target_name] = labelling_realisation_vectors
            label_fluorophore[target_name] = fluorophore_name
            plotting_params[target_name] = plotting_par
        # target_labeled_instance
        # then crete the instance constructor dictionary
        instance_constructor = dict(
            emitters_dictionary=targets_labeled_instance,
            emitters_vectors = targets_labeled_instance_vectors,
            ref_point=self._get_source_parameter("reference_pt"),
            scale=self._get_source_parameter("scale"),
            axis=self._get_source_parameter("axis"),
            label_fluorophore=label_fluorophore,
            plotting_params=plotting_params,
        )
        return instance_constructor
    

    def _generate_primary_epitopes(
            self,
            emitters_dictionary,
            emitters_vectors,
            **kwargs
        ):
        # label the soruces using only the primary
        # create another sequential source
        for target_name in emitters_dictionary.keys():
            col_vec = np.array([i for i in emitters_vectors[target_name]])
            self.primary["targets"][target_name] = dict(
                coordinates=emitters_dictionary[target_name],
                normals = col_vec
            )

    def _label_primary(self, target_name):
        emitters = None
        target_normals = self.primary["targets"][target_name]
        if target_name in self.secondary.keys():
            label4target = self.secondary[target_name]
            print(target_normals.keys())
            print(label4target)
            labelling_realisation, labelling_realisation_vectors = create_instance_label(
                    target_normals, "indirect", label4target
                )
            return emitters
        else:
            print(f"No secondary for {target_name}")

    def _generate_secondary_labelling(self):
        targets_labeled_instance = {}
        targets_labeled_instance_vectors = {}
        label_fluorophore = {}
        plotting_params = {}
        # for each target create the emitters
        for target_name in self.primary["targets"].keys():
            emitters = self._label_primary(target_name)

        pass


    def generate_instance(self):
        if self.sequential_labelling:
            primaries = self._generate_instance_constructor()
            self._generate_primary_epitopes(**primaries)
            self._generate_secondary_labelling()
            ##self.add_emitters_n_refpoint(**secondaries)
            return primaries
        elif self.status["source"] and self.status["labels"]:
            constructor = self._generate_instance_constructor()
            self.add_emitters_n_refpoint(**constructor)
        else:
            print("Missing source or Label info")

    def add_defects(
        self,
        deg_dissasembly=0.5,
        xmer_neigh_distance=100,
        fracture=-24,
        eps1=20,
        minsamples=1,
    ):
        """
        specify parameters for adding defects
        """
        d_cluster_params = dict(
            eps1=eps1,
            minsamples1=minsamples,
            eps2=xmer_neigh_distance,
            minsamples2=minsamples,
        )
        self.defects_params["d_cluster_params"] = d_cluster_params
        self.defects_params["deg_dissasembly"] = deg_dissasembly
        self.defects_params["xmer_neigh_distance"] = xmer_neigh_distance
        self.defects_params["fracture"] = fracture
        self.defects = True
        self.generate_instance()

    def _model_defects(self, target_normals_dictionary):
        """
        target_normals_dictionary: dictionary with keys: "coordinates"
        and "normals".
        valye of key Normals can be none.
        Returns:
            target_nomral_subset (dict):
                coordinates: numpy array of epitope coordinates
                normals: numpy array of unit vector to define nomral
                    from its corresponding epitope
        """
        # print(f"input for _model_defects: {target_normals_dictionary}")
        target_sites = target_normals_dictionary["coordinates"]
        # print(f"in model_defects: {target_sites.shape}")
        # print(target_sites, self.defects_params)
        boolean_subset = xmersubset_byclustering(
            epitopes_coords=target_sites, return_ids=True, **self.defects_params
        )
        coor = target_normals_dictionary["coordinates"][boolean_subset,]
        if target_normals_dictionary["normals"] is not None:
            nor = target_normals_dictionary["normals"][boolean_subset,]
        else:
            nor = None
        target_nomral_subset = dict(coordinates=coor, normals=nor)
        target_nomral_subset

        return target_nomral_subset

    # methods to create instance
    def _add_emitters(self, emitters_dictionary: dict):
        self.emitters = emitters_dictionary

    def _set_ref_point(self, newref_point: np.array):
        self.params["ref_point"] = newref_point

    def _set_scale(self, newscale: float):
        self.params["scale"] = newscale

    def _set_axis(self, axis: dict):
        self.axis = axis

    def _set_label_fluorophopre(self, label_fluo: dict):
        self.label2fluo = label_fluo

    def add_emitters_n_refpoint(
        self,
        emitters_dictionary: dict,
        ref_point: np.ndarray,
        scale: float,
        axis: dict,
        label_fluorophore: dict,
        plotting_params: dict,
        **kwargs
    ):
        """
        Emitters_dictionary: with the name of the label_name as keys
            and coordinates as values for each.
        """
        self.labelnames = list(emitters_dictionary.keys())
        self._add_emitters(dict(emitters_dictionary))  # making an explicit copy
        self._set_ref_point(ref_point)
        self._set_scale(scale)
        self._set_axis(dict(axis))
        self._set_label_fluorophopre(dict(label_fluorophore))
        self._gen_fluo2labels()
        self.plotting_params = dict(plotting_params)

    def get_axis(self):
        return self.axis

    def get_ref_point(self):
        return self.params["ref_point"]

    def transform_reorient(self, neworientation: np.array):
        if np.linalg.norm(neworientation) == 0:
            print(f"Norm for vector {neworientation} is 0. No reorientation done")
        else:
            thet = np.arccos(
                np.dot(self.axis["direction"], neworientation)
                / (
                    np.linalg.norm(self.axis["direction"])
                    * np.linalg.norm(neworientation)
                )
            )
            print(
                f"theta: {thet}, "
                f"new {neworientation}, "
                f"current: {self.axis['direction']}"
            )
            if np.absolute(thet) == 1:
                print(
                    f'input vector {neworientation} '
                    f'has same direction as {self.axis["direction"]}. '
                    f'No reorientation done'
                )
            else:
                nori = copy.copy(neworientation)
                # this function should take the new orientation and match the
                for trgt in self.labelnames:
                    reoriented = rotate_pts_by_vector(
                        self.get_emitter_by_target(trgt),
                        self.axis["direction"],
                        nori,
                        self.get_ref_point(),
                    )
                    self.emitters[trgt] = reoriented
                self.axis["direction"] = nori

    def transform_rotate_around_axis(self, degree: float):
        pass

    def transform_translate(self, newcenter):
        nref = copy.copy(newcenter)
        # move labels
        for labeltype in self.emitters.keys():
            self.emitters[labeltype], _ = transform_displace_set(
                self.emitters[labeltype], self.get_ref_point(), nref
            )
        # then replace reference point
        self._set_ref_point(nref)
        self.axis["pivot"] = nref

    def get_info(self):
        print(
            f"instance has {len(self.emitters.keys())} "
            f"types of emitters with names {list(self.emitters.keys())}"
        )

    def get_scale(self):
        return self.params["scale"]

    def scale_coordinates_system(self, new_scale: float):
        current_scale = self.get_scale()
        if current_scale == new_scale:
            print("Current scale and new scale are the same. No scaling performed")
        else:
            scaling_factor = current_scale / new_scale
            for labeltype in self.emitters.keys():
                # print(f'before: {self.emitters[labeltype]}')
                self.emitters[labeltype] = (
                    self.get_emitter_by_target(labeltype) * scaling_factor
                )
                # print(f'after: {self.emitters[labeltype]}')
            self.params["ref_point"] = self.params["ref_point"] * scaling_factor
            self.axis["pivot"] = self.axis["pivot"] * scaling_factor
            self._set_scale(new_scale)

    def get_emitter_by_target(self, targetname: str):
        return self.emitters[targetname]

    def _gen_fluo2labels(self):
        inv_map = {}
        for k, v in self.label2fluo.items():
            inv_map[v] = inv_map.get(v, []) + [k]
        self.fluo2labels = inv_map

    def get_emitters_by_fluorophore(self, fluoname: str):
        lab_target = self.fluo2labels[fluoname]
        pulled = np.empty((0, 3))
        for l in lab_target:
            pulled = np.vstack([pulled, self.get_emitter_by_target(l)])
        return pulled

    # Visualisation methods
    def show_instance(
        self,
        labelnames="All",
        reference_point=False,
        view_init=[20, 0, 0],
        axesoff=True,
        show_axis=False,
        with_sources=False,
        return_plot=False,
    ):
        if labelnames == "All":
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            for labs in self.labelnames:
                add_ax_scatter(
                    ax,
                    format_coordinates(
                        self.emitters[labs], **self._get_label_plotting_params(labs)
                    ),
                )
                if with_sources:
                    add_ax_scatter(
                        ax,
                        format_coordinates(
                            self._get_source_coords_normals(labs)["coordinates"]
                        ),
                    )
                # add_ax_scatter(ax, format_coordinates(self.emitters[labs]))
        if reference_point:
            ref = self.get_ref_point()
            print(f"Reference point at {ref}")
            ax.scatter(ref[0], ref[1], ref[2], c="k", label="ref", s=20, marker="x")
        if show_axis:
            print(f'current axis direction: {self.axis["direction"]}')
            draw1nomral_segment(self.axis, ax, lenght=150, colors=["g", "y"])
        ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        ax.set_box_aspect(
            [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
        )
        if axesoff:
            ax.set_axis_off()
        else:
            ax.set_xlabel("X (Angstroms)")
            ax.set_ylabel("Y (Angstroms)")
            ax.set_zlabel("Z (Angstroms)")
        if return_plot:
            return ax
        else:
            fig.show()

    def gen_axis_plot(
        self,
        labelnames="All",
        reference_point=False,
        view_init=[20, 0, 0],
        axesoff=True,
        show_axis=False,
        with_sources=False,
        source_plotsize=1,
        axis_object=None,
        emitter_plotsize=1,
    ):
        zlims = [0, 1]
        if labelnames == "All":
            for labs in self.labelnames:
                lab_plotparams = self._get_label_plotting_params(labs)
                lab_plotparams["plotsize"] = emitter_plotsize
                add_ax_scatter(
                    axis_object,
                    format_coordinates(self.emitters[labs], **lab_plotparams),
                )
                zlims[0] = np.min(self.emitters[labs])
                zlims[1] = np.max(self.emitters[labs])
                if with_sources:
                    if self.defects_target_normals is not None:
                        add_ax_scatter(
                            axis_object,
                            format_coordinates(
                                self.defects_target_normals["coordinates"]
                            ),
                        )
                        add_ax_scatter(
                            axis_object,
                            format_coordinates(
                                self._get_source_coords_normals(labs)["coordinates"],
                                plotcolour="#bbbbbb",
                                plotalpha=0.5,
                                plotsize=source_plotsize,
                            ),
                        )
                    else:
                        add_ax_scatter(
                            axis_object,
                            format_coordinates(
                                self._get_source_coords_normals(labs)["coordinates"],
                                plotsize=source_plotsize,
                            ),
                        )

                # add_ax_scatter(ax, format_coordinates(self.emitters[labs]))
        if reference_point:
            ref = self.get_ref_point()
            print(f"Reference point at {ref}")
            axis_object.scatter(
                ref[0], ref[1], ref[2], c="k", label="ref", s=20, marker="x"
            )
        if show_axis:
            print(f'current axis direction: {self.axis["direction"]}')
            draw1nomral_segment(self.axis, axis_object, lenght=150, colors=["g", "y"])
        axis_object.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        if axesoff:
            axis_object.set_axis_off()
        else:
            fontsize = 5
            axis_object.set_zlim(zlims[0], zlims[1])
            # axis_object.set_zticks(np.arange(-200, 200, 200))
            axis_object.set_xlabel("X (Angstroms)", size=fontsize)
            axis_object.set_ylabel("Y (Angstroms)", size=fontsize)
            axis_object.set_zlabel("Z (Angstroms)", size=fontsize)
            axis_object.tick_params(axis="both", which="major", labelsize=fontsize)
        axis_object.set_box_aspect(
            [
                ub - lb
                for lb, ub in (getattr(axis_object, f"get_{a}lim")() for a in "xyz")
            ]
        )


def create_particle(source_builder=None, label_params_list=None):
    particle = LabeledInstance()
    if source_builder is not None:
        particle.load_source(**source_builder)  # define the source
    if label_params_list is not None:
        particle = add_label_params_to_particle(particle, label_params_list)
    particle.generate_instance()
    return particle


def add_label_params_to_particle(particle: LabeledInstance, label_params):
    """
    Add source parameters to an initialised particle containing Target sites.

    """
    targets = dict()
    for lab in label_params:
        # print(lab)
        fluorophore = lab["fluorophore"]
        coordinates = None
        secondary=False
        if lab["target"]["type"] == "Primary":
            secondary=True
        if "coordinates" in lab.keys():
            # print(lab["coordinates"])
            coordinates = np.array(lab["coordinates"])
            # print(coordinates, type(coordinates))
        targets[fluorophore] = coordinates
        particle.load_label(
            targets=targets,
            label_name=lab["label_name"],
            labellig_efficiency=lab["labeling_efficiency"],
            secondary=secondary
        )
        if "epitope_site" in lab.keys():
            print(f'assigning epitope site : {lab["epitope_site"]}')
            particle.source["targets"]["epitope_site"] = copy.copy(lab["epitope_site"])
    return particle
