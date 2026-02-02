import matplotlib.pyplot as plt
import numpy as np
import yaml
import copy

from ..utils.visualisation.matplotlib_plots import (
    add_ax_scatter,
    draw1nomral_segment,
)
from ..utils.data_format.visualisation import format_coordinates
from ..utils.io.yaml_functions import load_yaml


class Label:
    def __init__(self):
        self.params = {}
        self.params["label_name"] = None
        self.params["labeling_efficiency"] = 1
        self.params["label_type"] = None
        self.params["fluorophore"] = None
        self.params["excitation"] = None
        self.params["emission"] = None
        self.params["atoms"] = None
        self.params["residues"] = None
        self.params["chains"] = None
        self.params["position"] = None
        self.params["plotcolour"] = None

        self.params["emitters_coords"] = None
        self.params["axis"] = dict(pivot=[0, 0, 0], direction=[0, 0, 1])

        self.params["target_sequence"] = None
        self.params["summary_method"] = "average"
        self.params["length"] = 0
        self.conjugation = dict()
        self.model = dict()
        self.binding = dict()
        self.epitope = dict()

    def set_params(self, **kwargs):
        """
        Method used to set the parameters using keyword arguments.
        :param kwargs: same as self.params.keys()
        """
        for key, value in kwargs.items():
            self.params[key] = value

    def _set_model_params(self, ID=None, format=None, database=None, **kwargs):
        self.model["ID"] = ID
        self.model["format"] = format
        self.model["database"] = database
        for key, value in kwargs.items():
            self.model[key] = value

    def _set_conjugation_params(self, target: dict, DoL=None, **kwargs):
        self.conjugation["target"] = target
        self.conjugation["DoL"] = DoL
        for key, value in kwargs.items():
            self.conjugation[key] = value

    def _set_binding_params(
        self,
        efficiency,
        orientation,
        wobble_range: dict,
        distance: dict,
        paratope: str,
        **kwargs,
    ):

        self.binding["efficiency"] = efficiency
        self.binding["orientation"] = orientation
        self.binding["wobble_range"] = wobble_range
        self.binding["distance"] = distance
        self.binding["paratope"] = paratope
        for key, value in kwargs.items():
            self.binding[key] = value

    def _set_epitope_params(
        self, target=None, normals=None, site=None, **kwargs
    ):
        self.epitope["target"] = target
        self.epitope["normals"] = normals
        self.epitope["site"] = site

    def set_axis(self, pivot: list, direction: list):
        self.params["axis"]["pivot"] = pivot
        self.params["axis"]["direction"] = direction

    def set_fluorophore(self, fluo_id):
        self.params["fluorophore"] = fluo_id

    def get_name(self):
        return self.params["label_name"]

    def get_label_type(self):
        return self.params["label_type"]

    def get_target_type(self):
        return self.params["target"]["type"]

    def get_plotcolour(self):
        return self.params["plotcolour"]

    def get_fluorophore(self):
        return self.params["fluorophore"]

    def get_efficiency(self):
        return self.params["labeling_efficiency"]

    def get_param(self, param: str):
        return self.params[param]

    def generate_linker(self):
        if self.params["length"] > 0:
            single_emitter = np.array([0, 0, self.params["length"]])
            single_emitter = single_emitter[np.newaxis, :]
            self.set_emitters(single_emitter)
        else:
            print("No linker created due to linker lenght == 0")

    def set_target_sequence(self, sequence: str):
        self.params["target_sequence"] = sequence

    def set_emitters(self, emitters):
        if (np.shape(emitters))[0] < 1:
            print("input need at least 1 points")
        else:
            self.params["emitters_coords"] = emitters

    def gen_labeling_entity(self):
        labeling_emitters = copy.copy(self.params["emitters_coords"])
        p1 = np.array(self.params["axis"]["pivot"])

        # Safety check for axis direction - use default if None
        direction = self.params["axis"]["direction"]
        if direction is None:
            direction = [0, 0, 1]  # Default axis direction

        p2 = p1 + np.array(direction)
        pivots = np.array([p1, p2])
        # print(f"pivots are: {pivots}")
        if len(labeling_emitters) < 1:
            print("there are no emitters specified")
        else:
            labeling_entity = np.vstack([pivots, labeling_emitters])
        return labeling_entity

    def export_as_file(
        self, filename="customlabel.yaml", write_dir=""
    ):
        writing_dir = write_dir + filename
        print(f"Writing file in: {writing_dir}")
        # define parameters to write
        if self.params["emitters_coords"] is None:
            coordinates = None
        else:
            coordinates = self.gen_labeling_entity()
            coordinates = coordinates.tolist()
        label_params = dict(
            label_type=self.get_param("label_type"),
            label_name="NAME",
            labeling_efficiency=self.get_efficiency(),
            coordinates=coordinates,
            target_sequence=self.get_param("target_sequence"),
            plotcolour=None,
            fluorophore=self.get_param("fluorophore"),
        )
        yaml_dict = dict(
            labels=[
                label_params,
            ]
        )
        with open(writing_dir, "w") as file:
            yaml.dump(yaml_dict, file)

    def plot_emitters(
        self,
        view_init=[0, 0, 0],
        axesoff=True,
        show_axis=True,
        return_plot=False,
    ):
        if self.params["emitters_coords"] is None:
            print(
                f"Label type is {self.params['label_type']} and does not contain "
                "emitter coordinates. To include them manually use method set_emitters"
            )
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            add_ax_scatter(
                ax, format_coordinates(self.params["emitters_coords"])
            )
            ax.view_init(
                elev=view_init[0], azim=view_init[1], roll=view_init[2]
            )
            if show_axis:
                draw1nomral_segment(
                    self.params["axis"], ax, lenght=150, colors=["g", "y"]
                )
            ax = plt.gca()
            emitters = self.params["emitters_coords"]
            minimum = np.min(emitters)
            maximum = np.max(emitters)
            ax.set_xlim(minimum, maximum)
            ax.set_ylim(minimum, maximum)
            ax.set_zlim(minimum, maximum)
            if axesoff:
                ax.set_axis_off()
            ax.set_box_aspect(
                [
                    ub - lb
                    for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")
                ]
            )
            if return_plot:
                plt.close()
                return fig
            else:
                fig.show

    def get_notNone_params(self):
        for key, value in self.params.items():
            if value is not None:
                print(f"{key}: {value}")


def construct_label(
    label_config_dictionary: None,
    lab_eff: float = None,
    target_info=None,
):
    """
    Construct an object of class Label.
    Assumes there exist a configuration file

    Args:
        label_id: (string) ID of Label configuration file
        fluorophore_id: (string) ID of fluorophore configuration file
    Output:
        label: (Label) Object of class Label()
        label_params: (dictionary)
        fluorophore_params: (dictionary)
    """
    label_params = copy.deepcopy(label_config_dictionary)

    # keynames reserved for parameters that can be iteratet over
    if "model_ID" in label_config_dictionary.keys():
        label_params["model"]["ID"] = label_config_dictionary["model_ID"]
    if "distance_to_epitope" in label_params.keys():
        label_params["binding"]["distance"]["to_target"] = label_params[
            "distance_to_epitope"
        ]
    else:
        print("No distance to epitope provided, using default value.")
    if "distance_between_epitope" in label_config_dictionary.keys():
        label_params["binding"]["distance"]["between_targets"] = label_config_dictionary[
            "distance_between_epitope"
        ]
    if "paratope" in label_config_dictionary.keys():
        label_params["binding"]["paratope"] = label_config_dictionary["paratope"]
    if "conjugation_target_info" in label_config_dictionary.keys():
        label_params["conjugation_sites"]["target"] = label_config_dictionary[
            "conjugation_target_info"
        ]
    #if "probe_DoL" in label_params.keys():
    #    label_params["conjugation_sites"]["DoL"] = kwargs[
    #        "probe_DoL"
    #    ]
    #else:
    #    print("No DoL provided, using default value. ########################")
    if "epitope_target_info" in label_config_dictionary.keys():
        label_params["epitope"]["target"] = label_config_dictionary["epitope_target_info"]
    if "wobble_theta" in label_params.keys():
        label_params["binding"]["wobble_range"]["theta"] = label_params[
            "wobble_theta"
        ]
    else:
        label_params["binding"]["wobble_range"]["theta"] = None
    ######## information about target
    if target_info:
        # expect label_params to have empty values on target type and value
        label_params["target"]["type"] = target_info["type"]
        label_params["target"]["value"] = target_info["value"]
    if label_params["target"]["type"] == "Atom_residue":
        label_params["atoms"] = [
            label_params["target"]["value"]["atoms"],
        ]
        label_params["residues"] = [
            label_params["target"]["value"]["residues"],
        ]
        try:
            label_params["position"] = label_params["target"]["value"][
                "position"
            ]
        except:
            label_params["position"] = None
    ######## Building the labelling entity: anribody, linker, direct...
    if "fluorophore_id" not in label_params.keys():
        fluorophore_id = "AF647"
    else:
        fluorophore_id = label_params["fluorophore_id"]
    label_params["fluorophore"] = fluorophore_id
    if lab_eff is not None:
        label_params["labelling_efficiency"] = lab_eff
    label = Label()
    label.set_params(**label_params)
    label.set_fluorophore(fluorophore_id)
    # check whether there is enough structural information for the labelling entity
    if (
        label_params["model"]["ID"]
        and label_params["conjugation_sites"]["target"]["type"]
    ):
        # building probe from PDB/CIF
        print("antibody")
        label._set_model_params(**label_params["model"])
        label._set_conjugation_params(**label_params["conjugation_sites"])
        label._set_binding_params(**label_params["binding"])
        if "epitope" in label_params.keys():
            label._set_epitope_params(**label_params["epitope"])
        else:
            label._set_epitope_params()
        # label.emitters_from_PDBCIF(**label_params)
        # label_params["coordinates"] = label.gen_labeling_entity()
    elif label_params["binding"]["distance"]["to_target"]:
        # not from a PDB, checking if at least distance to make a rigid linker
        print("rigid linker")
        label.set_params(
            length=label_params["binding"]["distance"]["to_target"]
        )
        label.set_axis(
            pivot=[0, 0, 0], direction=label_params["binding"]["orientation"]
        )
        label.generate_linker()
        label_params["coordinates"] = label.gen_labeling_entity()
    else:
        print("direct")
    # if none of these conditions exist, it will be treated as direct
    return label, label_params
