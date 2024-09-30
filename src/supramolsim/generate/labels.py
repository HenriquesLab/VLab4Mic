import matplotlib.pyplot as plt
import numpy as np
import yaml
import copy

# from simsmlm.utils.visualisation.plots import add_ax_scatter

from ..utils.visualisation.matplotlib_plots import add_ax_scatter, draw1nomral_segment
from ..utils.data_format.visualisation import format_coordinates


class Label():
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
        self.params["axis"] = dict(pivot=None, direction=None)

        self.params["target_sequence"] = None
        self.params["summary_method"] = "average"
        self.params["length"] = 0

    def set_params(self, **kwargs):
        """
        Method used to set the parameters using keyword arguments.
        :param kwargs: same as self.params.keys()
        """
        for key, value in kwargs.items():
            self.params[key] = value

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
        return self.params["target_type"]

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

    def set_emitters(self, emitters):  # TODO: test
        if (np.shape(emitters))[0] < 1:
            print(f'input need at least 1 points')
        else:
            self.params["emitters_coords"] = emitters

    def gen_labeling_entity(self):
        labeling_emitters = copy.copy(self.params["emitters_coords"])
        p1 = np.array(self.params["axis"]["pivot"])
        p2 = p1 + np.array(self.params["axis"]["direction"])
        pivots = np.array([p1, p2])
        print(f"pivots are: {pivots}")
        if len(labeling_emitters) < 1:
            print("there are no emitters specified")
        else:
            labeling_entity = np.vstack([pivots, labeling_emitters])
        return labeling_entity

    def export_as_file(self, filename="customlabel.yaml", write_dir=""):  #TODO: consider using get_notNone_params method
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
            fluorophore=self.get_param("fluorophore")
        )
        yaml_dict = dict(
            labels=[label_params, ]
        )
        with open(writing_dir, 'w') as file:
            yaml.dump(yaml_dict, file)

    def plot_emitters(self, view_init=[0,0,0], axesoff=True, show_axis=True):  # TODO: test
        if self.params["emitters_coords"] is None:
            print(f"Label type is {self.params['label_type']} and does not contain emitter coordinates. To include them manually use method set_emitters.")
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            add_ax_scatter(ax, format_coordinates(self.params["emitters_coords"]))
            ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
            if show_axis:
                draw1nomral_segment(self.params["axis"], ax, lenght=150, colors=['g', 'y'])
            ax = plt.gca()
            emitters = self.params["emitters_coords"]
            minimum = np.min(emitters)
            maximum = np.max(emitters)
            ax.set_xlim(minimum, maximum)
            ax.set_ylim(minimum, maximum)
            ax.set_zlim(minimum, maximum)
            if axesoff:
                ax.set_axis_off()
            ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
            fig.show

    def get_notNone_params(self):
        for key, value in self.params.items():
            if value is not None:
                print(f'{key}: {value}')
