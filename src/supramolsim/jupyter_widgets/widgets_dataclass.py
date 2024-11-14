from dataclasses import dataclass, field, fields
from typing import List, Dict
from ..experiments import ExperimentParametrisation
import easy_gui_jupyter
import os
from ..utils import data_format
from ..utils.io.yaml_functions import load_yaml
import matplotlib.pyplot as plt
from IPython.display import clear_output
from ..generate.molecular_structure import build_structure_cif


@dataclass
class jupyter_gui:
    my_experiment = ExperimentParametrisation()

    def __post_init__(self):
        self.demo_structures = []
        # get available structure IDs
        self.structures_info_list = dict()
        structure_dir = os.path.join(
            self.my_experiment.configuration_path, "structures"
        )
        for file in os.listdir(structure_dir):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                structure_params = load_yaml(os.path.join(structure_dir, file))
                struct_id = structure_params["id"]
                strict_title = structure_params["title"]
                id_title = struct_id + ": " + strict_title
                self.structures_info_list[id_title] = struct_id
                self.demo_structures.append(id_title)

    def structure_gui(self):
        structure_gui = easy_gui_jupyter.EasyGUI("Structure")
        active_widgets = dict()
        nostructure = True

        def select_structure(b):
            structure_id = self.structures_info_list[
                structure_gui["struct_dropdown"].value
            ]
            self.my_experiment.structure_id = structure_id
            self.my_experiment._build_structure()
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False

        def view_structure(b):
            plt.clf()
            clear_output()
            active_widgets["Elevation"] = None
            active_widgets["Azimuthal"] = None
            active_widgets["Roll"] = None
            for widgetname in active_widgets:
                display(structure_gui[widgetname])
            fraction = structure_gui["Fraction"].value
            el = structure_gui["Elevation"].value
            az = structure_gui["Azimuthal"].value
            rll = structure_gui["Roll"].value
            axes_off = structure_gui["hidegrid"].value
            self.my_experiment.structure.show_assembly_atoms(
                assembly_fraction=fraction, view_init=[el, az, rll], axesoff=axes_off
            )

        def upload_file(b):
            global structure, structure_param, structure_id
            filepath = structure_gui["File"].selected
            filename = structure_gui["File"].selected_filename
            structure = build_structure_cif(
                cif_file=filepath, struct_title=filename, cif_id=filename
            )
            structure_param = data_format.structural_format.struct_params_format()
            print("Structure Loaded!")
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False
            structure_id = filename.split(".")[0]

        def activate_demos(b):
            plt.clf()
            clear_output()
            active_widgets["struct_dropdown"] = None
            active_widgets["Select"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            active_widgets["hidegrid"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        def activate_upload(b):
            plt.clf()
            clear_output()
            active_widgets["File"] = None
            active_widgets["Upload"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            active_widgets["hidegrid"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        structure_gui.add_button("Demos", description="Load from demo structures")
        structure_gui.add_button("Fileupload", description="Upload CIF file")
        structure_gui.add_label("Select structure model to load")
        structure_gui.add_dropdown("struct_dropdown", options=self.demo_structures)
        structure_gui.add_button("Select", description="Load structure")
        structure_gui.add_file_upload(
            "File", description="Select from file", accept="*.cif", save_settings=False
        )
        structure_gui.add_button("Upload", description="Load from file")
        structure_gui.add_button(
            "View", description="Show structure", disabled=nostructure
        )
        structure_gui.add_float_slider(
            "Fraction",
            value=0.01,
            min=0,
            max=1,
            step=0.001,
            description="Fraction of total atoms to show",
            disabled=nostructure,
            readout_format=".3f",
        )

        structure_gui.add_int_slider(
            "Elevation",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Elevation",
            disabled=nostructure,
        )
        structure_gui.add_int_slider(
            "Azimuthal",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Azimuthal",
            disabled=nostructure,
        )
        structure_gui.add_int_slider(
            "Roll",
            value=0,
            min=-90,
            max=90,
            step=1,
            description="Roll",
            disabled=nostructure,
        )
        structure_gui.add_checkbox(
            "hidegrid", description="Hide plot grids", value=True, disabled=nostructure
        )
        structure_gui["Select"].on_click(select_structure)
        structure_gui["View"].on_click(view_structure)
        structure_gui["Upload"].on_click(upload_file)
        structure_gui["Demos"].on_click(activate_demos)
        structure_gui["Fileupload"].on_click(activate_upload)
        display(structure_gui["Demos"], structure_gui["Fileupload"])
