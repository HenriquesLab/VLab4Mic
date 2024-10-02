import easy_gui_jupyter
import supramolsim
from ..workflows import *
import matplotlib.pyplot as plt
from IPython.display import clear_output
import os

global structure, particle, configuration_path, coordinates_field, exported_field
particle = None
structure = None
configuration_path = None
coordinates_field = None
exported_field = None


def set_directory():
    config_gui = easy_gui_jupyter.EasyGUI("Config")

    pck_dir = os.path.dirname(os.path.abspath(supramolsim.__file__))
    local_dir = os.path.join(pck_dir, "configuration")

    def clear(b):
        global configuration_path
        configuration_path = []
        configuration_path.clear()
        configuration_path.append(local_dir)
        print(configuration_path)

    def select_custom(b):
        global configuration_path
        configuration_path = []
        configuration_path.clear()
        configuration_path.append(config_gui["custom_path"].value)
        config_gui.save_settings()
        print(configuration_path)

    config_gui.add_textarea("custom_path", remember_value=True)
    config_gui.add_button("Set", description="Set configuration path")
    config_gui["Set"].on_click(select_custom)
    config_gui.add_button("demos", description="Use demos")
    config_gui["demos"].on_click(clear)
    config_gui.show()


def select_structure():
    global structure, structure_param, configuration_path, nostructure, active_widgets
    active_widgets = dict()
    nostructure = True
    demo_structures = []
    structure_gui = easy_gui_jupyter.EasyGUI("Structure")
    if configuration_path is None:
        print("No configuration path")
    else:
        print(configuration_path[0])
        generic_labels = []
        # get available structure IDs
        structures_info_list = dict()
        # configuration_path = configuration_path[0]
        if os.path.isdir(configuration_path[0]):
            structure_dir = os.path.join(configuration_path[0], "structures")
            for file in os.listdir(structure_dir):
                if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                    structure_params = supramolsim.load_yaml(
                        os.path.join(structure_dir, file)
                    )
                    struct_id = structure_params["id"]
                    strict_title = structure_params["title"]
                    id_title = struct_id + ": " + strict_title
                    structures_info_list[id_title] = struct_id
                    demo_structures.append(id_title)

        def select_structure(b):
            structure_id = structures_info_list[structure_gui["struct_dropdown"].value]
            global structure, structure_param, struct_ready
            # structure_id = structure_gui["struct_dropdown"].value
            structure, structure_param = supramolsim.load_structure(
                structure_id, configuration_path[0]
            )
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False

        def view_structure(b):
            plt.clf()
            clear_output()
            # structure_gui.show()
            active_widgets["Elevation"] = None
            active_widgets["Azimuthal"] = None
            active_widgets["Roll"] = None
            for widgetname in active_widgets:
                display(structure_gui[widgetname])
            fraction = structure_gui["Fraction"].value
            el = structure_gui["Elevation"].value
            az = structure_gui["Azimuthal"].value
            rll = structure_gui["Roll"].value
            structure.show_assembly_atoms(
                assembly_fraction=fraction, view_init=[el, az, rll]
            )

        def upload_file(b):
            global structure, structure_param
            filepath = structure_gui["File"].selected
            filename = structure_gui["File"].selected_filename
            structure = build_structure_cif(
                cif_file=filepath, struct_title=filename, cif_id=filename
            )
            structure_param = struct_params_format()
            print("Structure Loaded!")
            structure_gui["View"].disabled = False
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False

        def activate_demos(b):
            global active_widgets
            plt.clf()
            clear_output()
            active_widgets["struct_dropdown"] = None
            active_widgets["Select"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        def activate_upload(b):
            global active_widgets
            plt.clf()
            clear_output()
            active_widgets["File"] = None
            active_widgets["Upload"] = None
            active_widgets["Fraction"] = None
            active_widgets["View"] = None
            for widgetname in active_widgets.keys():
                display(structure_gui[widgetname])

        structure_gui.add_button("Demos", description="Load from demo structures")
        structure_gui.add_button("Fileupload", description="Upload CIF file")
        structure_gui.add_label("Select structure model to load")
        structure_gui.add_dropdown("struct_dropdown", options=demo_structures)
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
        structure_gui["Select"].on_click(select_structure)
        structure_gui["View"].on_click(view_structure)
        structure_gui["Upload"].on_click(upload_file)
        structure_gui["Demos"].on_click(activate_demos)
        structure_gui["Fileupload"].on_click(activate_upload)
        display(structure_gui["Demos"], structure_gui["Fileupload"])
        # structure_gui.show()
