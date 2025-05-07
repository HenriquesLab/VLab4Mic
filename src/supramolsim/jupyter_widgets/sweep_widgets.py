import supramolsim.utils
from ..sweep_generator import sweep_generator
import ipywidgets as widgets
from .widget_generator import widgen
from .widgets_dataclass import jupyter_gui
import easy_gui_jupyter
import os
import supramolsim
from supramolsim.utils.io import yaml_functions
import copy

class Sweep_gui(jupyter_gui):
    sweep_gen = sweep_generator()
    wgen = widgen()
    selected_structure = None
    selected_probes = None
    selected_modalities = None
    vlab_probes = []
    models4probes = []
    probes_per_structure = {}
    probe_parameters = {}
    # widget parameters
    param_settings = None
    
    def __init__(self):
        super().__init__()
        self.vlab_probes = []
        param_settings_file = os.path.join(self.config_directories["base"], "parameter_settings.yaml")
        self.param_settings = yaml_functions.load_yaml(param_settings_file)
        for file in os.listdir(self.config_directories["probes"]):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                label_config_path = os.path.join(self.config_directories["probes"], file)
                label_parmeters = supramolsim.load_yaml(label_config_path)
                lablname = os.path.splitext(file)[0]
                #self.vlab_probes.append(lablname)
                #self.probe_parameters[lablname] = label_parmeters
                if "Mock" in label_parmeters["known_targets"]:
                    self.models4probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                else:
                    self.vlab_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                    for known_structures in label_parmeters["known_targets"]:
                        if known_structures in self.probes_per_structure.keys():
                            self.probes_per_structure[known_structures].append(lablname)
                        else:
                            self.probes_per_structure[known_structures] = [lablname,]
    
    
    
    def select_structure(self):
        ez_sweep_structure = easy_gui_jupyter.EasyGUI("structure")
        ez_sweep_structure.add_dropdown("structures", options=self.demo_structures)
        ez_sweep_structure.add_button("Select", description="Select")
        def select(b):
            self.selected_structure = self.structures_info_list[
                ez_sweep_structure["structures"].value
            ]
            ez_sweep_structure["structures"].disabled = True
        ez_sweep_structure["Select"].on_click(select)
        ez_sweep_structure.show()

    
    def select_probes_and_mods(self):
        ez_sweep = easy_gui_jupyter.EasyGUI("Sweep")
        probes2show = []
        if self.selected_structure:
            for p in self.probes_per_structure[self.selected_structure]:
                probes2show.extend(copy.copy(self.probes_per_structure[self.selected_structure]))
            #probes2show + copy.copy(self.probes_per_structure[self.selected_structure])
        probes2show.extend(copy.copy(self.models4probes))
        #ez_sweep.add_dropdown("probes", description="Probes", options=self.vlab_probes)
        ez_sweep._widgets["probes"]= widgets.SelectMultiple(
            description="probes", 
            options=probes2show
        )
        ez_sweep._widgets["modalities"]= widgets.SelectMultiple(
            description="Modalities", 
            options=self.modalities_default
        )
        # on clicks
        def select_str(b):
            self.selected_modalities = ez_sweep["modalities"].value
            self.selected_probes = ez_sweep["probes"].value
            ez_sweep["modalities"].disabled = True
            ez_sweep["probes"].disabled = True
        ez_sweep.add_button("Select", description="Select")
        ez_sweep["Select"].on_click(select_str)
        ez_sweep.show()
        
    def add_parameters_ranges(self):
        param_ranges = easy_gui_jupyter.EasyGUI("ranges")
        
        def change_param_list(change):
            new_options = list(self.param_settings[change.new].keys())
            param_ranges["parms_per_group"].options = new_options

        parameter_group_names = list(self.param_settings.keys())
        #print(self.param_settings)
        param_ranges.add_dropdown("groups", options = parameter_group_names)
        param_ranges.add_dropdown(
            "parms_per_group",
            options = list(self.param_settings[param_ranges["groups"].value].keys()))
        param_ranges["groups"].observe(change_param_list, names='value')
        param_ranges.show()

