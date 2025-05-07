from ..sweep_generator import sweep_generator
import ipywidgets as widgets
from .widget_generator import widgen
from .widgets_dataclass import jupyter_gui
import easy_gui_jupyter
import os
import supramolsim

class Sweep_gui(jupyter_gui):
    sweep_gen = sweep_generator()
    wgen = widgen()
    selected_structure = None
    selected_modalities = None
    vlab_probes = []
    model_probes = []
    probe_parameters = {}
    
    def __init__(self):
        super().__init__()
        self.vlab_probes = []
        for file in os.listdir(self.config_directories["probes"]):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                label_config_path = os.path.join(self.config_directories["probes"], file)
                label_parmeters = supramolsim.load_yaml(label_config_path)
                lablname = os.path.splitext(file)[0]
                #self.vlab_probes.append(lablname)
                #self.probe_parameters[lablname] = label_parmeters
                if "Mock" in label_parmeters["known_targets"]:
                    self.model_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
                else:
                    self.vlab_probes.append(lablname)
                    self.probe_parameters[lablname] = label_parmeters
    
    
    
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

    
    def static_selection(self):
        ez_sweep = easy_gui_jupyter.EasyGUI("Sweep")
        #ez_sweep.add_dropdown("probes", description="Probes", options=self.vlab_probes)
        ez_sweep._widgets["probes"]= widgets.SelectMultiple(
            description="Modalities", 
            options=self.vlab_probes
        )
        ez_sweep._widgets["modalities"]= widgets.SelectMultiple(
            description="Modalities", 
            options=self.modalities_default
        )
        # on clicks
        def select_str(b):

            self.selected_modalities = ez_sweep["modalities"].value
            ez_sweep["modalities"].disabled = True
            ez_sweep["probes"].disabled = True
        ez_sweep.add_button("Select", description="Select")
        ez_sweep["Select"].on_click(select_str)
        ez_sweep.show()
        
    

