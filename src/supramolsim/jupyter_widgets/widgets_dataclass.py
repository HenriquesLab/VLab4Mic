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
import supramolsim
import ipywidgets as widgets
from IPython.utils import io
from supramolsim.workflows import create_imaging_system
import numpy as np
import tifffile as tif
import copy
import mpl_toolkits.axes_grid1 as axes_grid1

@dataclass
class jupyter_gui:
    my_experiment = ExperimentParametrisation()
    structures_to_show = ["3J3Y", "7R5K", "1XI5"]
    modalities_default = ["Widefield", "Confocal", "STED", "SMLM"]
    def __post_init__(self):
        self.demo_structures = []
        # get available structure IDs
        self.structures_info_list = dict()
        structure_dir = os.path.join(
            self.my_experiment.configuration_path, "structures"
        )
        fluorophores_dir = os.path.join(
            self.my_experiment.configuration_path, "fluorophores"
        )
        probes_dir = os.path.join(self.my_experiment.configuration_path, "probes")
        modalities_dir = os.path.join(
            self.my_experiment.configuration_path, "modalities"
        )
        for file in os.listdir(structure_dir):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                structure_params = load_yaml(os.path.join(structure_dir, file))
                struct_id = structure_params["model"]["ID"]
                if struct_id in self.structures_to_show:
                    strict_title = structure_params["model"]["title"]
                    id_title = struct_id + ": " + strict_title
                    self.structures_info_list[id_title] = struct_id
                    self.demo_structures.append(id_title)
        self.config_directories = dict(
            structure=structure_dir,
            fluorophores=fluorophores_dir,
            probes=probes_dir,
            modalities=modalities_dir,
            base=self.my_experiment.configuration_path,
        )
    
    def _update_widgets(self, container, visibility):
        for widgetname in visibility.keys():
                if visibility[widgetname]:
                    container[widgetname].layout.display = 'inline-flex'
                else:
                    container[widgetname].layout.display = 'None'

    def structure_gui(self):
        structure_gui = easy_gui_jupyter.EasyGUI("Structure")
        active_widgets = dict()
        nostructure = True

        def select_structure(b):
            structure_gui["Select"].disabled = True
            widgets_visibility["label_2"] = True
            #plt.clf()
            #clear_output()
            self._update_widgets(structure_gui, widgets_visibility)
            
            structure_id = self.structures_info_list[
                structure_gui["struct_dropdown"].value
            ]
            self.my_experiment.structure_id = structure_id
            self.my_experiment._build_structure()
            widgets_visibility["label_3"] = True
            widgets_visibility["View"] = True
            widgets_visibility["label_2"] = False
            #plt.clf()
            #clear_output()
            self._update_widgets(structure_gui, widgets_visibility)
            structure_gui._widgets["label_2"].layout = widgets.Layout()
            structure_gui._widgets["label_2"].layout.display = "None"
            structure_gui["View"].disabled = False

            

        def view_structure(b):
            #plt.clf()
            #clear_output()
            structure_gui["Fraction"].disabled = False
            structure_gui["Elevation"].disabled = False
            structure_gui["Azimuthal"].disabled = False
            structure_gui["Roll"].disabled = False
            structure_gui["hidegrid"].disabled = False
            widgets_visibility["label_5"] = True
            widgets_visibility["Fraction"] = True
            widgets_visibility["label_4"] = True
            widgets_visibility["hidegrid"] = True
            widgets_visibility["Elevation"] = True
            widgets_visibility["Azimuthal"] = True
            widgets_visibility["Roll"] = True
            widgets_visibility["image_output"] = True
            self._update_widgets(structure_gui, widgets_visibility)
            fraction = structure_gui["Fraction"].value
            el = structure_gui["Elevation"].value
            az = structure_gui["Azimuthal"].value
            rll = structure_gui["Roll"].value
            axes_off = structure_gui["hidegrid"].value
            structure_gui["image_output"].clear_output()
            with structure_gui["image_output"]:
                #fig, ax = plt.subplots()
                display(
                    self.my_experiment.structure.show_assembly_atoms(
                        assembly_fraction=fraction, view_init=[el, az, rll], axesoff=axes_off
                    )
                )

        def upload_file(b):
            structure_gui["Upload"].disabled = True
            widgets_visibility["label_2"] = True
            #plt.clf()
            #clear_output()
            self._update_widgets(structure_gui, widgets_visibility)
            filepath = structure_gui["File"].selected
            filename = structure_gui["File"].selected_filename
            structure = build_structure_cif(
                cif_file=filepath, struct_title=filename, cif_id=filename
            )
            self.my_experiment.structure = structure
            self.my_experiment.objects_created["structure"] = True
            structure_param = data_format.structural_format.struct_params_format()
            widgets_visibility["label_3"] = True
            widgets_visibility["View"] = True
            #plt.clf()
            #clear_output()
            widgets_visibility["label_2"] = False
            self._update_widgets(structure_gui, widgets_visibility)
            structure_gui._widgets["label_2"].layout = widgets.Layout()
            structure_gui._widgets["label_2"].layout.display = "None"
            structure_gui["View"].disabled = False

        def activate_demos(b):
            #plt.clf()
            #clear_output()
            widgets_visibility["Demos"] = False
            widgets_visibility["Fileupload"] = False
            #
            widgets_visibility["label_1"] = True
            widgets_visibility["struct_dropdown"] = True
            widgets_visibility["Select"] = True
            #active_widgets["label_1"] = None
            #active_widgets["struct_dropdown"] = None
            #active_widgets["Select"] = None
            self._update_widgets(structure_gui, widgets_visibility)
                #display(structure_gui[widgetname])

        def activate_upload(b):
            #plt.clf()
            #clear_output()
            widgets_visibility["Demos"] = False
            widgets_visibility["Fileupload"] = False
            widgets_visibility["label_6"] = True
            widgets_visibility["File"] = True
            widgets_visibility["Upload"] = True
            self._update_widgets(structure_gui, widgets_visibility)

        structure_gui.add_button("Demos", description="Example structures")
        structure_gui.add_button("Fileupload", description="Upload CIF file")
        structure_gui.add_label("Select an example structure from the dropdown menu")
        structure_gui.add_label("Loading structure. This will take a few seconds...") # "Text_loading"]
        structure_gui.add_label("Structure Loaded! Use the panels below to visualise your structure.") # "Text_loading"]
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
            description="Fraction of atoms to show",
            disabled=nostructure,
            readout_format=".3f",
        )
        structure_gui.add_label("(A smaller fraction gets a faster rendering)")
        structure_gui.add_label("Click again in Show structure to update the figure")
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
        structure_gui.add_label("Locate your CIF file. Then click Load from file.")
        structure_gui["Select"].on_click(select_structure)
        structure_gui["View"].on_click(view_structure)
        structure_gui["Upload"].on_click(upload_file)
        structure_gui["Demos"].on_click(activate_demos)
        structure_gui["Fileupload"].on_click(activate_upload)
        structure_gui._widgets["image_output"] = widgets.Output()
        widgets_visibility = {}
        for wgt in structure_gui._widgets.keys():
            widgets_visibility[wgt] = False
            structure_gui._widgets[wgt].layout = widgets.Layout(width="50%", display="None")
        #
        structure_gui.show()
        widgets_visibility["Demos"] = True
        widgets_visibility["Fileupload"] = True
        self._update_widgets(structure_gui, widgets_visibility)
        #structure_gui["Demos"].layout.display = 'inline-flex'
        #structure_gui["Fileupload"].layout.display = 'inline-flex'
        #display(structure_gui["Demos"], structure_gui["Fileupload"])

    def structural_model_gui(self):
        # ensure that structure has no labels associated
        labels_gui = easy_gui_jupyter.EasyGUI("Labels")
        particle_created = False
        current_labels = dict()
        probe_widgets = dict()
        #generic_labels = []
        #structure_labels = []
        vlab_probes = []
        mock_labels = []
        fluorophores_list = []
        probe_parameters = {}
        structure_id = self.my_experiment.structure_id
        for fluoid in os.listdir(self.config_directories["fluorophores"]):
            if os.path.splitext(fluoid)[-1] == ".yaml" and "_template" not in fluoid:
                fluorophores_list.append(os.path.splitext(fluoid)[0])
        for file in os.listdir(self.config_directories["probes"]):
            if os.path.splitext(file)[-1] == ".yaml" and "_template" not in file:
                label_config_path = os.path.join(self.config_directories["probes"], file)
                label_parmeters = supramolsim.load_yaml(label_config_path)
                #print(label_parmeters)
                lablname = os.path.splitext(file)[0]
                if "Mock" in label_parmeters["known_targets"]:
                    mock_labels.append(lablname)
                    probe_parameters[lablname] = label_parmeters
                elif structure_id in label_parmeters["known_targets"]:
                    vlab_probes.append(lablname)
                    probe_parameters[lablname] = label_parmeters
                elif "Generic" in label_parmeters["known_targets"]:
                    vlab_probes.append(lablname)
                    probe_parameters[lablname] = label_parmeters

        def build_label(b):
            label_id = labels_gui["label_dropdown"].value
            if label_id == "<None>":
                print("Invalid label, no label added")
            else:
                fluorophore_id = labels_gui["fluo_dropdown"].value
                lab_eff = labels_gui["Labelling_efficiency"].value
                tmp_label = {
                    "probe_name": label_id,
                    "probe_fluorophore": fluorophore_id,
                    "labelling_efficiency": lab_eff
                }
                unique_name = label_id + "_conjugated_" + fluorophore_id
                if unique_name in current_labels.keys():
                    print("label already exist")
                else:
                    current_labels[unique_name] = tmp_label
                    print(f"label added: {unique_name}")
                labels_gui["Label"].disabled = False

        def build_mock_label(b):
            label_id = labels_gui["mock_label_dropdown"].value
            fluorophore_id = labels_gui["mock_fluo_dropdown"].value
            lab_eff = labels_gui["mock_Labelling_efficiency"].value
            probe_distance_to_epitope = labels_gui["mock_distance_to_epitope"].value 
            mock_type = options_dictionary[labels_gui["mock_type"].value]
            if mock_type == "Atom_residue":
                #atom, residue = labels_gui["mock_value"].value.split(".")
                residue = labels_gui["mock_type_options1"].value
                atom = labels_gui["mock_type_options2"].value
                target_value = dict(
                    atoms=atom,
                    residues=residue
                )
            elif mock_type == "Sequence": 
                #target_value = labels_gui["mock_value"].value
                chain_name = labels_gui["mock_type_options1"].value
                position = labels_gui["mock_type_options2"].value
                protein_name, _1, site, sequence = self.my_experiment.structure.get_peptide_motif(
                    chain_name=chain_name,
                    position=position)
                target_value = sequence
            elif mock_type == "Primary":
                target_value = labels_gui["mock_type_options1"].value
            target_info=dict(
                type=mock_type,
                value=target_value
            )
            as_linker = labels_gui["as_linker"].value
            if as_linker:
                options_per_type1["Primary_Probe"] = [label_id,]
                #labels_gui["mock_type_options1"].options = options_per_type1

            enable_wobble = labels_gui["wobble"].value
            wobble_theta = labels_gui["wobble_theta"].value

            tmp_label = {
                    "probe_name": label_id,
                    "probe_fluorophore": fluorophore_id,
                    "labelling_efficiency": lab_eff,
                    "probe_target_type": target_info["type"],
                    "probe_target_value": target_info["value"],
                    "probe_distance_to_epitope": probe_distance_to_epitope,
                    "as_primary": as_linker,
                    "probe_wobbling": enable_wobble,
                    "wobble_theta": wobble_theta
            }

            unique_name = label_id + "_conjugated_" + fluorophore_id
            if unique_name in current_labels.keys():
                print("label already exist")
            else:
                current_labels[unique_name] = tmp_label
                print(f"label added: {unique_name}")
            labels_gui["Label"].disabled = False

        def clear(b):
            current_labels.clear()

        def show(b):
            for lab in current_labels.keys():
                print(lab)

        def label_struct(b):
            if len(current_labels.keys()) > 0:
                self.nlabels = len(current_labels)
                for keys, values in current_labels.items():
                    self.my_experiment.add_probe(**values)
                self.my_experiment.build(modules=["particle",])
                probe_names = list(self.my_experiment.probe_parameters.keys())
                #display(probe_names)
                labels_gui["Label"].disabled = True
                labels_gui["Add"].disabled = True
                labels_gui["Add_mock"].disabled = True
            else:
                print("No label has been added")

        # buttons to select preset probes or customise
        def activate_demos(b):
            #plt.clf()
            #clear_output()
            probe_widgets["label_dropdown"] = True
            probe_widgets["label_message"] = True
            probe_widgets["label_1"] = True
            probe_widgets["Labelling_efficiency"] = True
            probe_widgets["Add"] = True
            probe_widgets["label_7"] = True
            probe_widgets["Label"] = True
            self._update_widgets(labels_gui, probe_widgets)
            
            #probe_widgets["label_dropdown"] = True
            #probe_widgets["fluo_dropdown"] = True

        def update_label_message(change):
            probe_type = probe_parameters[change.new]["target"]["type"]
            probe_value = probe_parameters[change.new]["target"]["value"]
            if probe_type == "Sequence":
                #print(f"This probe will label the amino acid sequence: {probe_value}")
                new_message = "This probe will label the amino acid sequence: " + str(probe_value)
            elif probe_type == "Atom_residue":
                new_message = "This probe will label all available residue(s): " + str(probe_value['residues'])
            elif probe_type == "Primary":
                new_message = "This is a secondary probe that labels this probe: " + str(probe_value)
            labels_gui["label_message"].value = new_message

        def activate_custom(b):
            #probe_widgets = dict()
            #plt.clf()
            #clear_output()
            for key, val in probe_widgets.items():
                probe_widgets[key] = False
            probe_widgets["label_3"] = True
            probe_widgets["mock_label_dropdown"] = True
            probe_widgets["label_4"] = True
            probe_widgets["as_linker"] = True
            probe_widgets["mock_type"] = True
            probe_widgets["mock_type_options1"] = True
            probe_widgets["mock_type_options2"] = True
            # widgets for advace setting
            #probe_widgets["label_5"] = True
            #probe_widgets["Suggestion"] = True
            #probe_widgets["mock_value"] = True
            #probe_widgets["radnom_value"] = True
            
            #probe_widgets["mock_fluo_dropdown"] = True
            probe_widgets["mock_Labelling_efficiency"] = True
            probe_widgets["mock_distance_to_epitope"] = True
            probe_widgets["wobble"] = True
            probe_widgets["wobble_theta"] = True
            probe_widgets["Add_mock"] = True
            probe_widgets["Label"] = True
            self._update_widgets(labels_gui, probe_widgets)
        
        def random_sequence(b):
            protein_name, _1, site, sequence = self.my_experiment.structure.get_peptide_motif()
            hint_message = "Our suggestion: C-terminal of " + protein_name + " (" + sequence + ")"
            labels_gui["Suggestion"].value = hint_message
            labels_gui["mock_value"].value = sequence

        def type_dropdown_change(change):
            labels_gui["mock_type_options1"].options = options_per_type1[change.new] 
            labels_gui["mock_type_options1"].value = options_per_type1[change.new][0]
            labels_gui["mock_type_options2"].options = options_per_type2[change.new]
            labels_gui["mock_type_options2"].value = options_per_type2[change.new][0]

        labels_gui.add_label("Use the dropdown menu to select an existing probe for your structure.")
        # initial buttons
        labels_gui.add_button("Demos", description="VLab4Mic probes")
        labels_gui.add_button("Customise", description="Customise your probe")
        # DEMOS
        labels_gui.add_label("Structure specific labels")
        if self.my_experiment.structure is not None:
            labels_gui.add_dropdown("label_dropdown", options=vlab_probes, value=vlab_probes[1])
            labels_gui._widgets["label_message"] = widgets.HTML(value = "")
            labels_gui["label_dropdown"].observe(update_label_message, names='value')
            labels_gui["label_dropdown"].value = vlab_probes[0]
            labels_gui.add_dropdown("fluo_dropdown", options=fluorophores_list)
            labels_gui.add_float_slider(
                "Labelling_efficiency",
                value=1,
                min=0,
                max=1,
                step=0.01,
                description="Labelling efficiency",
            )
            labels_gui.add_button(
                "Add", description="Add probe", disabled=(not vlab_probes)
            )
            labels_gui["Add"].on_click(build_label)
        # CUSTOM
        labels_gui.add_label("Customise your probe. First, select a probe model from the dropdown menu.")
        labels_gui.add_dropdown("mock_label_dropdown", options=mock_labels)
        labels_gui.add_label("Define the type of target for this probe.")
        options_dictionary = dict(
            Protein="Sequence",
            Residue="Atom_residue",
            Primary_Probe="Primary"
        )
        labels_gui.add_dropdown("mock_type", options=list(options_dictionary.keys()),description="I want this probe to target a: ")
        list_of_proteins = self.my_experiment.structure.list_protein_names()
        list_of_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        options_per_type1= dict(
            Protein=list_of_proteins,
            Residue=list_of_residues,
            Primary_Probe=[None,]
        )
        options_per_type2= dict(
            Protein=["cterminal", "nterminal"],
            Residue=["CA"],
            Primary_Probe=[None,]
        )

        labels_gui.add_dropdown("mock_type_options1", options=options_per_type1[labels_gui["mock_type"].value],description="Which one: ")
        labels_gui.add_dropdown("mock_type_options2", options=options_per_type2[labels_gui["mock_type"].value],description="Where: ")
        labels_gui["mock_type"].observe(type_dropdown_change, names='value')
        labels_gui.add_label("Type in a sequence motif, Residue names or name of a primary probe")
        protein_name = None
        sequence = None
        protein_name, _1, site, sequence = self.my_experiment.structure.get_peptide_motif()
        hint_message = "Our suggestion: C-terminal of " + protein_name + " (" + sequence + ")"
        labels_gui._widgets["Suggestion"] = widgets.HTML(value = hint_message)
        mock_value_suggestions = dict(
            Protein=sequence,
            Residue="LYS.CA",
            Primary_Probe=None
        )
        labels_gui.add_textarea("mock_value", description="At this position:", value=sequence)
        labels_gui.add_button("radnom_value", description="Give me another suggestion")
        labels_gui["radnom_value"].on_click(random_sequence)
        labels_gui.add_dropdown("mock_fluo_dropdown", options=fluorophores_list)
        labels_gui.add_float_slider(
            "mock_Labelling_efficiency",
            value=1,
            min=0,
            max=1,
            step=0.01,
            description="Labelling efficiency",
        )
        labels_gui.add_float_slider(
            "mock_distance_to_epitope",
            value=0,
            min=0,
            max=500,
            step=1,
            description="Distance from probe to epitope (A)",
        )
        labels_gui.add_label("Activate this option if you intent to use a secondary that recognises the current probe")
        labels_gui.add_checkbox("as_linker", description = "Model as primary with epitope for secondary probe",
                        value = False)
        labels_gui.add_checkbox("wobble", description = "Enable wobble",
                        value = False)
        labels_gui.add_float_slider(
            "wobble_theta",
            value=10,
            min=0,
            max=45,
            step=1,
            description="Wobble cone range (degrees)",
        )
        labels_gui.add_button("Add_mock", description="Add custom probe")
        labels_gui["Add_mock"].on_click(build_mock_label)


        labels_gui.add_button("Clear", description="Clear Labels")
        labels_gui.add_button("Show", description="Display current labels")
        labels_gui.add_label(
            "After adding your probes, click on 'Done' to continue"
        )
        labels_gui.add_button("Label", description="Done", disabled = True)
        labels_gui["Demos"].on_click(activate_demos)
        labels_gui["Customise"].on_click(activate_custom)
        labels_gui["Clear"].on_click(clear)
        labels_gui["Show"].on_click(show)
        labels_gui["Label"].on_click(label_struct)
        probe_widgets = {}
        for wgt in labels_gui._widgets.keys():
            probe_widgets[wgt] = False
            labels_gui._widgets[wgt].layout = widgets.Layout(width="50%", display="None")
        labels_gui.show()
        probe_widgets["Demos"] = True
        probe_widgets["Customise"] = True
        if self.my_experiment.structure is None:
            print("No structure has been loaded")
        else:
            self.my_experiment.remove_probes()
            probe_widgets["Demos"] = True
            probe_widgets["Customise"] = True
            self._update_widgets(labels_gui, probe_widgets)
            probe_widgets["Demos"] = False
            probe_widgets["Customise"] = False
            #display(labels_gui["Demos"], labels_gui["Customise"])

    def refine_model_gui(self):
        width = "50%"
        style = {"description_width": "initial"}
        structural_model_gui = easy_gui_jupyter.EasyGUI("StructuralModel",width=width)
        def show_model(b):
            #clear_output()
            # plot
            emitter_plotsize = structural_model_gui["emitterplotsize"].value
            source_size = structural_model_gui["sourceplotsize"].value
            plot_size = structural_model_gui["plot_size"].value
            structural_model_gui.show()
            if self.my_experiment.particle:
                particle = self.my_experiment.particle
                fig, axs = plt.subplots(1, 3, subplot_kw={"projection": "3d"}, figsize = [plot_size, plot_size/2])
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[0, 0, 0],
                    axis_object=axs[0],
                    emitter_plotsize=emitter_plotsize,
                )
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[30, 0, 0],
                    axis_object=axs[1],
                    emitter_plotsize=emitter_plotsize,
                )
                particle.gen_axis_plot(
                    with_sources=structural_model_gui["WTarget"].value,
                    source_plotsize=source_size,
                    axesoff=structural_model_gui["Axes"].value,
                    view_init=[90, 0, 0],
                    axis_object=axs[2],
                    emitter_plotsize=emitter_plotsize,
                )
                plt.subplots_adjust(wspace=0.5)
                plt.close()
                structural_model_gui["image_output"].clear_output()
                with structural_model_gui["image_output"]:
                    display(fig)
            else:
                print(
                    "You have not created a labelled structure. "
                    "Make sure you select 'Label structure' button on previous cell"
                )
        def activate_visualisation_options(b):
            # show hidden widgets
            structural_model_gui["show_preview"].layout.display = "None"
            structural_model_gui["label_1"].layout.display = "block"
            structural_model_gui["plot_size"].layout.display = 'inline-flex'
            structural_model_gui["emitterplotsize"].layout.display = 'inline-flex'
            structural_model_gui["sourceplotsize"].layout.display = 'inline-flex'
            structural_model_gui["WTarget"].layout.display = "block"
            structural_model_gui["Axes"].layout.display = "block"
            structural_model_gui["update_plot"].layout.display = "block"
            structural_model_gui["Relabel"].layout.display = "block"
            structural_model_gui["label_2"].layout.display = "block"
            structural_model_gui["enable_defects"].layout.display = "block"
            show_model(b)


        def activate_defects(b):
            structural_model_gui._widgets["eps1"].layout.display = "block"
            structural_model_gui._widgets["xmer_neigh_distance"].layout.display = "block"
            structural_model_gui._widgets["Defect"].layout.display = "block"
            structural_model_gui._widgets["use_defects"].layout.display = "block"
            structural_model_gui["enable_defects"].disabled = True


        def add_defects(b):
            if self.nlabels == 1:
                self.my_experiment.particle.add_defects(
                    eps1=structural_model_gui["eps1"].value,
                    xmer_neigh_distance=structural_model_gui[
                        "xmer_neigh_distance"
                    ].value,
                    deg_dissasembly=structural_model_gui["Defect"].value,
                )
                message = "Defects added"
            else:
                message = (
                    "Defect modelling is currently unsupported for more than one label"
                )
            structural_model_gui.save_settings()
            show_model(b)
            print(message)

        def update_plot(b):
            show_model(b)

        def relabel(b):
            self.my_experiment.particle.generate_instance()
            show_model(b)

        structural_model_gui.add_button("show_preview", description="Preview your model")
        structural_model_gui.add_label("Visualisation parameters")
        structural_model_gui["label_1"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_int_slider(
            "plot_size",
            value=10,
            min=5,
            max=24,
            step=1,
            description="Figure size",
        )
        structural_model_gui["plot_size"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_float_slider(
            "emitterplotsize",
            value=24,
            min=0,
            max=50,
            step=1,
            description="Emitter size",
        )
        structural_model_gui["emitterplotsize"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui["emitterplotsize"].style = style
        structural_model_gui.add_float_slider(
            "sourceplotsize", value=1, min=0, max=50, step=1, description="Target size"
        )
        structural_model_gui["sourceplotsize"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui["sourceplotsize"].style = style
        structural_model_gui.add_checkbox(
            "WTarget", description="With target site", value=True
        )
        structural_model_gui["WTarget"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_checkbox("Axes", description="Hide Axes", value=True)
        structural_model_gui["Axes"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_button("update_plot", description="Update figure")
        structural_model_gui["update_plot"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_button(
            "Relabel", description="Relabel particle and update figure"
        )
        structural_model_gui["Relabel"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_label("Model defects parameters (optional):")
        structural_model_gui["label_2"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui._widgets["eps1"] = widgets.BoundedIntText(
            value=300,
            min=0,
            max=100000,
            description="Short distance cluster",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui._widgets["eps1"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui._widgets["xmer_neigh_distance"] = widgets.BoundedIntText(
            value=600,
            min=0,
            max=100000,
            description="Long distance cluster",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui._widgets["xmer_neigh_distance"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui._widgets["Defect"] = widgets.BoundedFloatText(
            value=0.5,
            min=0,
            max=1,
            description="percentage of defect",
            layout=structural_model_gui._layout,
            style=structural_model_gui._style,
            remember_value=True,
        )
        structural_model_gui._widgets["Defect"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui.add_button(
            "use_defects", description="Apply defects and update figure"
        )
        structural_model_gui._widgets["use_defects"].layout = widgets.Layout(width=width, display = "None")
        #
        structural_model_gui.add_button(
            "enable_defects", description="Enable particle defects"
        )
        structural_model_gui["enable_defects"].layout = widgets.Layout(width=width, display = "None")
        structural_model_gui["enable_defects"].on_click(activate_defects)
        structural_model_gui["show_preview"].on_click(activate_visualisation_options)
        structural_model_gui["update_plot"].on_click(update_plot)
        structural_model_gui["use_defects"].on_click(add_defects)
        structural_model_gui["Relabel"].on_click(relabel)
        structural_model_gui._widgets["image_output"] = widgets.Output()
        if self.my_experiment.particle:
            structural_model_gui.show()
        else:
            print("No particle has been created")

    def create_field(self):
        field_gui = easy_gui_jupyter.EasyGUI("field")
        width = "50%"
        style = {"description_width": "initial"}
        def createmin(b):
            random_pl = field_gui["random"].value
            nparticles = field_gui["nparticles"].value
            min_distance = field_gui["mindist"].value
            # set parameters for virtual samle
            #self.virtualsample_params
            # then create the field
            self.my_experiment._build_coordinate_field(
                keep=True,
                nparticles=nparticles,
                use_self_particle=False,
                random_placing=random_pl,
                minimal_distance=min_distance
            )
            field_gui["Use_Upload"].disabled = True
            field_gui["show"].layout.display = "block"
            #field_gui["show"].disabled = False

        def createmin_particles(b):
            random_pl = field_gui["random"].value
            nparticles = field_gui["nparticles"].value
            if field_gui["distance_from_particle"].value:
                min_distance = None
            else:
                min_distance = field_gui["mindist"].value
            random_orientations = field_gui["random_orientations"].value
            self.my_experiment._build_coordinate_field(
                keep=True,
                nparticles=nparticles,
                random_placing=random_pl,
                minimal_distance=min_distance,
                random_orientations = random_orientations
                )
            field_gui["Use_Upload"].disabled = True
            field_gui["show"].layout.display = "block"

        def upload(b):
            filepath = field_gui["File"].selected
            img = tif.imread(filepath)
            xyz_relative, image_physical_size = field.gen_positions_from_image(
                img,
                mode=field_gui["Mode"].value,
                sigma=field_gui["sigma"].value,
                background=field_gui["background"].value,
                threshold=field_gui["threshold"].value,
                pixelsize=field_gui["pixelsize"].value,
                min_distance=field_gui["min_distance"].value)
            coordinates_field = field.create_min_field(relative_positions=xyz_relative)
            if particle is not None:
                coordinates_field.create_molecules_from_InstanceObject(particle)
            coordinates_field.construct_static_field()
            exported_field = coordinates_field.export_field()
            #display(field_gui["show"])
            field_gui["minimal"].disabled = True
            field_gui["minimal_particles"].disabled = True
            field_gui["show"].layout.display = "block"
        
        def useimage(b):
            field_widgets["File"] = True
            field_widgets["pixelsize"] = True
            field_widgets["background"] = True
            field_widgets["min_distance"] = True
            field_widgets["sigma"] = True
            field_widgets["threshold"] = True
            field_widgets["Mode"] = True
            field_widgets["Upload"] = True
            field_widgets["show"] = True
            self._update_widgets(field_gui, field_widgets)

        def showfield(b):
            #plt.clf()
            plt.close()
            field_widgets["image_output"] = True
            self._update_widgets(field_gui, field_widgets)
            field_gui["show"].disabled = True
            field_gui["image_output"].clear_output()
            with field_gui["image_output"]:
                display(
                    self.my_experiment.coordinate_field.show_field(
                        view_init=[90, 0, 0], 
                        initial_pos=False, 
                        return_fig=True)
                )
            plt.close()

        field_gui.add_file_upload(
                "File", description="Select from file", accept="*.tif", save_settings=False
            )
        field_gui._widgets["pixelsize"] = widgets.BoundedIntText(
            value=15,
            description="Image pixelsize",
            layout=field_gui._layout,
            style=field_gui._style,
            remember_value=True,
        )
        field_gui._widgets["background"] = widgets.BoundedIntText(
            value=0,
            max=100000,
            description="Background",
            layout=field_gui._layout,
            style=field_gui._style,
            remember_value=True,
        )
        field_gui._widgets["min_distance"] = widgets.BoundedIntText(
            value=1,
            description="min_distance",
            layout=field_gui._layout,
            style=field_gui._style,
            remember_value=True,
        )
        field_gui._widgets["sigma"] = widgets.BoundedIntText(
            value=2,
            description="sigma",
            layout=field_gui._layout,
            style=field_gui._style,
            remember_value=True,
        )
        field_gui._widgets["threshold"] = widgets.BoundedIntText(
            value=2,
            description="threshold",
            layout=field_gui._layout,
            style=field_gui._style,
            remember_value=True,
        )
        field_gui.add_dropdown("Mode", options=["mask", "localmaxima"]) 
        field_gui.add_checkbox("random", value=True, description="Randomise positions (enforced when there is more than one particle)")
        field_gui.add_checkbox("random_orientations", value=True, description="Randomise orientation")
        field_gui.add_checkbox("distance_from_particle", value=True, description="Use minimal distance of particle")
        
        #field_gui["show"].layout = widgets.Layout(width=width, display = "None")
        field_gui.add_int_slider(
                "nparticles",
                value=1,
                min=1,
                max=20,
                step=1,
                description="Number of Particles in field",
                disabled=False,
            )
        field_gui._widgets["mindist"] = widgets.BoundedIntText(
            value=100,
            min=1,
            max=1000,
            description="Minimal distance from particle center (nanometers)",
            layout=field_gui._layout,
            style=field_gui._style
        )
        field_gui.add_button("minimal_particles", description="Create field with particle")
        field_gui.add_button("Use_Upload", description="Use image for positioning")
        field_gui.add_button("Upload", description="Load image")
        field_gui.add_button("minimal", description="Create field")
        field_gui.add_button("show", description="Show field", disabled=False)
        field_gui["minimal"].on_click(createmin)
        field_gui["show"].on_click(showfield)
        field_gui["minimal_particles"].on_click(createmin_particles)
        field_gui["Use_Upload"].on_click(useimage)
        field_gui["Upload"].on_click(upload)
        field_gui._widgets["image_output"] = widgets.Output()
        field_widgets = {}
        for wgt in field_gui._widgets.keys():
            field_widgets[wgt] = False
            field_gui._widgets[wgt].layout = widgets.Layout(width="50%", display="None")
        field_gui.show()
        if self.my_experiment.generators_status("particle"):
            print("Using particle model")
            field_widgets["nparticles"] = True
            field_widgets["random"]= True
            field_widgets["random_orientations"]= True
            field_widgets["distance_from_particle"]= True
            field_widgets["mindist"]= True
            field_widgets["minimal_particles"]= True
            field_widgets["Use_Upload"]= True
            #field_widgets["show"]= True
        else:
            print("No particle available")
            field_widgets["nparticles"] = True
            field_widgets["random"] = True
            field_widgets["mindist"] = True
            field_widgets["minimal"] = True
            field_widgets["Use_Upload"] = True
            #field_widgets["show"] = True
        self._update_widgets(field_gui, field_widgets)
        field_widgets["Use_Upload"] = False
        field_widgets["minimal_particles"] = False
        field_widgets["minimal"] = False
        


    def set_image_modalities(self, mode = "default"):
        imager_created = False
        modalities_options = []
        if mode == "default":
            modalities_list = self.modalities_default
        else:
            modalities_list = self.my_experiment.local_modalities_names
        modalities_options = copy.copy(modalities_list)
        modalities_options.append("All")
        imaging_gui = easy_gui_jupyter.EasyGUI("imaging", width="70%")
        selected_mods = []
        modality_info = {}
        for mod in modalities_list:
            mod_info = data_format.configuration_format.compile_modality_parameters(
                mod, self.config_directories["base"]
            )
            modality_info[mod] = mod_info
        with io.capture_output() as captured:
            temp_imager, tmp_modality_parameters = create_imaging_system(
                modalities_id_list=modalities_list,
                config_dir=self.config_directories["base"],
            )

        def add_mod(b):
            if "All" == imaging_gui["modalities_dropdown"].value:
                for mod_names in modalities_list:
                    self.my_experiment.add_modality(
                        modality_name=mod_names,
                        save=True
                    )
                imaging_gui["Add"].disabled=True
            else:
                self.my_experiment.add_modality(
                    modality_name=imaging_gui["modalities_dropdown"].value,
                    save=True
                )
            
            imaging_gui["label_4"].value = str(list(self.my_experiment.imaging_modalities.keys()))

        def clear(b):
            #selected_mods.clear()
            mods = list(self.my_experiment.imaging_modalities.keys())
            for mod_name in mods:
                self.my_experiment.update_modality(mod_name, remove=True)
            imaging_gui["label_4"].value = "No modality selected"

        def create_imager(b):
            if self.my_experiment.exported_coordinate_field is not None:
                self.my_experiment.build(modules=["imager",])
                print("Imaging system created")
                imaging_gui["Show"].disabled = False
                imaging_gui["modalities_dropdown"].disabled = True
                imaging_gui["Create"].disabled = True
                imaging_gui["Clear"].disabled = True
                imaging_gui["Add"].disabled = True
            else:
                print("No field info")

        def preview(b):
            def get_info(imaging_gui):
                def preview_info(message, Modality):
                    if Modality != "All":
                        pixelsize = modality_info[Modality]["detector"]["pixelsize"]
                        pixelsize_nm = pixelsize * 1000
                        psf_sd = np.array(modality_info[Modality]["psf_params"]["std_devs"])
                        psf_voxel = np.array(
                            modality_info[Modality]["psf_params"]["voxelsize"]
                        )
                        psf_sd_metric = np.multiply(psf_voxel, psf_sd)
                        print(f"Detector pixelsize (nm): {pixelsize_nm}")
                        print(f"PSF sd (nm): {psf_sd_metric}")
                        print(f"PSF preview (on a 1x1 µm field of view): ")
                        # show PSF
                        fig, axs = plt.subplots()
                        modality_preview = temp_imager.modalities[Modality]["psf"][
                            "psf_stack"
                        ]
                        psf_shapes = modality_preview.shape
                        stack_max = np.max(modality_preview)
                        axs.imshow(
                            modality_preview[:, :, int(psf_shapes[2] / 2)],
                            cmap="gray",
                            interpolation="none",
                            vmin=0,
                            vmax=stack_max,
                        )
                        axs.set_xticks([])
                        axs.set_yticks([])

                widgets.interact(
                    preview_info, message=imaging_gui["label_1"], Modality=imaging_gui["modalities_dropdown"]
                )

            get_info(imaging_gui)

        def showfov(b):
            self.my_experiment.imager.show_field()

        imaging_gui.add_label("Use the dropdown menu to get a preview of the selected imaging modality")
        imaging_gui.add_dropdown("modalities_dropdown", options=modalities_options)
        imaging_gui.add_label("Click on Add modality to include the selected modality to the virtual microscope")
        imaging_gui.add_button("Add", description="Add modality")
        imaging_gui.add_label("Selected modalities:")
        imaging_gui.add_label("No modality selected")
        imaging_gui.add_button("Clear", description="Reset selected modalities")
        if self.my_experiment.generators_status("exported_coordinate_field") is None:
            disable_button = True
        else:
            disable_button = False
        imaging_gui.add_label("Click 'Done!' after selecting your modalities ")
        imaging_gui.add_button(
            "Create", description="Done!", disabled=disable_button
        )
        imaging_gui.add_button("Show", description="Show field of view (yellow)", disabled=True)
        imaging_gui.add_int_slider(
            "PSF_nslice", min=0, max=400, continuous_update=False
        )
        imaging_gui["Add"].on_click(add_mod)
        imaging_gui["Clear"].on_click(clear)
        imaging_gui["Create"].on_click(create_imager)
        imaging_gui["Show"].on_click(showfov)
        preview(True)
        display(
            imaging_gui["label_2"],
            imaging_gui["Add"],
            imaging_gui["label_3"],
            imaging_gui["label_4"],
            imaging_gui["Clear"],
            imaging_gui["label_5"],
            imaging_gui["Create"],
            imaging_gui["Show"],
        )


    def set_acq_params(self):
        acquisition_gui = easy_gui_jupyter.EasyGUI("acquisition_params")
        imager_channels = []
        anymod = list(self.my_experiment.imager.modalities.keys())[0]
        for chann in self.my_experiment.imager.modalities[anymod]["filters"].keys():
            print(chann)
            imager_channels.append(chann)
        nchannels = len(imager_channels)

        def set_params(b):
            mod_id = acquisition_gui["modalities_dropdown"].value
            exp_time = acquisition_gui["Exposure"].value
            noise = acquisition_gui["Noise"].value
            nframes = acquisition_gui["Frames"].value
            if acquisition_gui["Channels"].value:
                channels = []
                for chann in self.my_experiment.imager.modalities[mod_id]["filters"].keys():
                    channels.append(chann)
                print(f"using all channels: {channels}")
            else:
                channels = [
                    "ch0",
                ]
            self.my_experiment.set_modality_acq(
                modality_name=mod_id,
                exp_time=exp_time,
                noise=noise,
                save=True,
                nframes=nframes,
                channels=channels
            )
            acq_per_modalities[mod_id] = "Custom"
            output_str = '<br>'.join(f"{name}: {val}" for name, val in self.my_experiment.selected_mods.items())
            acquisition_gui["message"].value = output_str
            acquisition_gui.save_settings()

        def preview_mod(b):
            preview_image = None

            def get_preview(imaging_system, acq_gui):

                def preview_exposure(message, Modality, Exposure, Noise):
                    fig = plt.figure()
                    grid = axes_grid1.AxesGrid(
                        fig,
                        111,
                        nrows_ncols=(1, nchannels),
                        axes_pad=1,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="10%",
                        cbar_pad="20%",
                    )
                    i = 0
                    for single_channel in imager_channels:
                        single_mod_acq_params = dict(
                            exp_time=Exposure,
                            noise=Noise,
                            save=False,
                            nframes=1,
                            channel=single_channel,
                        )
                        with io.capture_output() as captured:
                            timeseries, calibration_beads = imaging_system.generate_imaging(
                                modality=Modality, **single_mod_acq_params
                            )
                            min_val = np.min(timeseries[0])
                            max_val = np.max(timeseries[0])
                        preview_image = grid[i].imshow(
                            timeseries[0],
                            cmap="gray",
                            interpolation="none",
                            vmin=min_val,
                            vmax=max_val,
                        )
                        grid[i].set_xticks([])
                        grid[i].set_yticks([])
                        grid[i].set_title("preview channel:" + single_channel)
                        grid.cbar_axes[i].colorbar(preview_image)
                        i = i + 1

                def preview_parameters(Settings):
                    pass

                widgets.interact(
                    preview_exposure,
                    message = acq_gui["label_1"],
                    Modality=acq_gui["modalities_dropdown"],
                    Exposure=acq_gui["Exposure"],
                    Noise=acq_gui["Noise"]
                )
                widgets.interact(
                    preview_parameters,
                    Settings = acq_gui["message"]
                )


            get_preview(self.my_experiment.imager, acquisition_gui)

        def clear(b):
            print("Acquisition parameters cleared")
            self.my_experiment.reset_to_defaults(module="acquisitions",save=True)
            output_str = '<br>'.join(f"{name}: {val}" for name, val in self.my_experiment.selected_mods.items())
            acquisition_gui["message"].value = output_str
            acquisition_gui.save_settings()

        acquisition_gui.add_label("Set acquisition parameters")
        acquisition_gui.add_label("Acquisition parameters per modality:")
        acq_per_modalities = dict()
        for mods_selected in self.my_experiment.selected_mods.keys():
            acq_per_modalities[mods_selected] = "Default"
        output_str = ""
        output_str = '<br>'.join(f"{name}: {val}" for name, val in self.my_experiment.selected_mods.items())
        acquisition_gui._widgets["message"] = widgets.HTML(value = output_str)
        selected_mods = list(self.my_experiment.imaging_modalities.keys())
        acquisition_gui.add_dropdown("modalities_dropdown", options=selected_mods)
        acquisition_gui.add_checkbox("Noise", description="Use Noise", value=True)
        acquisition_gui.add_checkbox("Channels", description="Use all channels", value=True)
        ## bounded int Text
        acquisition_gui._widgets["Frames"] = widgets.BoundedIntText(
            value=1,
            min=1,
            max=100000,
            description="Frames (not used for preview)",
            layout=acquisition_gui._layout,
            style=acquisition_gui._style,
            remember_value=True,
        )
        acquisition_gui._widgets["Exposure"] = widgets.BoundedFloatText(
            value=0.001,
            min=0.000000,
            step=0.0001,
            description="exposure (sec)",
            layout=acquisition_gui._layout,
            style=acquisition_gui._style,
            remember_value=True,
        )
        acquisition_gui.add_button("Set", description="Update acquisition parameters")
        acquisition_gui.add_button("Preview", description="Preview (Expermiental feature)")
        acquisition_gui.add_button("Clear", description="Reset params")
        acquisition_gui["Preview"].on_click(preview_mod)
        acquisition_gui["Set"].on_click(set_params)
        acquisition_gui["Clear"].on_click(clear)
        display(acquisition_gui["Frames"])
        preview_mod(True)
        display(acquisition_gui["Set"], acquisition_gui["Clear"])


    def acquire_images(self):
        experiment_gui = easy_gui_jupyter.EasyGUI("experiment")

        def run_simulation(b):
            experiment_gui["Acquire"].disabled = True
            sav_dir = experiment_gui["saving_directory"].value
            if sav_dir is not None:
                self.my_experiment.output_directory = sav_dir
                save=True
            exp_name = experiment_gui["experiment_name"].value
            self.my_experiment.run_simulation(
                name=exp_name,
                save=save
            )
            experiment_gui.save_settings()

        experiment_gui.add_label("Set experiment name")
        experiment_gui.add_textarea(
            "experiment_name", value="Exp_name", remember_value=True
        )
        experiment_gui.add_label("Set saving directory")
        experiment_gui.add_textarea("saving_directory", remember_value=True)
        experiment_gui.add_button("Acquire", description="Run Simulation")
        experiment_gui["Acquire"].on_click(run_simulation)
        experiment_gui.show()