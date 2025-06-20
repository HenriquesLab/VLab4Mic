from ezinput import EZInput
import matplotlib.pyplot as plt
from IPython.display import display, clear_output


def ui_show_structure(experiment):
    gui = EZInput(title="Structure")

    def show_structure(widget_elements):
        widget_elements["preview_structure"].clear_output()
        total = experiment.structure.num_assembly_atoms
        widget_elements["n_atoms"].disabled = False
        atoms_number = widget_elements["n_atoms"].value
        if total > atoms_number:
            fraction = atoms_number/total
        else:
            fraction = 1.0
        with widget_elements["preview_structure"]:
            display(
                experiment.structure.show_assembly_atoms(assembly_fraction=fraction)
            )
            plt.close()
    
    gui.add_callback(
        "button",
        show_structure,
        gui.elements,
        description="Show structure",
    )

    def update_plot(value):
        gui["preview_structure"].clear_output()
        total = experiment.structure.num_assembly_atoms
        if total > value.new:
            fraction = value.new/total
        else:
            fraction = 1.0
        with gui["preview_structure"]:
            display(
                experiment.structure.show_assembly_atoms(assembly_fraction=fraction)
            )
            plt.close()


    gui.add_int_slider("n_atoms", description="Atoms to use", min=1, max=10000, step=1, value = 1000, on_change=update_plot, continuous_update=False, disabled=True)

    gui.add_output("preview_structure")
    gui["preview_structure"].clear_output()
    return gui