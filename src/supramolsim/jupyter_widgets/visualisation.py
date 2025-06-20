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


def ui_show_labelled_structure(experiment):
    gui = EZInput(title="Labelled Structure")

    def show_particle(
                    emitter_plotsize = 1, 
                    source_plotsize = 1, 
                    hview=0,
                    vview=0):
        #with io.capture_output() as captured:   
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        experiment.particle.gen_axis_plot(
            axis_object=ax,
            with_sources=True, 
            axesoff=True,
            emitter_plotsize=emitter_plotsize,
            source_plotsize=source_plotsize,
            view_init=[vview, hview, 0]
            )
        plt.close()
        return fig


    def show_labelled_structure(values):
        gui["preview_labelled_structure"].clear_output()
        with gui["preview_labelled_structure"]:
            display(show_particle())
    
    
    gui.add_callback(
        "button",
        show_labelled_structure,
        gui.elements,
        description="Show labelled structure",
    )
    
    gui.add_output("preview_labelled_structure")
    
    return gui