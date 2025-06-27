from ezinput import EZInput
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from supramolsim.utils.visualisation.matplotlib_plots import slider_normalised
import ipywidgets as widgets
import mpl_toolkits.axes_grid1 as axes_grid1
import io
import numpy as np
from IPython.utils import io

def update_widgets_visibility(ezwidget, visibility_dictionary):
    for widgetname in visibility_dictionary.keys():
        if visibility_dictionary[widgetname]:
            ezwidget[widgetname].layout.display = "inline-flex"
        else:
            ezwidget[widgetname].layout.display = "inline-flex"



def ui_show_structure(experiment):
    gui = EZInput(title="Structure")
    def show_structure(widget_elements):
        if experiment.objects_created["structure"]:
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
        else:
            widget_elements["preview_structure"].clear_output()
            with widget_elements["preview_structure"]:
                print("Structure not created yet, please create it first.")
    
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

    def show_labelled_structure(change):
        gui["preview_labelled_structure"].clear_output()
        with gui["preview_labelled_structure"]:
            if not experiment.objects_created["particle"]:
                display("Particle not created yet, please create it first.")
            else:
                gui["emitter_plotsize"].disabled = False
                gui["source_plotsize"].disabled = False
                gui["hview"].disabled = False
                gui["vview"].disabled = False
                display(show_particle(
                    emitter_plotsize=gui["emitter_plotsize"].value,
                    source_plotsize=gui["source_plotsize"].value,
                    hview=gui["hview"].value,
                    vview=gui["vview"].value
            ))
    
    gui.add_int_slider(
        "emitter_plotsize",
        description="Emitter size",
        min=0,
        max=30,
        step=1,
        value=1,
        continuous_update=False,
        on_change=show_labelled_structure,
        disabled=True
    )
    gui.add_int_slider(
        "source_plotsize",
        description="Epitope size",
        min=0,
        max=30,
        step=1,
        value=1,
        continuous_update=False,
        on_change=show_labelled_structure,
        disabled=True
    )
    gui.add_int_slider(
        "hview",
        description="Horizontal view",
        min=-90,
        max=90,
        step=1,
        value=0,
        continuous_update=False,
        on_change=show_labelled_structure,
        disabled=True
    )
    gui.add_int_slider(
        "vview",
        description="Vertical view",
        min=-90,
        max=90,
        step=1,
        value=0,
        continuous_update=False,
        on_change=show_labelled_structure,
        disabled=True
    )
    gui.add_button(
        "show_labelled_structure",
        description="Show labelled structure",
    )
    gui["show_labelled_structure"].on_click(show_labelled_structure)
    gui.add_output("preview_labelled_structure")
    return gui


def ui_show_virtual_sample(experiment):
    gui = EZInput(title="Virtual Sample")

    def update_plot(change):
        gui["preview_virtual_sample"].clear_output()
        hview = gui["horizontal_view"].value
        vview = gui["vertical_view"].value
        with gui["preview_virtual_sample"]:
            display(experiment.coordinate_field.show_field(
                view_init=[vview, hview, 0],
                return_fig=True))
            plt.close()

    gui.add_int_slider(
        "horizontal_view",
        description="Rotation angle (degrees)",
        min=-90,
        max=90,
        step=1,
        value=0,
        continuous_update=False,
        on_change=update_plot,
    )
    gui.add_int_slider(
        "vertical_view",
        description="Tilt angle (degrees)",
        min=-90,
        max=90,
        step=1,
        value=90,
        continuous_update=False,
        on_change=update_plot,
    )

    gui.add_output("preview_virtual_sample")
    gui["preview_virtual_sample"].clear_output()
    update_plot(True)
    return gui


def ui_show_modality(experiment):
    gui = EZInput(title="Modality")
    xy_zoom_in = 0.5
    def update_plot(change):
        mod_name = gui["modality"].value
        psf_stack = experiment.imager.get_modality_psf_stack(mod_name)
        psf_shape = psf_stack.shape
        half_xy = int(psf_shape[0] / 2)
        half_z = int(psf_shape[2] / 2)
        psf_stack = psf_stack[
            half_xy - int(half_xy * xy_zoom_in) : half_xy + int(half_xy * xy_zoom_in),
            half_xy - int(half_xy * xy_zoom_in) : half_xy + int(half_xy * xy_zoom_in),
            :]
        dimension_plane = gui["dimension_slice"].value
        if dimension_plane == "YZ plane":
            dimension = 0
        elif dimension_plane == "XZ plane":
            dimension = 1
        elif dimension_plane == "XY plane":
            dimension = 2
        gui["preview_modality"].clear_output()
        with gui["preview_modality"]:
            display(slider_normalised(
                psf_stack,
                dimension=dimension,
                cbar=False,))

    current_modalities = list(experiment.imaging_modalities.keys())
    if len(current_modalities) == 0:
        value_modalities = None
    else:
        value_modalities = current_modalities[0]
    gui.add_dropdown(
        "modality",
        description="Modality",
        options=current_modalities,
        value=value_modalities,
        on_change=update_plot,
    )
    gui.add_custom_widget(
        "dimension_slice",
        widgets.ToggleButtons,
        options=["YZ plane", "XZ plane", "XY plane"],
        value="XY plane",
        on_change=update_plot,
        style={"description_width": "initial"},
        description="Plane of view: ",
    )
    gui.add_output("preview_modality")
    gui["preview_modality"].clear_output()
    if value_modalities is not None:
        update_plot(True)
    return gui


def ui_set_acq_params(experiment):
    acquisition_gui = EZInput(title="acquisition_params")
    imager_channels = []
    anymod = list(experiment.imager.modalities.keys())[0]
    for chann in experiment.imager.modalities[anymod]["filters"].keys():
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
            for chann in experiment.imager.modalities[mod_id][
                "filters"
            ].keys():
                channels.append(chann)
            print(f"using all channels: {channels}")
        else:
            channels = [
                "ch0",
            ]
        experiment.set_modality_acq(
            modality_name=mod_id,
            exp_time=exp_time,
            noise=noise,
            save=True,
            nframes=nframes,
            channels=channels,
        )

    def preview_mod(b):
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
                        timeseries, calibration_beads = (
                            imaging_system.generate_imaging(
                                modality=Modality, **single_mod_acq_params
                            )
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
                    # grid[i].set_visible(False)
                return fig

            figure = preview_exposure(
                message=acq_gui["label_1"].value,
                Modality=acq_gui["modalities_dropdown"].value,
                Exposure=acq_gui["Exposure"].value,
                Noise=acq_gui["Noise"].value,
            )
            plt.close()
            acquisition_gui["image_output"].clear_output()
            with acquisition_gui["image_output"]:
                display(figure)

        get_preview(experiment.imager, acquisition_gui)

    def clear(b):
        print("Acquisition parameters cleared")
        experiment.reset_to_defaults(module="acquisitions", save=True)

    def preview_params_chage(change):
        preview_mod(True)

    acquisition_gui.add_label("Set acquisition parameters")
    selected_mods = list(experiment.imaging_modalities.keys())
    acquisition_gui.add_dropdown("modalities_dropdown", options=selected_mods)
    acquisition_gui.add_checkbox("Noise", description="Use Noise", value=True)
    acquisition_gui.add_checkbox(
        "Channels", description="Use all channels", value=True
    )
    ## bounded int Text
    acquisition_gui.add_bounded_int_text(
        "Frames",
        description="Frames (not used for preview)",
        vmin=1,
        vmax=100000,
        value=1,
        step=1,
    )
    acquisition_gui.add_bounded_float_text(
        "Exposure",
        description="Exposure (sec)",
        vmin=0.000000,
        vmax=10.0,
        step=0.0001,
        value=0.01,
    )   
    acquisition_gui["modalities_dropdown"].observe(
        preview_params_chage, names="value"
    )
    acquisition_gui["Noise"].observe(preview_params_chage, names="value")
    acquisition_gui["Exposure"].observe(preview_params_chage, names="value")
    acquisition_gui.add_button("Set", description="Update acquisition parameters")
    acquisition_gui.add_button("Clear", description="Reset params")
    acquisition_gui["Set"].on_click(set_params)
    acquisition_gui["Clear"].on_click(clear)
    acquisition_gui.elements["image_output"] = widgets.Output()
    acq_widgets = {}
    for wgt in acquisition_gui.elements.keys():
        acq_widgets[wgt] = False
        acquisition_gui.elements[wgt].layout = widgets.Layout(
            width="50%", display="None"
        )
    acquisition_gui.show()
    acq_widgets["Frames"] = True
    acq_widgets["Set"] = True
    acq_widgets["image_output"] = True
    acq_widgets["label_1"] = True
    acq_widgets["modalities_dropdown"] = True
    acq_widgets["Exposure"] = True
    acq_widgets["Noise"] = True
    acq_widgets["Clear"] = True
    update_widgets_visibility(acquisition_gui, acq_widgets)
    preview_mod(True)
    return acquisition_gui