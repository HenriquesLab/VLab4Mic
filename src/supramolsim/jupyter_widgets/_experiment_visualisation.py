"""
_experiment_visualisation
-------------------------

This module provides EZInput-based widget interfaces for visualising experiment components and results
in the virtual microscopy simulation workflow. It enables interactive visualisation of structures, labelled structures,
virtual samples, imaging modalities, acquisition parameters, and experiment results within Jupyter notebooks.

Functions
---------

- update_widgets_visibility(ezwidget, visibility_dictionary):
    Utility to show/hide widgets in an EZInput widget based on a visibility dictionary.

- ui_show_structure(experiment):
    Returns a widget for visualising the experiment structure.

- ui_show_labelled_structure(experiment):
    Returns a widget for visualising the labelled structure (with probes).

- ui_show_virtual_sample(experiment):
    Returns a widget for visualising the virtual sample.

- ui_show_modality(experiment):
    Returns a widget for visualising the point spread function (PSF) of selected imaging modalities.

- ui_set_acq_params(experiment):
    Returns a widget for setting and previewing acquisition parameters.

- ui_preview_results(experiment):
    Returns a widget for previewing experiment results.


Each function returns an EZInput-based widget for use in a Jupyter notebook.
"""

from ezinput import EZInput
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from supramolsim.utils.visualisation.matplotlib_plots import slider_normalised
import ipywidgets as widgets
import mpl_toolkits.axes_grid1 as axes_grid1
import io
import numpy as np
from IPython.utils import io
from ..utils.visualisation.matplotlib_plots import plot_projection


select_colour = "#4daf4ac7"
remove_colour = "#ff8000da"
update_colour = "#00bfffda"

select_icon = "fa-check"
add_icon = "fa-plus"
show_icon = "fa-eye"
remove_icon = "fa-minus"
loding_icon = "fa-spinner fa-spin"
update_icon = "fa-wrench" # create
toggle_icon = "fa-eye-slash"
upload_icon = "fa-upload"
reset_icon = "fa-undo"

def update_widgets_visibility(ezwidget, visibility_dictionary):
    """
    Show or hide widgets in an EZInput widget based on a visibility dictionary.

    Parameters
    ----------
    ezwidget : EZInput
        The EZInput widget containing elements to show/hide.
    visibility_dictionary : dict
        Dictionary mapping widget names to booleans (True to show, False to hide).

    Returns
    -------
    None
    """
    for widgetname in visibility_dictionary.keys():
        if visibility_dictionary[widgetname]:
            ezwidget[widgetname].layout.display = "inline-flex"
        else:
            ezwidget[widgetname].layout.display = "None"

def _unstyle_widgets(ezwidget, visibility_dictionary):
    for wgt in ezwidget.elements.keys():
        visibility_dictionary[wgt] = True
        if isinstance(ezwidget[wgt], widgets.Button):
            ezwidget.elements[wgt].layout = widgets.Layout(width="50%", display="inline-flex", align_items="center", justify_content="center")
        else:
            ezwidget.elements[wgt].layout = widgets.Layout(width="50%", display="inline-flex")   

def ui_show_structure(experiment):
    """
    Create a widget for visualising the experiment structure.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing the structure.

    Returns
    -------
    EZInput
        Widget for structure visualisation.
    """
    gui = EZInput(title="Structure")
    def show_structure(widget_elements):
        if experiment.objects_created["structure"]:
            enable_view_widgets(True)
            widget_elements["preview_structure"].clear_output()
            total = experiment.structure.num_assembly_atoms
            hview = widget_elements["hview"].value
            vview = widget_elements["vview"].value
            widget_elements["n_atoms"].disabled = False
            atoms_number = widget_elements["n_atoms"].value
            if total > atoms_number:
                fraction = atoms_number/total
            else:
                fraction = 1.0
            with widget_elements["preview_structure"]:
                display(
                    experiment.structure.show_assembly_atoms(
                        assembly_fraction=fraction,
                        view_init=[vview, hview, 0]
                    )
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
        icon=show_icon,
    )

    def update_plot(value):
        gui["preview_structure"].clear_output()
        show_structure(gui.elements)

    gui.add_int_slider("n_atoms", description="Atoms to use", min=1, max=10000, step=1, value = 1000, on_change=update_plot, continuous_update=False, disabled=True)
    gui.add_int_slider(
        "hview",
        description="Horizontal view",
        min=-90,
        max=90,
        step=1,
        value=0,
        on_change=update_plot,
        continuous_update=False,
    )
    gui.add_int_slider(
        "vview",
        description="Vertical view",
        min=-90,
        max=90,
        step=1,
        value=0,
        on_change=update_plot,
        continuous_update=False,
    )
    def enable_view_widgets(b):
        widgets_visibility["n_atoms"] = True
        widgets_visibility["hview"] = True
        widgets_visibility["vview"] = True
        update_widgets_visibility(gui, widgets_visibility)

    #gui.add_button("show_structure", description="Show structure")
    gui.add_output("preview_structure")
    widgets_visibility = {}
    _unstyle_widgets(gui, widgets_visibility)
    gui["button"].layout = widgets.Layout(width='50%', align_items='center', justify_content='center')
    widgets_visibility["n_atoms"] = False
    widgets_visibility["hview"] = False
    widgets_visibility["vview"] = False
    update_widgets_visibility(gui, widgets_visibility)

    gui["preview_structure"].clear_output()
    return gui

def ui_show_labelled_structure(experiment):
    """
    Create a widget for visualising the labelled structure (with probes).

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing the labelled structure.

    Returns
    -------
    EZInput
        Widget for labelled structure visualisation.
    """
    gui = EZInput(title="Labelled Structure")

    def show_particle(
                    emitter_plotsize = 1, 
                    source_plotsize = 1, 
                    hview=0,
                    vview=0):
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
                enable_view_widgets(True)
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
    gui.add_button(
            "show_labelled_structure",
            description="Show labelled structure",
            icon=show_icon,
        )
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
    
    def enable_view_widgets(b):
        widgets_visibility["emitter_plotsize"] = True
        widgets_visibility["source_plotsize"] = True
        widgets_visibility["hview"] = True
        widgets_visibility["vview"] = True
        update_widgets_visibility(gui, widgets_visibility)
    gui["show_labelled_structure"].on_click(show_labelled_structure)
    gui.add_output("preview_labelled_structure")
    widgets_visibility = {}
    _unstyle_widgets(gui, widgets_visibility)
    widgets_visibility["emitter_plotsize"] = False
    widgets_visibility["source_plotsize"] = False
    widgets_visibility["hview"] = False
    widgets_visibility["vview"] = False
    update_widgets_visibility(gui, widgets_visibility)
    return gui

def ui_show_virtual_sample(experiment):
    """
    Create a widget for visualising the virtual sample.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing the virtual sample.

    Returns
    -------
    EZInput
        Widget for virtual sample visualisation.
    """
    gui = EZInput(title="Virtual Sample")

    def update_plot(change):
        gui["preview_virtual_sample"].clear_output()
        with gui["preview_virtual_sample"]:
            if experiment.generators_status("coordinate_field"):
                widgets_visibility["horizontal_view"] = True
                widgets_visibility["vertical_view"] = True
                update_widgets_visibility(gui, widgets_visibility)
                gui["horizontal_view"].disabled = False
                gui["vertical_view"].disabled = False
                hview = gui["horizontal_view"].value
                vview = gui["vertical_view"].value
                with gui["preview_virtual_sample"]:
                    display(experiment.coordinate_field.show_field(
                        view_init=[vview, hview, 0],
                        return_fig=True))
                    plt.close()
            else:
                display("Virtual sample not created yet, please create it first.")
    gui.add_button(
        "show_virtual_sample",
        description="Show virtual sample",
        icon=show_icon,
    )
    gui.add_int_slider(
        "horizontal_view",
        description="Rotation angle (degrees)",
        min=-90,
        max=90,
        step=1,
        value=0,
        continuous_update=False,
        on_change=update_plot,
        disabled=True,
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
        disabled=True,
    )
    widgets_visibility = {}
    _unstyle_widgets(gui, widgets_visibility)
    widgets_visibility["horizontal_view"] = False
    widgets_visibility["vertical_view"] = False
    update_widgets_visibility(gui, widgets_visibility)
    gui["show_virtual_sample"].on_click(update_plot)
    gui.add_output("preview_virtual_sample")
    return gui

def ui_show_modality(experiment):
    """
    Create a widget for visualising the point spread function (PSF) of selected imaging modalities.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing imaging modalities.

    Returns
    -------
    EZInput
        Widget for modality PSF visualisation.
    """
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
    """
    Create a widget for setting and previewing acquisition parameters.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing acquisition settings.

    Returns
    -------
    EZInput
        Widget for acquisition parameter selection and preview.
    """
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
        acquisition_gui["current_parameters"].value = _mods_text_update(
            mods_text_base="Current acquisition parameters for modalities: ",
            mod_acq_params=experiment.selected_mods,
            keys_to_use=["exp_time", "noise", "nframes"],
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

    def preview_volume(b):
        def get_vol_preview(imaging_system, acq_gui):
            def get_volume(Modality, Exposure, Noise):
                fig = plt.figure()
                for single_channel in imager_channels:
                    single_mod_acq_params = dict(
                        exp_time=Exposure,
                        noise=Noise,
                        save=False,
                        nframes=1,
                        channel=single_channel,
                        convolution_type="raw_volume"
                    )
                    with io.capture_output() as captured:
                        images_volumes, _beads = (
                            imaging_system.generate_imaging(
                                modality=Modality, **single_mod_acq_params
                            )
                        )
                        volume = images_volumes[0]
                    return volume

            volume = get_volume(
                Modality=acq_gui["modalities_dropdown"].value,
                Exposure=acq_gui["Exposure"].value,
                Noise=acq_gui["Noise"].value,
            )
            acquisition_gui["image_output"].clear_output()
            slider = widgets.IntSlider(
                value=0,
                min=0,
                max=180,
                step=20,
                description='Angle:',
                continuous_update=False,
            )
            toggle_buttons = widgets.ToggleButtons(
                options=['XY', 'XZ', 'YZ'],
                description='plane:',
                disabled=False,
            )
            with acquisition_gui["image_output"]:
                display(
                    widgets.interactive(
                        plot_projection,
                        stack=widgets.fixed(volume),
                        angle=slider,
                        plane=toggle_buttons,
                        method=widgets.fixed("max"),
                    )   
                )
        get_vol_preview(experiment.imager, acquisition_gui)             

    def clear(b):
        print("Acquisition parameters cleared")
        with io.capture_output() as captured:
            experiment.reset_to_defaults(module="acquisitions", save=True)
        acquisition_gui["current_parameters"].value = _mods_text_update(
            mods_text_base="Current acquisition parameters for modalities: ",
            mod_acq_params=experiment.selected_mods,
            keys_to_use=["exp_time", "noise", "nframes"],
        )
        

    def preview_params_chage(change):
        if acquisition_gui["show_preview"].value:
            if acquisition_gui["show_as_volume"].value:
                acquisition_gui["image_output"].clear_output()
                preview_volume(True)
            else:
                acquisition_gui["image_output"].clear_output()
                preview_mod(True)
        else:
            acquisition_gui["image_output"].clear_output()
            

    def _mods_text_update(mods_text_base, mod_acq_params, keys_to_use = ["exp_time", "noise"]):
        mods_text = "<b>"+mods_text_base+"</b>" + "<br>"
        for modality_name, acq_params in mod_acq_params.items():
            if acq_params is None:
                acq_params = "Default"
            else:
                keys_subset = {key: acq_params[key] for key in keys_to_use}
                acq_params["exp_time"] = round(keys_subset["exp_time"], 3)
            mods_text +=  modality_name + ": " + "&emsp;" +  str(keys_subset) + "<br>"
        return mods_text

    acquisition_gui.add_label("Set acquisition parameters")
    selected_mods = list(experiment.imaging_modalities.keys())
    acquisition_gui.add_dropdown("modalities_dropdown", options=selected_mods)
    acquisition_gui.add_checkbox("Noise", description="Use Noise", value=True)
    acquisition_gui.add_checkbox(
        "Channels", description="Use all channels", value=True)
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
        value=0.001,
    )   
    acquisition_gui["modalities_dropdown"].observe(
        preview_params_chage, names="value"
    )
    acquisition_gui["Noise"].observe(preview_params_chage, names="value")
    acquisition_gui["Exposure"].observe(preview_params_chage, names="value")
    acquisition_gui.add_HTML(
        "current_parameters",
        _mods_text_update(
            mods_text_base="Current acquisition parameters for modalities: ",
            mod_acq_params=experiment.selected_mods,
            keys_to_use=["exp_time", "noise", "nframes"],
        ),
    )
    acquisition_gui.add_button("Set", description="Update acquisition parameters", 
                               icon=update_icon, style={"button_color": update_colour})
    acquisition_gui.add_button("Clear", description="Reset params", 
                               icon=reset_icon)
    acquisition_gui["Set"].on_click(set_params)
    acquisition_gui["Clear"].on_click(clear)
    acquisition_gui.add_checkbox(
        "show_preview",
        description="Show interactive acquisition preview",
        value=False,
        on_change=preview_mod,
        continuous_update=False,
    )
    acquisition_gui.add_checkbox(
        "show_as_volume",
        description="Show interactive acquisition convolved volume (experimental feature)",
        value=False,
        on_change=preview_volume,
        continuous_update=False,
    )
    acquisition_gui["show_preview"].observe(preview_params_chage, names="value")
    acquisition_gui["show_as_volume"].observe(preview_params_chage, names="value")
    acquisition_gui.add_output("image_output")
    
    acq_widgets = {}
    _unstyle_widgets(acquisition_gui, acq_widgets)
    acquisition_gui.show()
    acq_widgets["Frames"] = True
    acq_widgets["Set"] = True
    acq_widgets["image_output"] = True
    acq_widgets["label_1"] = True
    acq_widgets["modalities_dropdown"] = True
    acq_widgets["Exposure"] = True
    acq_widgets["Noise"] = True
    acq_widgets["Clear"] = True
    acq_widgets["show_preview"] = True
    acq_widgets["show_as_volume"] = True
    
    update_widgets_visibility(acquisition_gui, acq_widgets)
    return acquisition_gui

def ui_preview_results(experiment):
    """
    Create a widget for previewing experiment results.

    Parameters
    ----------
    experiment : ExperimentParametrisation
        The experiment object containing results.

    Returns
    -------
    EZInput
        Widget for previewing experiment results.
    """
    gui = EZInput(title="Preview Results")

    if len(experiment.results.keys()) == 0:
        results_options = ["No results available",]
    else:
        results_options = list(experiment.results.keys())
    value_modalities = results_options[0]  # Set initial value
    def show_results(b):
        gui["preview_results"].clear_output()
        if len(experiment.results.keys()) == 0:
            with gui["preview_results"]:
                display("Experiment not run yet, please run it first.")
        else:
            gui["modality"].disabled = False
            gui["modality"].options = list(experiment.results.keys())
            gui["modality"].value = list(experiment.results.keys())[0]  # Set initial value
            gui["modality"].layout.display = "inline-flex"
            gui["modality"].observe(update_plot, names="value")
            with gui["preview_results"]:
                update_plot(True)
    gui.add_HTML(
        "section_header",
        "Preview Results of the Experiment",
        style={"text-align": "center", "margin-bottom": "20px"}
    )
    #gui.add_label("Preview Results of the Experiment")
    gui.add_button(
        "show_results",
        description="Show Results",
    )
    
    def update_plot(change):
        modality = gui["modality"].value
        image = experiment.results[modality]
        if image.ndim == 3:
            image = image[0]
        figure, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image, cmap='gray')
        ax.axis('off')
        ax.set_title(f"Preview of {modality} results")
        plt.close()
        gui["preview_results"].clear_output()  
        with gui["preview_results"]:
            display(figure)

    gui.add_dropdown(
        "modality",
        description="Modality",
        options=results_options,
        value=value_modalities,
        disabled=True,
    )
    gui["modality"].layout = widgets.Layout(width="50%", display="None")
    gui.add_output("preview_results") 
    gui["show_results"].on_click(show_results)
    return gui
