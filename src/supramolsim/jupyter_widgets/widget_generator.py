import ipywidgets as widgets
from dataclasses import dataclass, field, fields
import numpy as np


@dataclass
class widgen:
    int_minmaxstep = [0, 10, 1]
    float_minmaxstep = [0, 10, 0.1]
    layout_config = {'width': '50%'}

    def gen_range_slider(
        self,
        slidertype="int",
        minmaxstep=None,
        options: list = None,
        description="range_slider",
        disabled=False,
        orientation="horizontal",
        layout=None,
        **kwargs
    ):
        if layout is None:
            layout = self.layout_config
        if slidertype == "int":
            if minmaxstep is None:
                minmaxstep = self.int_minmaxstep
            range_slider = widgets.IntRangeSlider(
                value=[minmaxstep[0], minmaxstep[1]],
                min=minmaxstep[0],
                max=minmaxstep[1],
                step=minmaxstep[2],
                description=description,
                disabled=False,
                orientation=orientation,
                layout=layout,
                **kwargs
            )
        if slidertype == "float":
            if minmaxstep is None:
                minmaxstep = self.float_minmaxstep
            range_slider = widgets.FloatRangeSlider(
                value=[minmaxstep[0], minmaxstep[1]],
                min=minmaxstep[0],
                max=minmaxstep[1],
                step=minmaxstep[2],
                description=description,
                disabled=False,
                orientation=orientation,
                layout=layout,
                **kwargs
            )
        return range_slider

    def gen_logicals(self, layout=None, **kwargs):
        if layout is None:
            layout = self.layout_config
        logicals =  widgets.RadioButtons(
            options=["True", "False", "Both"],
            layout=layout, # If the items' names are long
            description="",
            disabled=False
            )
        return logicals

    def gen_bound_int(self, value = 2, max = 100, **kwargs):
        bound_int = widgets.BoundedIntText(
            value=value,
            min=0,
            max=max,
            step=1,
            **kwargs
        )
        return bound_int
    
    def gen_box(self, widget1=None, widget2=None, orientation="horizontal", **kwargs):
        items = [widget1, widget2]
        if orientation == "horizontal":
            box = widgets.HBox(items, **kwargs)
        else:
            box = widgets.VBox(items, **kwargs)
        return box
    
    def gen_dropdown(self, options=None):
        menu = widgets.Dropdown(options=options)
        return menu
    
    def gen_box_linked(self, w1=None, w2=None, dependant=None, observed=None, orientation="horizontal", update_method = None, update_params = None, **kwargs):
        def update(change):
            print("updating")
            update_method(change.new, update_params)
        observed.observe(update, names="value")
        items = [w1, w2]
        if orientation == "horizontal":
            box = widgets.HBox(items, **kwargs)
        else:
            box = widgets.VBox(items, **kwargs)
        return box

    def gen_interactive_dropdown(self,
                                 options=None,
                                 orientation="horizontal",
                                 routine=None,
                                 height = '400px',
                                 **kwargs
                                 ):
        params_widgets = dict()
        list_of_paramwidgets = []
        for keyname, val in kwargs.items():
            wtype = val[0]
            wparams = val[1]
            if wtype == "float_slider":
                params_widgets[keyname] = widgets.FloatSlider(
                    value=wparams[0],
                    min=wparams[1],
                    max=wparams[2],
                    step=wparams[3],
                    description = keyname,
                    continuous_update=False
                )
            if wtype == "int_slider":
                params_widgets[keyname] = widgets.IntSlider(
                    value=wparams[0],
                    min=wparams[1],
                    max=wparams[2],
                    step=wparams[3],
                    description = keyname,
                    continuous_update=False
                )
        drop = self.gen_dropdown(options=options)
        list_of_paramwidgets.append(drop)
        #
        def func(dropdown, **kwargs2):
            r_out = routine(dropdown, **kwargs2)
        
        funct_dictionary = {'dropdown': drop}
        #
        if len(params_widgets.keys()) > 0:
            for key, wid in params_widgets.items():
                funct_dictionary[key] = wid
                list_of_paramwidgets.append(wid)
        params = widgets.VBox(list_of_paramwidgets)
        params.layout = widgets.Layout(height='auto', width='auto')
        out = widgets.interactive_output(func, funct_dictionary) 
        out.layout = widgets.Layout(height='auto', width='auto')
        box = self.gen_box(widget1=params, widget2=out, orientation=orientation)
        return box
    

    def gen_action_with_options(self, param_widget = None, options = None, routine=None, action_name = "action", height = '400px', **kwargs):
        params_widgets = dict()
        action_result = None
        for keyname, val in kwargs.items():
            wtype = val[0]
            wparams = val[1]
            if wtype == "float_slider":
                params_widgets[keyname] = widgets.FloatSlider(
                    value=wparams[0],
                    min=wparams[1],
                    max=wparams[2],
                    step=wparams[3],
                    description = keyname,
                    readout_format='.4f',
                    style = {'description_width': 'initial'}
                )
            if wtype == "int_slider":
                params_widgets[keyname] = widgets.IntSlider(
                    value=wparams[0],
                    min=wparams[1],
                    max=wparams[2],
                    step=wparams[3],
                    description = keyname,
                    style = {'description_width': 'initial'}
                )
            if wtype == "checkbox":
                params_widgets[keyname] = widgets.Checkbox(
                    value=wparams,
                    description = keyname,
                )
            if wtype == "button":
                params_widgets[keyname] = widgets.Button(
                    description = wparams[0],
                    layout=widgets.Layout(width='100%')
                )
                params_widgets[keyname].on_click(wparams[1])
        list_of_paramwidgets = []
        if len(params_widgets.keys()) > 0:
            for key, wid in params_widgets.items():
                #funct_dictionary[key] = wid
                list_of_paramwidgets.append(wid)
        
        def action(b):
            widget_values = {}
            if len(params_widgets.keys()) > 0:
                for key, wid in params_widgets.items():
                    if wid._model_name != 'ButtonModel':
                        widget_values[key] = wid.value
            out_static.clear_output()
            with out_static:
                if options:
                    r_out = routine(param_widget, options, **widget_values)
                    display(r_out)
                else:
                    r_out = routine(param_widget, **widget_values)
                    display(r_out)
        out_static = widgets.Output()
        out_static.layout=widgets.Layout(height='auto', width='auto')
        trigger = widgets.Button(description=action_name, layout=widgets.Layout(width='100%'))
        list_of_paramwidgets.append(trigger)
        params = widgets.VBox(list_of_paramwidgets)
        params.layout=widgets.Layout(height='auto', width='auto')
        trigger.on_click(action)
        preview_area = self.gen_box(
            widget1=params,
            widget2=out_static,
            orientation="vertical")
        return preview_area