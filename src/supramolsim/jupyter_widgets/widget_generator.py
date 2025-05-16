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