import numpy as np


def format_coordinates(
    coords, plotsize=1, plotalpha=1, plotmarker="o", plotcolour="k", name="", **kwargs
):
    formatted = {
        "coordinates": coords,
        "plotsize": plotsize,
        "plotalpha": plotalpha,
        "plotmarker": plotmarker,
        "plotcolour": plotcolour,
        "label_name": name,
    }
    for key, value in kwargs.items():
        formatted[key] = value
    return formatted


def set_colorplot(plotting_params):
    # this function takes the dictionary that is stored in self.label_targets
    if len(plotting_params) == 0:
        defaultcolour = "#984ea3"
    else:
        unavailable = []
        for labelname in plotting_params.keys():
            # print(f'colours taken: {plotting_params[labelname]["plotcolour"]}')
            unavailable.append(plotting_params[labelname]["plotcolour"])
        CB_color_cycle = [
            "#377eb8",
            "#ff7f00",
            "#4daf4a",  # https://gist.github.com/thriveth/8560036
            "#f781bf",
            "#a65628",
            "#984ea3",
            "#999999",
            "#e41a1c",
            "#dede00",
        ]
        defaultcolour = "#984ea3"
        while defaultcolour in unavailable:
            defaultcolour = np.random.choice(CB_color_cycle)
    return defaultcolour
