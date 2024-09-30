import numpy as np


def builder_format(targets: dict,
                reference_point: np.ndarray,
                scale: float = 1e-10,
                axis: dict = dict(pivot=None, direction=None),
                info: str = "No info",
                **kwargs):
    """
    targets is a dictionary with at least two keys: "coordinates" and "normals"
    """
    instance_builder = dict(
        targets=targets,
        reference_point=reference_point,
        scale=scale,
        axis=axis,
        info=info,
    )
    for key, value in kwargs.items():
        instance_builder[key] = value
    return instance_builder


def label_builder_format(label_id,
                        fluorophore_id,
                        labelling_efficiency=1):
    label_info = dict(label_id=label_id, fluorophore_id=fluorophore_id, labelling_efficiency=labelling_efficiency)
    return label_info


def struct_params_format(struct_ID=None,
                        title=None,
                        labels=list(["<None>"])):
    structure_params = dict(
        id=struct_ID,
        title=title,
        labels=labels
    )
    return structure_params