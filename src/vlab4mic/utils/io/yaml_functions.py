import yaml
import os
import copy
import numpy as np


def load_yaml(yaml_file: str):
    """
    Returns the associated dictionary from a yaml
    """
    with open(yaml_file, "r") as f:
        config_dictionary = yaml.safe_load(f)
    return config_dictionary


def save_yaml(data=None, name="NONAME", output_directory=None):
    """
    Saves a dictionary to a yaml file
    """
    file_name = name + ".yml"
    file_path = os.path.join(output_directory, file_name)
    #yaml.Dumper.ignore_aliases = lambda self, data: True
    with open(file_path, "w") as f:
        yaml.safe_dump(data, f, default_flow_style=False)

def simulation_report(positions):
    """
    Parameters actually used to simulate each image, with defaults and notes.

    Parameters
    ----------
    positions : dict
        {modality: {channel: {fluorophore: record}}}, as returned by
        Imager.get_positions().

    Returns
    -------
    dict
        {modality: {channel: {fluorophore: {"parameters", "defaults_used",
        "notes"}}}}. "defaults_used" lists config keys that were not set and
        took a default; "notes" lists computations that were skipped.
        Coordinate arrays are left out, and so are records without any of
        these entries.
    """
    keys = ("parameters", "defaults_used", "notes")
    report = {}
    for modality, channels in (positions or {}).items():
        for channel, fluorophores in channels.items():
            for fluo, record in fluorophores.items():
                entry = {k: copy.deepcopy(record.get(k)) for k in keys if record.get(k)}
                if entry:
                    report.setdefault(modality, {}).setdefault(channel, {})[fluo] = entry
    return report
