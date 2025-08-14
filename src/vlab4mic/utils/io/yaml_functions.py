import yaml
import os


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