import yaml
import os
import copy


def load_yaml(yaml_file: str):
    """
    Returns the associated dictionary from a yaml
    """
    with open(yaml_file, "r") as f:
        config_dictionary = yaml.safe_load(f)
    return config_dictionary


def save_yaml(data=None, name="NONAME", output_directory=None, as_str=False):
    """
    Saves a dictionary to a yaml file
    """
    file_name = name + ".yml"
    file_path = os.path.join(output_directory, file_name)
    #yaml.Dumper.ignore_aliases = lambda self, data: True
    if as_str:
        copy_of_params = copy.deepcopy(data)
        for combination_id, list_of_parameters in copy_of_params.items():
            for parameter in list_of_parameters:
                if type(parameter) is dict:
                    for parameter_name in parameter.keys():
                        parameter[parameter_name] = str(parameter[parameter_name])
        with open(file_path, "w") as f:
            yaml.safe_dump(copy_of_params, f, default_flow_style=False)
    else:
        with open(file_path, "w") as f:
            yaml.safe_dump(data, f, default_flow_style=False)