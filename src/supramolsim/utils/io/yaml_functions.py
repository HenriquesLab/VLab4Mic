import yaml


def load_yaml(yaml_file: str):
    """
    Returns the associated dictionary from a yaml
    """
    with open(yaml_file,'r') as f:
        config_dictionary = yaml.safe_load(f)
    return config_dictionary
