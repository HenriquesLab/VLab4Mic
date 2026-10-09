from importlib.metadata import PackageNotFoundError, version

from .utils.io.yaml_functions import load_yaml

try:
    __version__ = version("vlab4mic")
except PackageNotFoundError:  # running from a source checkout
    __version__ = "unknown"
