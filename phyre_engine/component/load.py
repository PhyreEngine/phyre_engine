"""Module containing components for loading a pipeline state from a file."""
import json
import pickle
try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeLoader as SafeLoader
except ImportError:
    from yaml import SafeLoader

from phyre_engine.component import Component
from phyre_engine.tools.util import Stream

class Loader(Component):
    """
    Base class of loaders. Loaded data is merged into the pipeline state,
    overwriting any fields already present.

    :param input_source: A file name, stream or :py:class:`pathlib.Path` object.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, input_source):
        self.input_source = input_source

class Json(Loader):
    """Load the pipeline state from JSON data."""
    def run(self, data, config=None, pipeline=None):
        """Dump pipeline state."""
        with Stream(self.input_source, "r") as in_fh:
            loaded = json.load(in_fh)
            data.update(loaded)
        return data

class Yaml(Loader):
    """Load pipeline state from YAML data."""
    def run(self, data, config=None, pipeline=None):
        with Stream(self.input_source, "r") as in_fh:
            loaded = yaml.load(in_fh, SafeLoader)
            data.update(loaded)
        return data

class Pickle(Loader):
    """Load pipeline state from a pickled state."""

    def run(self, data, config=None, pipeline=None):
        """Load pipeline state from pickle."""
        with Stream(self.input_source, "rb") as in_fh:
            loaded = pickle.load(in_fh)
            data.update(loaded)
        return data
