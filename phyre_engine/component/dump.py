"""
This module contains components for dumping the pipeline state in various
formats. By default, each dumper will write to standard output. The ``output``
parameter may be set to either a file name or output stream to determine where
the output will be written.
"""
import json
import sys
from phyre_engine.component.component import Component
from phyre_engine.tools.util import Stream

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeDUmper as SafeDumper
except ImportError:
    from yaml import SafeDumper

class Dumper(Component):
    """
    Base class for dumpers.

    :param output: A file name, stream or :py:class:`pathlib.Path` object.
    """
    # Disable warnings for this half-implemented class.
    # pylint: disable=abstract-method

    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, output=sys.stdout):
        self.output = output

class Json(Dumper):
    """Dump the pipeline state in JSON."""
    def run(self, data, config=None, pipeline=None):
        """Dump pipeline state."""
        with Stream(self.output, "w") as out_fh:
            json.dump(data, out_fh)
        return data

class Yaml(Dumper):
    """Dump pipeline state in YAML."""
    def run(self, data, config=None, pipeline=None):
        with Stream(self.output, "w") as out_fh:
            yaml.dump(data, out_fh, SafeDumper)
