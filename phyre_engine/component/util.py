"""
This module contains components that do system tasks rather than bioinformatics
tasks.
"""
import os
from pathlib import Path

from phyre_engine.component.component import Component


class ChangeDir(Component):
    """
    Change the process working directory. The new directory must be given in the
    ``working_dir`` key of the pipeline state.
    """
    REQUIRED = ["working_dir"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Change directory to ``working_dir``."""
        working_dir = self.get_vals(data)
        os.chdir(working_dir)
        return data

class MakeDir(Component):
    """
    Make a new directory. The new directory must be given in the ``working_dir``
    key of the pipeline state.

    :param bool exist_ok: If `False`, an error will be raised if the target
        directory already exists.

    :param bool parents: If `True`, create all parents of the directory.
    """
    REQUIRED = ["working_dir"]
    ADDS = []
    REMOVES = []

    def __init__(self, exist_ok=True, parents=True):
        self.exist_ok = exist_ok
        self.parents = parents

    def run(self, data, config=None, pipeline=None):
        """Create a new directory."""
        working_dir = self.get_vals(data)
        Path(working_dir).mkdir(parents=self.parents, exist_ok=self.exist_ok)
        return data
