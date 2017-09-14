"""
This module contains components that do system tasks rather than bioinformatics
tasks.
"""
import os
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
