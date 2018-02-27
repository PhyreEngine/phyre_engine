"""
This module contains components that do system tasks rather than bioinformatics
tasks.
"""
import os
from pathlib import Path

from phyre_engine.component.component import Component


class SaveRootDir(Component):
    """Save the current directory in the ``root_directory`` key."""

    ADDS = ["root_directory"]
    REMOVES = []
    REQUIRED = []

    def run(self, data, config=None, pipeline=None):
        """Save current directory."""
        data["root_directory"] = os.getcwd()
        return data


class RestoreRootDir(Component):
    """Change back to the ``root_directory``."""

    ADDS = []
    REMOVES = []
    REQUIRED = ["root_directory"]

    def run(self, data, config=None, pipeline=None):
        """Restore root directory."""
        os.chdir(data["root_directory"])
        return data


class ChangeDir(Component):
    """
    Change the process working directory. The new directory must be given in the
    ``working_dir`` key of the pipeline state.

    This component will push the previous working directory onto the
    ``directory_stack`` list in the pipeline state. The previous working
    directory can be returned to with :py:class:`.PopDir`.
    """
    REQUIRED = ["working_dir"]
    ADDS = ["directory_stack"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Change directory to ``working_dir``."""
        cwd = os.getcwd()
        data.setdefault("directory_stack", []).append(cwd)
        working_dir = self.get_vals(data)
        os.chdir(working_dir)
        return data


class PopDir(Component):
    """
    Reset working directory, or if no change in working directory was recorded
    do nothing.

    This component will attempt to pop a directory off of the
    ``directory_stack`` list. If that list does not exist or is empty, the
    working directory is not changed. If an element could be popped, then the
    working directory is changed to that directory.  If this leaves
    ``directory_stack`` empty, it will be removed from the pipeline state.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def run(self, data, config=None, pipeline=None):
        """Go back to a previous directory."""
        if "directory_stack" not in data or not data["directory_stack"]:
            self.logger.info("Directory stack empty -- staying put.")
        else:
            prev_dir = data["directory_stack"].pop()
            self.logger.info("Moving to '%s'", prev_dir)
            os.chdir(prev_dir)

        if "directory_stack" in data and not data["directory_stack"]:
            del data["directory_stack"]
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
