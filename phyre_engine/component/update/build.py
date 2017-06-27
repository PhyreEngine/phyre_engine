from phyre_engine.component.component import Component
import logging
import subprocess
from pathlib import Path

import yaml
try:
    from yaml import CSafeLoader as SafeLoader, CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader, SafeDumper

log = lambda: logging.getLogger(__name__)

class Build(Component):
    """
    Build all tools in the ``update_required`` list contained within the
    pipeline state. The following keys may be used in each
    command with standard python string formatting:

    archive_name:
        Name (without directories) of the saved archive.

    archive_path:
        Full path of the saved archive.

    install_dir:
        Directory in which to install the software (the ``install_dir``
        argument to the constructor).

    new_version:
        Version to be installed.

    .. note::
        Python string formatting is applied to each command, so if you wish
        to use a brace (``{`` or ``}``), you must double them, so ``{``
        would become ``{{``.

    Commands should be supplied as a dictionary of lists of strings, indexed by
    the name of program to build.

    .. warning::
        For flexibility, commands are executed with a shell. You should *never*
        supply build commands from untrusted sources, as it will allow an
        attacker to fake full control of your system.

    :param str install_dir: Installation directory. This may be used as the
        ``install_dir`` key in each build command.
    :param build_commands: Commands to run.
    :type build_commands: Dictionary indexed by program name, containing a list
        of build commands to run.
    """
    REQUIRED = ["update_required"]
    ADDS = []
    REMOVES = []

    def __init__(self, install_dir, build_commands):
        self.install_dir = Path(install_dir)
        self.build_commands = build_commands

    def build(self, tool):
        """Execute build commands for tool."""
        # Start executing build commands.
        name = tool["name"]
        if name not in self.build_commands:
            return

        for cmd in self.build_commands[name]:
            formatted_cmd = cmd.format(**{k: str(v) for k, v in tool.items()})
            log().info("Running `%s' for tool `%s'", formatted_cmd, name)
            subprocess.run([formatted_cmd], shell=True, check=True)

    def run(self, data, config=None, pipeline=None):
        for tool in data["update_required"]:
            tool["install_dir"] = str(self.install_dir)
            self.build(tool)
        return data

class ChangeVersion(Component):
    """
    Update the version database with the new version number of the Installed
    software.

    The pipeline state must contain a list in the ``update_required`` key, each
    element of which must have a ``name`` and ``new_version`` element.

    :param str version_db: File path of the version database.
    """
    ADDS = []
    REMOVES = ["update_required"]
    REQUIRED = ["update_required"]

    def __init__(self, version_db):
        self.version_db = Path(version_db)

    def update_db(self, tools):
        """Update the version database with the new version."""
        with self.version_db.open("r+") as db_fh:
            version_db = yaml.load(db_fh, SafeLoader)
            if version_db is None:
                version_db = {}
            for tool in tools:
                version_db[tool["name"]] = tool["new_version"]
            db_fh.truncate()
            yaml.dump(version_db, db_fh, SafeDumper, default_flow_style=False)

    def run(self, data, config=None, pipeline=None):
        self.update_db(data["update_required"])
        del data["update_required"]
        return data
