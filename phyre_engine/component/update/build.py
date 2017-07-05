from phyre_engine.component.component import Component
import logging
import subprocess
from pathlib import Path
import collections.abc

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
                version_db[tool["name"]] = str(tool["new_version"])
            db_fh.seek(0)
            db_fh.truncate()
            yaml.dump(version_db, db_fh, SafeDumper, default_flow_style=False)

    def run(self, data, config=None, pipeline=None):
        self.update_db(data["update_required"])
        del data["update_required"]
        return data

class UpdateConfigFile(Component):
    """
    Update a configuration file for each element of the pipeline's
    ``update_required`` element.

    The configuration file is assumed to be a YAML file. Each element in the
    pipeline state's ``update_required`` list is examined and the fields
    specified in the constructor's ``update_fields`` parameter are updated.

    For example, the pipeline state may look like this:

    ::
        state = {
            # Unrelated pipeline state...
            update_required: [
                {"name": "foo", "new_version": "1.2", "install_dir": "/opt"}
                # Any other tools
            ]
        }

    An :py:class:`~.UpdateConfigFile` object may be created to update the
    details of the ``foo`` tool in ``config.yml``:

    ::
        update_cfg = UpdateConfigFile(
            "config.yml", {
                "foo": {
                    "dirs": {
                        "bin": "{install_dir}/foo_tool/{new_version}/bin",
                        "lib": "{install_dir}/foo_tool/{new_version}/lib"
                    },
                    "version": "{new_version}"
                }
            })

    When the :py:meth:`~.run` method of ``update_cfg`` is called, it will update
    the file ``config.yml`` to contain the following:

    .. code-block:: yaml
        # Definitions for other tools...
        foo:
            dirs:
                bin: /opt/foo_tool/1.2/bin
                lib: /opt/foo_tool/1.2/lib
            version: 1.2
            # Any remaining data for "foo" with different keys.

    Any strings supplied in the ``update_fields`` data structure will have
    python's string formatting applied. The available values will be the
    contents of the matching element of the ``update_required`` pipeline field.

    .. warning::

        Because each string in the ``update_fields`` section is processed by
        python's string formatting, you *must* escape any literal braces by
        doubling them: ``{{hello world}}`` will be escaped to ``{hello world}``
        in the output file.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["update_required"]

    def __init__(self, config_file, update_fields):
        self.config_file = Path(config_file)
        self.update_fields = update_fields

    def run(self, data, config=None, pipeline=None):
        # Read the existing config, then mutate it as we iterate over tools
        config = {}
        if self.config_file.exists():
            with self.config_file.open("r") as conf_in:
                config = yaml.load(conf_in, SafeLoader)

        for tool in data["update_required"]:
            name = tool["name"]
            if name in self.update_fields:
                formatted = _format_nested_strings(
                    self.update_fields[name],
                    tool)
                # Set the section if it doesn't exist or is not a map
                if name not in config or not isinstance(config[name], dict):
                    config[name] = formatted
                else:
                    config[name].update(formatted)

        # Write config file
        with self.config_file.open("w") as conf_out:
            yaml.dump(config, conf_out, SafeDumper)

        return data

def _format_nested_strings(fields, replacements):
    """Find all string values and format them with kwargs."""

    if isinstance(fields, str):
        fields = fields.format(**replacements)
    elif isinstance(fields, dict):
        # Call to list() because we are going to mutating "fields"
        for k, v in list(fields.items()):
            fields[k] = _format_nested_strings(v, replacements)
    elif isinstance(fields, collections.abc.Iterable):
        fields = [_format_nested_strings(v, replacements) for v in fields]
    return fields
