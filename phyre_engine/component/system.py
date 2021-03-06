"""
This module contains components that allow arbitrary system commands to be run.
"""
from phyre_engine.component.component import Component
import subprocess

class System(Component):
    """

    Run an arbitrary system command. The specified command will be formatted
    using Python's string formatting, making available all elements of the
    pipeline state.

    To help avoid shell injection and provide a vague semblance of security,
    this component does not run commands within a shell. Individual arguments
    must be supplied as separate elements of a list. If you really need to run
    commands with a shell, set the `shell` parameter to `True`.

    .. note::

        * The commands processed by this component will be processed by Python's
          :py:func:`str.format` method. To include literal braces (``{``) you
          must escape them by doubling: ``{{surrounded by literal braces}}``.

        * If you find yourself consistently using this component for a task,
          consider getting in touch (ideally with some code) with the project
          maintainers so that the task can be included within Phyre Engine.

    .. warning::

        This command will run arbitrary commands. It is extremely dangerous. If
        you are using the pipeline state to run commands, be aware that this is
        indirectly user-generate data and should *not be trusted!*

    >>> import phyre_engine.component.system as system
    >>> sys_cmpt = system.System([
    ...     ["/usr/bin/do_thing", "arg1", "arg2"],
    ...     # other commands
    ... ])
    >>> sys_cmpt = system.System(
    ...     ["echo 'shell features' >/dev/null"],
    ...     shell=True)

    :param list[list[str]] commands: A list of commands to run.
    :param bool shell: Whether to run each command in a shell.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, commands, shell=False):
        self.commands = commands
        self.shell = shell

    def run(self, data, config=None, pipeline=None):
        """Running system commands."""
        for cmd in self.commands:
            if isinstance(cmd, str):
                cmd = [cmd]
            formatted_cmd = [arg.format(**data) for arg in cmd]
            subprocess.run(formatted_cmd, shell=self.shell)
        return data

class Python(Component):
    """
    Evaluate arbitrary Python code.

    The `code` parameter is executed as though it is the :py:meth:`~.run`
    method of this component. The local variables ``data``, ``config`` and
    ``pipeline`` are set from the arguments of :py:meth:`~.run`.

    If the optional parameters `adds`, `removes` or `required` are set, each is
    interpreted as Python code describing the fields added, removed or required
    by the pipeline state.
    """

    @property
    def REQUIRED(self):
        return eval(self.required) if self.required is not None else []

    @property
    def ADDS(self):
        return eval(self.adds) if self.adds is not None else []

    @property
    def REMOVES(self):
        return eval(self.removes) if self.removes is not None else []

    def __init__(self, code, adds=None, removes=None, required=None):
        self.code = code
        self.adds = adds
        self.removes = removes
        self.required = required

    def run(self, data, config=None, pipeline=None):
        """Evaluating Python code."""
        exec(self.code, {}, {
            "data": data,
            "config": config,
            "pipeline": pipeline
        })
        return data
