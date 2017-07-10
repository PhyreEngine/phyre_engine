"""
Module containing a base class for running external tools.
"""
from pathlib import Path
from enum import Enum

class Positional(Enum):
    """Possible locations for positional arguments."""
    FIRST = 1
    LAST = 2

class ExternalTool:
    """
    This class provides provides factories for building command lines for
    external tools.

    This class attempts to handle the following:

    *   Finding the correct executable when passed an executable name and an
        optional directory. If no directory is given, then the system PATH is
        searched as usual; if a directory is supplied, then the executable must
        be in that directory.

    *   Mapping human-readable argument names to the actual names. For example,
        mapping the string ``database`` to the flag ``-d``.

    *   Automatically adding leading dashes to options.

    :param dict[str, str] flag_map: Mapping of human-readable names to argument
        names.
    :param str short_prefix: Prefix used for short flags.
    :param str long_prefix: Prefix used for long flags.
    :param str join_options: Join flags and values with this. For example, set
        this to "``=``" to get flags like ``--foo=bar``. The default value of
        ``False`` indicates that flags and values are separate arguments.
    :param Positional positional_location: Where to place positional arguments
        in relation to the optional arguments.

    For example, let's create a simple wrapper around the ``ls`` command. Let's
    say I find the name of the ``--no-preserve-root`` flag to be too cryptic,
    so I will alias it to ``really-unsafe``:

    >>> from phyre_engine.tools.external import ExternalTool
    >>> ls = ExternalTool(
    ...     {"really-unsafe": "no-preserve-root"},
    ...     join_options="=")
    >>> ls(("/bin", "ls"), ["/tmp", "/etc"],
    ...    ["really-unsafe"], {"color": "always"})
        ["/bin/ls", "--no-preserve-root", "--color=always", "/tmp", "/etc"]

    For broken commands that don't use consistent prefixes for long and short
    options, you should process the broken flags in your application code and
    pass them directly with leading hyphens.
    """

    def __init__(
            self, flag_map=None, short_prefix="-", long_prefix="--",
            join_options=False, positional_location=Positional.LAST):

        self.flag_map = flag_map if flag_map is not None else {}
        self.short_prefix = short_prefix
        self.long_prefix = long_prefix
        self.join_options = join_options
        self.positional_location = positional_location

    def __call__(self, executable, positional=None, flags=None, options=None):
        """
        Get a command line for the given tool.

        :param tuple(str, str) executable: Optional directory path, and name of
            the executable to run. This can be a single string, in which case it
            should either be an absolute path or the name of a command on the
            system PATH. If given as a tuple, the first parameter (the
            directory) may be ``None`` to search the system path.

        :param list[str] positional: Positional arguments to append to the
            command line.

        :param list[str] flags: Flags that will be appended to the command line
            with the correct prefix.

        :param dict[str, str] options: Mapping of flags to values that will be
            appended to the command line.
        """

        positional = positional if positional is not None else []
        flags = flags if flags is not None else []
        options = options if options is not None else {}

        command_line = []
        positional = [str(arg) for arg in positional]

        # Executable name
        if isinstance(executable, str):
            command_line.append(executable)
        elif executable[0] is None:
            command_line.append(executable[1])
        else:
            command_line.append(str(Path(executable[0], executable[1])))

        optional_args = []
        # Apply flag_map to all non-positional arguments and add them to the
        # command line.
        for flag in flags:
            optional_args.append(self._remap_flag(flag))

        for flag, value in options.items():
            value = str(value)
            flag = self._remap_flag(flag)
            if self.join_options:
                optional_args.append(self.join_options.join([flag, value]))
            else:
                optional_args.extend([flag, value])

        # Add positional arguments
        if self.positional_location == Positional.LAST:
            command_line.extend(optional_args + positional)
        elif self.positional_location == Positional.FIRST:
            command_line.extend(positional + optional_args)

        return command_line

    def _remap_flag(self, flag):
        if (flag.startswith(self.short_prefix)
                or flag.startswith(self.long_prefix)):
            return flag
        elif flag in self.flag_map:
            flag = self.flag_map[flag]
        return self._prefix(flag)

    def _prefix(self, flag):
        if len(flag) == 1:
            return self.short_prefix + flag
        else:
            return self.long_prefix + flag
