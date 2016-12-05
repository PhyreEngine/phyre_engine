"""Module containing classes for running hh-suite tools."""
import subprocess

class HHSuiteTool:
    """Base class for hh-suite tools responsible for running processes.

    Each tool in hh-suite accepts several command-line flags. Some of these
    flags have non-obvious names, so this method will convert between
    human-readable flags and the actual command line flags. Subclasses
    should define the dictionary ``flags`` mapping a human-readable version
    to the actual flag. Flags not found in the ``flags`` map will be passed
    straight to the tool after prepending a single hyphen.

    >>> class HHSearchExample(HHSuiteTool)
    >>>     program = "hhsearch"
    >>>     flags = {"database": "d", "input": "i"}
    >>>     pass
    >>> hhs = HHSearchExample(database="uniprot_2012_03",
    >>>        input="in.fasta", cpu=10, nocons=True)
    >>> hhs.command_line
    ["hhsearch", "-d" "uniprot_2012_03", "-i", "in.fasta",
            "-cpu", "10", "-nocons"]
    """

    def __init__(self, **flags):
        """Set up a tool with the given command-line flags."""
        long_flags = flags
        self.command_line = [self.program]
        self._map_flags(long_flags)


    def _map_flags(self, flags):
        """Map human-readable option names to command-line flags."""
        for long_flag, value in flags.items():
            flag = "-"

            if long_flag in self.flags:
                #Convert long names to short names
                flag += self.flags[long_flag]
            else:
                #Or fall back to the short name
                flag += long_flag

            if type(value) == bool:
                #Bool flags do not take a value
                self.command_line.append(flag)
            else:
                #Non-bool flags take an argument
                self.command_line.append(flag)
                self.command_line.append(value)


    def run(self):
        """Execute the tool.

        :raises subprocess.CalledProcessError: A non-zero exit code was returned
                by the tool. Standard error is returned in the ``output``
                attribute.
        """

        #Don't capture standard output. Do capture standard error so that we can
        #package it into an exception if a non-zero exit code was returned.
        process = subprocess.Popen(self.command_line, stdout=None,
                stderr=subprocess.PIPE)
        _, stderr = process.communicate()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                    self.command_line, stderr)


class HHBlits(HHSuiteTool):
    """Wrapper for the hhblits tool included within the hh-suite package.

    As well as the standard short-form flags of hhblits (see ``hhblits -help``),
    the following flags are accepted:

    :param database: HHsearch database, equivalent to the ``-d`` flag.
    :param input: Input alignment file, equivalent to the ``-i`` flag.
    :param iterations: Number of search iterations, equivalent to the
        ``-n`` flag.
    :param evalue_cutoff: e-value cutoff, equivalent to the ``-e`` flag.
    :param output: output file, equivalent to the ``-o`` flag.
    """

    program = "hhblits"
    flags = {
            "database":"d",
            "input": "i",
            "iterations": "n",
            "evalue_cutoff": "e",
            "output": "o"
            }

    def __init__(self, **flags):
        """Set up command-line flags for hhblits."""
        super().__init__(**flags)


class HHSearch(HHSuiteTool):
    """Wrapper for the hhsearch tool included within the hh-suite package.

    As well as the standard short-form flags of hhsearch (see ``hhsearch
    -help``), the following flags are accepted:

    :param database: HHsearch database, equivalent to the ``-d`` flag.
    :param input: input alignment file, equivalent to the ``-i`` flag.
    :param evalue_cutoff: e-value cutoff, equivalent to the ``-e`` flag.
    :param output: output file, qquivalent to the ``-o`` flag.
    """

    program = "hhblits"
    flags = {
            "database":"d",
            "input": "i",
            "evalue_cutoff": "e",
            "output": "o"
            }

    def __init__(self, **flags):
        """Set up command-line flags for hhblits."""
        super().__init__(**flags)


