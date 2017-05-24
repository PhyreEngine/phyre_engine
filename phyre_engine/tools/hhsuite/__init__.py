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

    Extra positional arguments will be appended to the command line verbatim.

    >>> class HHSearchExample(HHSuiteTool)
    >>>     flags = {"database": "d", "input": "i"}
    >>>     def __init__(self, program="hhsearch", **args):
    >>>         super().__init__(program, args)
    >>> hhs = HHSearchExample(database="uniprot_2012_03",
    >>>        input="in.fasta", cpu=10, nocons=True)
    >>> hhs.command_line
    ["hhsearch", "-d" "uniprot_2012_03", "-i", "in.fasta",
            "-cpu", "10", "-nocons"]

    :ivar list command_line: Command line of the program, separated into
    individual arguments. The first element is the name of the program.
    """

    def __init__(self, program, *args, **flags):
        """Set up a tool with the given command-line flags."""
        self.program = program
        long_flags = flags
        self.command_line = [self.program]
        self._map_flags(long_flags)
        self.command_line.extend(args)

    def _map_flags(self, flags):
        """Map human-readable option names to command-line flags."""
        for long_flag, value in flags.items():
            # Convert flag if a mapping was supplied
            if long_flag in self.flags:
                flag = self.flags[long_flag]
            else:
                flag = long_flag

            # Guess appropriate number of dashes if they weren't specified.
            if not flag.startswith("-"):
                if len(flag) > 1:
                    flag = "--" + flag
                else:
                    flag = "-" + flag

            if type(value) == bool:
                #Bool flags do not take a value
                self.command_line.append(flag)
            else:
                #Non-bool flags take an argument
                self.command_line.append(flag)
                self.command_line.append(str(value))


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

    flags = {
            "database":"d",
            "input": "i",
            "iterations": "n",
            "evalue_cutoff": "e",
            "output": "o"
            }

    def __init__(self, program="hhblits", **flags):
        """Set up command-line flags for hhblits."""
        super().__init__(program, **flags)


class HHSearch(HHSuiteTool):
    """Wrapper for the hhsearch tool included within the hh-suite package.

    As well as the standard short-form flags of hhsearch (see ``hhsearch
    -help``), the following flags are accepted:

    :param database: HHsearch database, equivalent to the ``-d`` flag.
    :param input: input alignment file, equivalent to the ``-i`` flag.
    :param evalue_cutoff: e-value cutoff, equivalent to the ``-e`` flag.
    :param output: output file, qquivalent to the ``-o`` flag.
    """

    flags = {
            "database":"d",
            "input": "i",
            "evalue_cutoff": "e",
            "output": "o"
            }

    def __init__(self, program="hhsearch", **flags):
        """Set up command-line flags for hhblits."""
        super().__init__(program, **flags)

class HHMake(HHSuiteTool):
    """Wrapper for the hhmake tool included within the hh-suite package.

    All options to ``hhmake`` may be included in the constructor. The input file
    (the ``-i`` option) is mandatory. The following aliases are accepted as
    well as the standard options.

    :param str input: Query alignment file.
    :param str output: HMM file to be written (``-o`` flag)
    :param str append: HMM file to be appended (``-a`` flag)
    :param int verbose: Increase verbosity from 0 (silent) to 2 (verbose). This
        is the ``-v`` flag.
    """

    flags = {
            "input": "i",
            "output": "o",
            "append": "a",
            "verbose": "v",
            }

    def __init__(self, input, program="hhmake", **flags):
        """Set up command-line flags for hhmake."""
        super().__init__(program, input=input, **flags)

class CSTranslate(HHSuiteTool):
    """
    Wrapper around hhsuite's cstranslate.

    The following options are available. Hyphens may be replaced with
    underscores to make valid python identifiers.

    .. code-block:: none
        infile <file>            Input file with alignment or sequence
        outfile <file>           Output file for generated abstract state sequence (def: <infile>.as)
        append <file>            Append generated abstract state sequence to this file
        informat prf|seq|fas|... Input format: prf, seq, fas, a2m, or a3m (def=auto)
        outformat seq|prf        Outformat: abstract state sequence or profile (def=seq)
        match-assign [0:100]     Make all FASTA columns with less than X% gaps match columns
                                 (def: make columns with residue in first sequence match columns)
        alphabet <file>          Abstract state alphabet consisting of exactly 219 states (def=off)
        context-data <file>      Add context-specific pseudocounts using given context-data (def=off)
        pc-admix [0,1]           Pseudocount admix for context-specific pseudocounts (def=0.90)
        pc-ali [0,inf[           Constant in pseudocount calculation for alignments (def=12.0)
        weight [0,inf[           Weight of abstract state column in emission calculation (def=1000.00)
        binary                   Write binary instead of character sequence (def=off)
        ffindex                  Read from -i <ffindex>, write to -o <ffindex>; enables openmp if possible (def=off)

    """

    flags = {
        "match_assign": "match-assign",
        "pc_admix": "pc-admix",
        "pc_ali": "pc-ali",
    }

    def __init__(self, program="cstranslate", **flags):
        """Set up command line for cstranslate."""
        super().__init__(program, **flags)

class FFIndexBuild(HHSuiteTool):
    """
    Wrapper around ffindex_build.

    :param str data: Output data file.

    :param str index: Output index file.

    :param bool append: (``-a``) Append files/indexes, also needed for sorting
        an already existing ffindex.

    :param str data_file_2: (``-d``) A second ffindex data file for
        inserting/appending.

    :param str index_file_2: (``-i``) A second ffindex index file for
        inserting/appending.

    :param str file_list: (``-f``) File containing a list of file names, one per
        line.

    :param str sort_file: (``-s``)  Sort index file, so that the index can queried.
        Another append operations can be done without sorting.
    """

    flags = {
        "append": "a",
        "data_file_2": "d",
        "index_file_2": "i",
        "file_list": "f",
        "sort_file": "s",
    }

    def __init__(self, data, index, program="ffindex_build", **flags):
        """Set up ffindex_build command line."""
        super().__init__(program, data, index, **flags)
