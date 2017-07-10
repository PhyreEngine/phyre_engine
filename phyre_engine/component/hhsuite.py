"""Components for running parts of hh-suite."""
from phyre_engine.component import Component
import phyre_engine.tools.hhsuite.tool as tools
import phyre_engine.tools.hhsuite.parser as parser

import Bio.SeqIO
from enum import Enum
import os
import subprocess
import contextlib
import tempfile
import re

class QueryType(Enum):
    """
    The different query types that may be read by classes in this module.

    An enumeration value from this object may be passed as the ``input_type``
    parameter to some of the classes in this module. The corresponding key from
    the pipeline state is then treated as the query (the ``-i`` flag for most
    hh-suite tools).
    """

    #: Pass the file pointed to by the ``input`` key of the pipeline state
    #: directly to the ``-i`` key.
    FILE = "input"

    #: Write the :py:class:`Bio.Seq.SeqRecord` object pointed to by the
    #: ``sequence`` key of the pipeline state to a temporary file and use that
    #: as the query sequence.
    SEQUENCE = "sequence"

    #: Use the multiple sequence alignment saved in the file named by the
    #: ``msa`` key as input.
    MSA = "msa"

    #: Use the HMM pointed in the file named by the ``hmm`` key as input.
    HMM = "hmm"

class HHSuiteTool(Component):  #pylint: disable=abstract-method
    """
    Base class for HH-Suite tools.

    :param tuple(str, phyre_engine.tools.external.ExternalTool): Name and
        callable for generating command-line arguments.
    :param list[str] flags: List of flags without values to pass to tool.
    :param dict[str, str] options: Flags with values to pass to the tool.
    :param str bin_dir: Optional directory containing executable.
    :param str HHLIB: Optional HHLIB environment variable.
    :param QueryType input_type: Query type.
    """
    def __init__(self, tool, flags, options, bin_dir, HHLIB, input_type):
        self.name, self.tool = tool
        self.flags = flags
        self.options = options
        self.bin_dir = bin_dir
        self.HHLIB = HHLIB
        self.input_type = QueryType(input_type)

    @contextlib.contextmanager
    def set_input(self, data):
        """
        Context manager used to set the ``input`` (``-i``) flag for hh-suite
        tools. This will write a temporary file if necessary.
        """
        temp_query = None
        try:
            if "input" not in self.options and "i" not in self.options:
                options = self.options.copy()
                if self.input_type == QueryType.SEQUENCE:
                    # Write temporary file and set option
                    temp_query = tempfile.NamedTemporaryFile(
                        "w", prefix="query-", suffix=".fasta")
                    Bio.SeqIO.write(
                        data[self.input_type.value], temp_query, "fasta")
                    temp_query.flush()
                    options["input"] = temp_query.name
                    yield options
                else:
                    options["input"] = str(data[self.input_type.value])
                    yield options
            else:
                yield self.options
        finally:
            if temp_query is not None:
                temp_query.close()

    def _find_option(self, regex):
        """
        Find an option in ``self.options`` matching the regex. Returns the value
        of the matching option, or ``None`` if the option was not found.
        """
        compiled_regex = re.compile(regex)
        for key, value in self.options.items():
            if compiled_regex.search(key):
                return value
        return None


    def set_output_key(self, data):
        """
        Set keys in the pipeline state depending on the options supplied to
        hhblits. The following keys may be set:

        ``msa``
            Set if the ``oa3m`` option is supplied.

        ``report``
            Set if the ``output`` (``-o``) option is supplied.

        ``atab``
            Set if the ``atab`` option is supplied.
        """
        msa = self._find_option(r"^-?oa3m$")
        report = self._find_option(r"^(?:-?o|output)$")
        atab = self._find_option(r"^-?atab$")

        if msa: data["msa"] = msa
        if report: data["report"] = report
        if atab: data["atab"] = atab

        return data

    def execute(self, data):
        """
        Execute tool, writing query sequence to temporary file if
        necessary.
        """
        with self.set_input(data) as options:
            cmd_line = self.tool(
                (self.bin_dir, self.name),
                positional=None,
                flags=self.flags,
                options=options)

            # Set HHLIB environment variable if provided. Emulate default
            # behaviour by inheriting system environment.
            env = None
            if self.HHLIB is not None:
                env = os.environ
                env["HHLIB"] = str(self.HHLIB)

            subprocess.run(cmd_line, env=env, check=True)

class HHBlits(HHSuiteTool):
    """
    Run hhblits. This can be used to generate an MSA from a query sequence or
    to search a fold library given an MSA.

    :param str database: Path to an hhblits database.
    :param *args: Flags to pass to hhblits (e.g. ``all`` to pass the ``-all``
        option).
    :param **kwargs: Flag/option pairs to pass to hhblits (e.g. ``n=3`` to pass
        the option ``-n 3`` to hhblits).
    :param str bin_dir: Directory containing the hh-suite binaries.
    :param str HHLIB: Set the ``HHLIB`` environment variable to this value
        before calling hhblits.
    :param QueryType input_type: Input type.
    """

    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(
            self, database, *args, bin_dir=None, HHLIB=None,
            input_type=QueryType.SEQUENCE, **kwargs):

        kwargs["database"] = database
        super().__init__(
            ("hhblits", tools.hhblits),
            args, kwargs, bin_dir, HHLIB, input_type)


    def run(self, data, config=None, pipeline=None):
        """Build a sequence profile using hhblits."""
        super().execute(data)
        return self.set_output_key(data)

class HHSearch(HHSuiteTool):
    """
    Run hhsearch to match a query profile (or HMM) against a database of
    profiles.

    .. seealso::

        `~.HHBlits`
            For parameters.
    """

    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(
            self, database, *args, bin_dir=None, HHLIB=None,
            input_type=QueryType.MSA, **kwargs):

        kwargs["database"] = database
        super().__init__(
            ("hhsearch", tools.hhsearch),
            args, kwargs, bin_dir, HHLIB, input_type)

    def run(self, data, config=None, pipeline=None):
        """Search a profile database using hhsearch."""
        super().execute(data)
        return self.set_output_key(data)

class ReportParser(Component):
    """Parse hhsearch reports."""


    #: :param str hhsearch_atab: Location of a file produced by the ``-atab``
    #:      option of hhsearch.
    #: :param str hhsearch_report: Location of a report file produced by
    #:      hhsearch(``-o`` option).
    REQUIRED = ["hhsearch_atab", "hhsearch_report"]
    #: :param hits: Hits parsed from report.
    #: :type hits: List of :class:`phyre_engine.tools.hhsuite.parser.Hit`
    #:      objects.
    ADDS = ["hits"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Parse report file and atab file into an array of hits."""

        atab_file, report_file = self.get_vals(data)

        #Parse report files to get hits, then combine the hits so that we have
        #the summary information and alignments in the same hit.
        report = parser.Report(report_file)
        atab   = parser.Tabular(atab_file)
        hits = []
        for r, a in zip(report.hits, atab.hits):
            hit = parser.Hit()
            hit.info = r.info
            hit.aln  = a.aln
            hits.append(hit)

        data["hits"] = hits
        return data
