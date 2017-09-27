"""Components for running parts of hh-suite."""
from phyre_engine.component import Component
import phyre_engine.tools.hhsuite.tool as tools
import phyre_engine.tools.hhsuite.parser as parser

from enum import Enum
import os
import subprocess
import contextlib
import tempfile
import textwrap
import re
from phyre_engine.tools.template import Template

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

    #: Write the sequence pointed to by the
    #: ``sequence`` key of the pipeline state to a temporary file and use that
    #: as the query sequence.
    SEQUENCE = "sequence"

    #: Use the multiple sequence alignment saved in the file named by the
    #: ``a3m`` key as input.
    A3M = "a3m"

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

    CONFIG_SECTION = "hhsuite"

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
                    print(
                        ">{name}\n{sequence}".format(**data),
                        file=temp_query)
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

        ``a3m``
            Set if the ``oa3m`` option is supplied.

        ``report``
            Set if the ``output`` (``-o``) option is supplied.

        ``atab``
            Set if the ``atab`` option is supplied.

        ``pairwise_fasta``
            Set if the ``Ofas`` option is supplied.
        """
        a3m = self._find_option(r"^-?oa3m$")
        report = self._find_option(r"^(?:-?o|output)$")
        atab = self._find_option(r"^-?atab$")
        ofas = self._find_option(r"^-?Ofas$")

        if a3m: data["a3m"] = a3m
        if report: data["report"] = report
        if atab: data["atab"] = atab
        if ofas: data["pairwise_fasta"] = ofas

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
            input_type=QueryType.A3M, **kwargs):

        kwargs["database"] = database
        super().__init__(
            ("hhsearch", tools.hhsearch),
            args, kwargs, bin_dir, HHLIB, input_type)

    def run(self, data, config=None, pipeline=None):
        """Search a profile database using hhsearch."""
        super().execute(data)
        return self.set_output_key(data)

class ReportParser(Component):
    """
    Parse hhsearch/hhblits report files.

    This component will read a report file produced by the ``-o`` option of
    hhblits or hhsearch. Only hits in the summary section of the hhsearch report
    will be parsed, so you may need to set the ``-Z`` option of hhsearch if you
    want to increase the number of reported hits.

    Results are added to the ``templates`` list of the pipeline state.

    The key ``report`` must be present in the pipeline state. This key should
    contain the location of the report file to be parsed.

    .. warning::

        This component is considered a primary source of templates. If a
        ``templates`` key is already present in the pipeline state, it will be
        *overwritten* by this component.
    """

    REQUIRED = ["report"]
    ADDS = ["templates"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Parse report file into an array of templates."""
        report_file = self.get_vals(data)
        report = parser.Report(report_file)
        templates = []

        # For each hit returned by the parser, just get each score and convert
        # it to a template dict.
        for hit in report.hits:
            template = dict(hit._asdict())
            templates.append(template)
        data["templates"] = templates
        return data

class TabularParser(Component):
    """
    Parse tabular file generated by hhsearch/hhblits using the ``-atab`` option.

    This component reads a tabular file generated by hhsearch or hhblits and
    parses the alignment for each template. Tabular files may be generated using
    the ``-atab`` option of hhblits or hhsearch (version 3 and above), and
    contain alignment information for individual residues.

    The alignment will be stored in the ``alignment`` key of each element of the
    ``templates`` list. Alignments are represented as lists of residue pairs.
    Each residue pair will have an ``i`` and ``j`` element indicating the
    query-template mapping, as well as any extra fields such as pair confidence.

    .. warning::

        This component is a secondary source of template information. If a
        ``templates`` key is not already present in the pipeline state it will
        be created. Otherwise, this component assumes that the hits in the
        tabular file appear in the same order as the templates in the
        ``templates`` list. If a mismatch is detected, a :py:exc:`ValueError`
        will be raised.
    """

    REQUIRED = ["atab"]
    ADDS = ["templates"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Parse tabular file into an array of templates with alignments."""
        atab_file = self.get_vals(data)
        atab = parser.Tabular(atab_file)

        if "templates" in data:
            for existing_hit, atab_hit in zip(data["templates"], atab.hits):
                if existing_hit["name"] != atab_hit.name:
                    raise ValueError(
                        "Hit name '{}' does not match '{}'".format(
                            existing_hit["name"], atab_hit.name))
                existing_hit["alignment"] = atab_hit.alignment
        else:
            templates = []
            for hit in atab.hits:
                templates.append(hit._asdict())
            data["templates"] = templates
        return data

class FastaParser(Component):
    """
    Parse a pairwise FASTA file generated by hhsearch or hhblits using the
    ``-Ofas`` option.

    Alignments will be stored in the ``sequence_alignment`` mapping of each
    element in the ``templates`` list. The ``sequence_alignment`` mapping will
    contain the sequences returned by hhblits. If the residue confidences are
    known, they will be included in the ``confidence`` field. The template
    sequences will be in the ``template`` field. All sequences will be given
    relative to the query, i.e. with gaps added at the beginning and end if
    necessary.

    .. code-block:: python
        “sequence_alignments”: {
            “ss_pred”:  “-----CCHH...",
            “ss_conf”:  “-----8899....”,
            “ss_dssp”:  “-----CCTH...",  #Optionally
            “template”: “-----IDLE...”,
            “confidence”: “-----8899..."
        }

    If no alignment was generated for a hit or that
    alignment could not be parsed, ``sequence_alignments`` will be empty.

    """
    REQUIRED = ["templates", "sequence", "name", "pairwise_fasta"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        templates, sequence, qname, fasta = self.get_vals(data)
        parsed_alignments = parser.Fasta(fasta, qname)

        for template, hit in zip(templates, parsed_alignments.hits):
            if hit is None:
                continue

            template["sequence_alignments"] = {}
            for seq_name, seq in hit["template"].items():
            # Pad query and template according to start/end of query alignment
                seq = (
                    "-" * (template["query_range"].start - 1)
                    + seq
                    + "-" * (len(sequence) - template["query_range"].stop))
                template["sequence_alignments"][seq_name] = seq
        return data

class AlignmentToFasta(Component):
    """
    Convert an alignment--that is, a list of *i, j* pairs--into FASTA format.
    This component must know the query sequence, and so must access the pipeline
    state as a whole, rather than an individual hit. Each hit must contain a
    ``template`` key pointing at a template file.

    The FASTA alignment will be added to the ``fasta_alignment`` field of the
    given query.

    :param str list_key: The list containing each hit.
    :param str q_name: The query name will be extracted from this field. By
        default, the query name is simply "Query".
    :param str t_name: The template name will be taken from this field. By
        default, the template name is "Template".
    """
    REQUIRED = ["sequence"]
    ADDS = []
    REMOVES = []

    def __init__(self, list_key="templates", q_name=None, t_name=None):
        self.list_key = list_key
        self.q_name = q_name
        self.t_name = t_name


    @staticmethod
    def build_template_aln(query_seq, alignment, template):
        # Store the template residue index, indexed by query index
        aln_mapping = {}
        for res_pair in alignment:
            q_index, t_index = res_pair[0:2]
            # The -1s are to convert from seq to array numbering
            aln_mapping[q_index - 1] = t_index - 1

        # Generate alignment
        t_seq = ["-"] * len(query_seq)
        for q_index, _ in enumerate(query_seq):
            if q_index in aln_mapping:
                t_index = aln_mapping[q_index]
                t_seq[q_index] = template.canonical_seq[t_index]
        return "".join(t_seq)

    def _query_name(self, data):
        return data[self.q_name] if self.q_name is not None else "Query"

    def _template_name(self, hit):
        return hit[self.q_name] if self.t_name is not None else "Template"

    def fasta(self, data, hit, hit_seq):
        return textwrap.dedent("""\
        >{query_name}
        {query_seq}
        >{template_name}
        {template_seq}
        """).format(
            query_name=self._query_name(data),
            template_name=self._template_name(hit),
            query_seq=data["sequence"],
            template_seq=hit_seq)

    def run(self, data, config=None, pipeline=None):
        """Convert alignments to FASTA."""
        for hit in data[self.list_key]:
            # Skip hits that don't have the information we need.
            if "alignment" not in hit or "template" not in hit:
                continue

            # Read template. The template sequence is the canonical sequence
            template = Template.load(hit["template"])

            template_seq = self.build_template_aln(
                data["sequence"], hit["alignment"], template)
            hit["fasta_alignment"] = self.fasta(data, hit, template_seq)
        return data
