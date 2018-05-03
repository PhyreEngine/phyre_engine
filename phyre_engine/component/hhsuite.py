"""Components for running parts of hh-suite."""
import Bio.SeqUtils
import jmespath

from phyre_engine.component import Component
import phyre_engine.tools.hhsuite.tool as tools
import phyre_engine.tools.hhsuite.parser as parser
from phyre_engine.tools.external import ExternalTool
from phyre_engine.tools.jmespath import JMESExtensions

from enum import Enum
import os
import subprocess
import contextlib
import fileinput
import tempfile
import textwrap
import re
from phyre_engine.tools.template import Template
import math
from pathlib import Path, PurePath
import shutil

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
    HMM = "hhm"

class HHSuiteTool(Component):  #pylint: disable=abstract-method
    """
    Base class for HH-Suite tools.

    :param tuple(str, phyre_engine.tools.external.ExternalTool): Name and
        callable for generating command-line arguments.
    :param list[str] flags: List of flags without values to pass to tool.
    :param dict[str,str] options: Flags with values to pass to the tool.
    :param str bin_dir: Optional directory containing executable.
    :param str HHLIB: Optional HHLIB environment variable.
    :param QueryType input_type: Query type.
    :param str cache_dir: If all output files are already present in this
        directory, the tool will not be run. Set this to `None` to disable
        caching.
    """

    CONFIG_SECTION = "hhsuite"

    def __init__(self, tool, flags, options, bin_dir, HHLIB, input_type,
                 cache_dir="."):
        self.name, self.tool = tool
        self.flags = flags if flags is not None else []
        self.options = options if options is not None else {}
        self.bin_dir = bin_dir
        self.HHLIB = HHLIB
        self.input_type = QueryType(input_type)
        self.cache_dir = cache_dir

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


    def output_keys(self):
        """
        Parse the command-line parameters of the tool to set the keys on the
        command line pointing to output files. The following keys may be set:

        ``a3m``
            Set if the ``oa3m`` option is supplied.

        ``report``
            Set if the ``output`` (``-o``) option is supplied.

        ``atab``
            Set if the ``atab`` option is supplied.

        ``pairwise_fasta``
            Set if the ``Ofas`` option is supplied.
        """
        output = {}
        a3m = self._find_option(r"^-?oa3m$")
        report = self._find_option(r"^(?:-?o|output)$")
        atab = self._find_option(r"^-?atab$")
        ofas = self._find_option(r"^-?Ofas$")

        if a3m: output["a3m"] = a3m
        if report: output["report"] = report
        if atab: output["atab"] = atab
        if ofas: output["pairwise_fasta"] = ofas

        return output

    def cached(self):
        """
        Return `True` if caching is enabled and all output files already exist.
        """
        if self.cache_dir is None:
            return False

        for out_file in self.output_keys().values():
            if not Path(self.cache_dir, out_file).exists():
                return False

        return True

    def execute(self, data):
        """
        Execute tool, writing query sequence to temporary file if
        necessary.
        """
        if self.cached():
            return

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

            self.logger.debug("Running %s", cmd_line)
            subprocess.run(cmd_line, env=env, check=True)

class HHBlits(HHSuiteTool):
    """
    Run hhblits. This can be used to generate an MSA from a query sequence or
    to search a fold library given an MSA.

    :param str database: Path to an hhblits database.
    :param flags: Flags to pass to hhblits (e.g. ``all`` to pass the ``-all``
        option).
    :param options: Flag/option pairs to pass to hhblits (e.g. ``n=3`` to pass
        the option ``-n 3`` to hhblits).
    :param str bin_dir: Directory containing the hh-suite binaries.
    :param str HHLIB: Set the ``HHLIB`` environment variable to this value
        before calling hhblits.
    :param QueryType input_type: Input type.

    .. seealso::

        :py:data:`phyre_engine.tools.hhsuite.tool.hhblits`
            Command-line interface to ``hhblits``, including aliased
            parameters.
    """
    REMOVES = []

    @property
    def REQUIRED(self):
        if self.input_type == QueryType.SEQUENCE:
            return ["name", "sequence"]
        elif self.input_type == QueryType.A3M:
            return ["a3m"]
        elif self.input_type == QueryType.HMM:
            # sic: hhsuite uses "hhm" files for HMMs, so the "hhm" key in the
            # pipeline state points to the file containing the HMM.
            return ["hhm"]

    @property
    def ADDS(self):
        return list(self.output_keys())

    def __init__(self, database, flags=None, bin_dir=None, HHLIB=None,
                 input_type=QueryType.SEQUENCE, options=None,
                 cache_dir="."):

        if options is None:
            options = {}
        options["database"] = database

        super().__init__(
            ("hhblits", tools.hhblits),
            flags, options, bin_dir, HHLIB, input_type, cache_dir)


    def run(self, data, config=None, pipeline=None):
        """Build a sequence profile using hhblits."""
        super().execute(data)
        data.update(self.output_keys())
        return data

class HHSearch(HHSuiteTool):
    """
    Run hhsearch to match a query profile (or HMM) against a database of
    profiles.

    .. seealso::

        :py:class:`~.HHBlits`
            Uses the same class parameters.

        :py:data:`phyre_engine.tools.hhsuite.tool.hhsearch`
            Command-line interface to ``hhsearch``, including aliased
            parameters.
    """

    REMOVES = []

    @property
    def REQUIRED(self):
        if self.input_type == QueryType.A3M:
            return ["a3m"]
        elif self.input_type == QueryType.HMM:
            return ["hhm"]
        else:
            raise AttributeError(
                "Input type {} not supported".format(self.input_type))

    @property
    def ADDS(self):
        return list(self.output_keys())

    def __init__(
            self, database, flags=None, bin_dir=None, HHLIB=None,
            input_type=QueryType.A3M, options=None,
            cache_dir="."):

        if options is None:
            options = {}
        options["database"] = database
        super().__init__(
            ("hhsearch", tools.hhsearch),
            flags, options, bin_dir, HHLIB, input_type, cache_dir)

    def run(self, data, config=None, pipeline=None):
        """Search a profile database using hhsearch."""
        super().execute(data)
        data.update(self.output_keys())
        return data

class HHMake(HHSuiteTool):
    """
    Run hhmake to convert an ``a3m`` file into an ``hhm`` file.

    .. seealso::

        :py:class:`~.HHBlits`
            Uses the same class parameters.

        :py:data:`phyre_engine.tools.hhsuite.tool.hhsearch`
            Command-line interface to ``hhsearch``, including aliased
            parameters.
    """

    REQUIRED = ["a3m"]
    REMOVES = []
    ADDS = ["hhm"]

    def __init__(self, flags=None, bin_dir=None, HHLIB=None, options=None,
                 cache_dir="."):
        super().__init__(
            ("hhmake", tools.hhmake),
            flags, options, bin_dir, HHLIB, QueryType.A3M, cache_dir)


    def run(self, data, config=None, pipeline=None):
        """Build an hhm file from an a3m file."""
        super().execute(data)
        hhm_file = self._find_option(r"^(?:-?o|output)$")
        data["hhm"] = hhm_file
        return data

class CSTranslate(HHSuiteTool):
    """
    Run cstranslate to build column state sequences (cs219) for an a3m file.

    .. seealso::

        :py:class:`~.HHBlits`
            Uses the same class parameters.

        :py:data:`phyre_engine.tools.hhsuite.tool.hhsearch`
            Command-line interface to ``hhsearch``, including aliased
            parameters.
    """

    REQUIRED = ["a3m"]
    REMOVES = []
    ADDS = ["cs219"]

    def __init__(self, flags=None, bin_dir=None, HHLIB=None, options=None,
                 cache_dir="."):
        super().__init__(
            ("cstranslate", tools.cstranslate),
            flags, options, bin_dir, HHLIB, QueryType.A3M, cache_dir)

    @contextlib.contextmanager
    def set_input(self, data):
        """Sets --infile param for cstranslate."""
        options = {"infile": data["a3m"]}
        options.update(self.options)
        yield options

    def output_keys(self):
        """Output keys, extracted from outfile parameter."""
        cs219_file = self._find_option(r"^(?:-?o|outfile)$")
        return {"cs219": cs219_file}

    def run(self, data, config=None, pipeline=None):
        """Build 219-state sequence representing sequence profile."""
        super().execute(data)
        data.update(self.output_keys())
        return data


class AddPsipred(Component):
    """
    Add predicted secondary structure to MSAs.

    This component is not very smart, and simply calls the ``addss.pl`` script
    included in hh-suite. This means that paths to PSIPRED and legacy BLAST must
    be set up in ``HHPaths.pm``.
    """

    REQUIRED = ["a3m"]
    ADDS = []
    REMOVES = []
    CONFIG_SECTION = "hhsuite"

    @classmethod
    def config(cls, params, pipeline_config):
        if "HHLIB" in params:
            return params

        if "hhsuite" in pipeline_config:
            if "HHLIB" in pipeline_config["hhsuite"]:
                params["HHLIB"] = pipeline_config["hhsuite"]["HHLIB"]
        return params

    def __init__(self, HHLIB=None):
        if HHLIB is None and "HHLIB" not in os.environ:
            raise ValueError(
                "HHLIB not set as parameter or environment variable.")
        self.HHLIB = HHLIB

        hhlib = HHLIB if HHLIB is not None else os.environ["HHLIB"]
        self.addss = str(Path(hhlib, "scripts/addss.pl"))

    def run(self, data, config=None, pipeline=None):
        a3m = self.get_vals(data)
        # Skip if an ss_pred line is already present
        with open(a3m, "r") as a3m_in:
            for line in a3m_in:
                if line.startswith(">ss_pred"):
                    return data

        cmd_line = [self.addss, "-i", a3m]
        self.logger.debug("Running %s", cmd_line)
        tools.run(cmd_line, check=True, HHLIB=self.HHLIB)
        return data


class AddDssp(Component):
    """
    Add secondary structure information to MSAs produced by hhblits using
    `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_.

    This sets the ``>ss_dssp`` and ``>aa_dssp`` fields in the MSA pointed to by
    the ``a3m`` key. If those lines are already present in the MSA nothing is
    done.

    This component requires the ``secondary_structure`` key of the pipeline
    state as the source of secondary structure data. It also requires the
    ``template_obj`` field to point to a
    :py:class:`phyre_engine.tools.template.Template` object so that the
    secondary structure can be matched to the canonical sequence.
    """

    REQUIRED = ["a3m", "secondary_structure", "template_obj"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Add ``>ss_dssp`` and ``>aa_dssp`` fields."""
        a3m, sec_struc, template = self.get_vals(data)
        if self.contains_dssp(a3m):
            return data

        sec_struc = sec_struc["dssp"]

        # Index secondary structure states by residue ID
        ss_struc_dict = {ss["res_id"]: ss for ss in sec_struc}

        # Default everything to a gap, so we have sequences of the correct
        # length even if DSSP is missing an assignment.
        ss_dssp = ["-"] * len(template.canonical_indices)
        aa_dssp = ["-"] * len(template.canonical_indices)

        # For each residue in the canonical sequence (which is what should be in
        # the a3m file), get the SS state and set the corresponding char in the
        # sequences.
        for i, canonical_id in enumerate(template.canonical_indices):
            if canonical_id in ss_struc_dict:
                ss_state = ss_struc_dict[canonical_id]["assigned"]
                residue = template.chain[canonical_id]
                ss_dssp[i] = ss_state
                aa_dssp[i] = Bio.SeqUtils.seq1(residue.get_resname())

        aa_dssp = "".join(aa_dssp)
        ss_dssp = "".join(ss_dssp)
        self.update_a3m(a3m, ss_dssp, aa_dssp)
        return data

    @staticmethod
    def contains_dssp(a3m):
        """
        A `True` if both the ``ss_dssp`` and ``aa_dssp`` lines are present in
        the `a3m` file.
        """
        ss_dssp = False
        aa_dssp = False
        with open(a3m, "r") as a3m_in:
            for line in a3m_in:
                if line.startswith(">ss_dssp"):
                    ss_dssp = True
                elif line.startswith(">aa_dssp"):
                    aa_dssp = True

                if aa_dssp and ss_dssp:
                    return True

        return ss_dssp and aa_dssp

    def update_a3m(self, a3m, ss_dssp, aa_dssp):
        """
        If `a3m` does not contain ``ss_dssp`` or ``aa_dssp`` lines, add them.
        """
        # First, check if the "ss_dssp" or "aa_dssp" lines are present
        with open(a3m, "r") as a3m_in:
            a3m_lines = a3m_in.readlines()

        with open(a3m, "w") as a3m_out:
            if not any(ln.startswith(">aa_dssp") for ln in a3m_lines):
                a3m_out.write(">aa_dssp\n{}\n".format(aa_dssp))
            if not any(ln.startswith(">ss_dssp") for ln in a3m_lines):
                a3m_out.write(">ss_dssp\n{}\n".format(ss_dssp))
            a3m_out.writelines(a3m_lines)


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

    See :py:class:`phyre_engine.tools.hhsuite.parser.Report.Hit` for a full
    list of fields added to each element of the ``templates`` list.

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

        "sequence_alignments": {
            "ss_pred":  "-----CCHH...",
            "ss_conf":  "-----8899....",
            "ss_dssp":  "-----CCTH...",  #Optionally
            "template": "-----IDLE...",
            "confidence": "-----8899..."
        }

    If no alignment was generated for a hit or that
    alignment could not be parsed, ``sequence_alignments`` will be empty.

    If the ``alignment`` key is present in the template, the residue-pair
    confidences will be used to generate a string of confidence values.

    To ignore certain sequences, specify them in the ``ignore`` key.
    """
    REQUIRED = ["templates", "sequence", "name", "pairwise_fasta"]
    ADDS = []
    REMOVES = []

    def __init__(self, ignore=()):
        self.ignore = ignore

    @staticmethod
    def _add_confidence(query_seq, template):
        # First, index the alignment by query index. Offset by 1 to
        # convert between residue numbering and array numbering
        aln_map = {aln.i - 1: aln for aln in template["alignment"]}

        # Start with gaps everywhere
        confidence_str = ["-"] * len(query_seq)

        # query_index gives the residue index, and query_str_index
        # gives the index in the query string.
        query_index = 0
        for query_str_index, residue in enumerate(query_seq):
            if residue == "-":
                continue
            if query_index in aln_map:
                residue_pair = aln_map[query_index]
                conf_value = math.floor(residue_pair.probab * 10)
                if conf_value >= 10:
                    conf_value = 9
                confidence_str[query_str_index] = str(conf_value)
            query_index += 1
        template["sequence_alignments"]["confidence"] = "".join(confidence_str)

    def run(self, data, config=None, pipeline=None):
        templates, sequence, qname, fasta = self.get_vals(data)
        parsed_alignments = parser.Fasta(fasta, qname)

        for template, hit in zip(templates, parsed_alignments.hits):
            if hit is None:
                continue
            template["sequence_alignments"] = {}

            try:
                query_seq = hit["query"]["sequence"]
                query_seq = (
                    sequence[:template["query_range"].start - 1]
                    + query_seq
                    + sequence[template["query_range"].stop:])
                template["sequence_alignments"]["query"] = query_seq

                # Write confidence string
                if "alignment" in template:
                    self._add_confidence(query_seq, template)

                for seq_name, seq in hit["template"].items():
                    if seq_name in self.ignore:
                        continue

                    # Pad query and template according to start/end of query
                    # alignment
                    seq = (
                        "-" * (template["query_range"].start - 1)
                        + seq
                        + "-" * (len(sequence) - template["query_range"].stop))
                    template["sequence_alignments"][seq_name] = seq
            except KeyError as key_err:
                self.logger.warning("Error reading sequences for %s",
                                    template, exc_info=key_err)
        return data

class A3MSSParser(Component):
    """
    Parse secondary structure information from an a3m file. Instead of adding a
    ``secondary_structure`` element, like the components in
    :py:mod:`phyre_engine.component.secstruc`, this will add the
    ``secondary_structure_sequence`` key, which will contain the secondary
    structure lines from the parsed a3m file, indexed by the sequence name.

    >>> from phyre_engine.component.hhsuite import ParseA3MSS
    >>> parser = ParseA3MSS()
    >>> results = parser.parse({"a3m": "foo.a3m"})
    >>> results["secondary_structure_sequence"]
    {
        "ss_conf": "88889999...",
        "ss_pred": "CCCCHHHH...",
        "ss_dssp": "CCCTHHHH..."
    }

    Only the lines included in the a3m are added, so be sure to check which
    sequences are present.
    """
    REQUIRED = ["a3m"]
    ADDS = ["secondary_structure_sequence"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Parse secondary structure from a3m file."""
        a3m_file = self.get_vals(data)
        secstruc = {}
        with open(a3m_file, "r") as a3m_in:
            for line in a3m_in:
                if line.startswith(">"):
                    seq_id = line.split(" ")[0][1:].strip()
                    if seq_id.startswith("ss_"):
                        sequence = a3m_in.readline().strip()
                        secstruc[seq_id] = sequence
        data["secondary_structure_sequence"] = secstruc
        return data

class AlignmentToFasta(Component):
    """
    Convert an alignment--that is, a list of *i, j* pairs--into FASTA format.

    .. note::

        The alignments generated by this component are quite simplistic. They
        always include the entire query sequence as consecutive (non-gapped)
        residues, and map the template sequence onto that.

    This component reads the following keys:

    ``query_sequence``
        The sequence of the query protein.

    ``template_obj``
        A :py:class:`~phyre_engine.tools.template.Template` object from which
        the template sequence is extracted.

    ``alignment``
        The actual alignment, expressed as a tuple of ``(i, j)`` indices
        that indicate an alignment between query residue *i* and template
        residue *j*. See :py:class:`.TabularParser` for an example of a
        component that produces an ``alignment`` key.

    The FASTA alignment will be added to the ``fasta_alignment`` field of the
    given query. It will also be written to `file_name` (formatted with
    elements of the pipeline state), if `file_name` is not `None`.

    :param str q_name: Name of the query sequence written to the FASTA file.
        This is a Python format string that is passed each element in the
        pipeline state (e.g. ``{name}.fasta``). The default is ``Query``.

    :param str t_name: Like `q_name`, but for the template sequence. The
        default is ``Template``.

    :param str file_name: Template used to generate the file name to be
        written, or `None` to disable writing.
    """
    # TODO: This should operate on a single hit
    REQUIRED = ["query_sequence", "alignment", "template_obj"]
    ADDS = ["fasta_alignment"]
    REMOVES = []

    def __init__(self, q_name="Query", t_name="Template",
                 file_name="{rank:02d}-{PDB}_{chain}.fasta"):
        self.q_name = q_name
        self.t_name = t_name
        self.file_name = file_name

    @staticmethod
    def build_template_aln(query_seq, alignment, template):
        """Build template alignment string."""
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
        return self.q_name.format(**data)

    def _template_name(self, data):
        return self.t_name.format(**data)

    def fasta(self, data, hit_seq):
        return textwrap.dedent("""\
        >{query_name}
        {query_seq}
        >{template_name}
        {template_seq}
        """).format(
            query_name=self._query_name(data),
            template_name=self._template_name(data),
            query_seq=data["query_sequence"],
            template_seq=hit_seq)

    @staticmethod
    def _write_file(file_name, file_contents):
        """Write FASTA sequence to file."""
        with open(file_name, "w") as fasta_out:
            fasta_out.write(file_contents)

    def run(self, data, config=None, pipeline=None):
        """Convert alignments to FASTA."""

        query_sequence, alignment, template = self.get_vals(data)

        template_seq = self.build_template_aln(
            query_sequence,
            alignment,
            template)
        data["fasta_alignment"] = self.fasta(data, template_seq)

        if self.file_name is not None:
            self._write_file(
                self.file_name.format(**data),
                data["fasta_alignment"])

        return data

class PSSM(Component):
    """
    Using PSI-BLAST from the (legacy) NCBI BLAST toolkit, generate a PSSM.

    This component builds a PSSM in the same way as the ``addss.pl`` script
    that is packaged with hh-suite:

    1. If the query sequence contains gaps, calculate a consensus sequence
       using ``hhconsensus`` and use that as the query sequence for PSI-BLAST.

    2. Filter the query MSA down to a :math:`N_\text{eff}` of 7 using
       ``hhfilter``.

    3. Reformat the filtered query sequence into PSI-BLAST format. Remove the
       pseudo-sequence such as secondary structure in this step.

    4. Use ``blastpgp`` to generate a PSSM using a dummy database (bundled with
       hh-suite).

    This component requires the fields ``name`` and ``a3m`` in the pipeline
    state. The ``name`` is used to ensure that the correct query sequence is
    extracted from the MSA.

    The field ``pssm`` will be added to the pipeline state. The field will be
    a dict, the keys of which indicate different formats of PSSM and the
    values of which will point to the file containing the PSSM. The following
    formats will be added:

    ``mtx``
        IMPALA-format PSSM, generated from the PSI-BLAST checkpoint using
        ``makemat``

    ``chk``
        Not technically a PSSM, this is the ``blastpgp`` checkpoint file,
        generated using the ``-B`` option.

    ``ascii``
        The ASCII PSSM written by ``blastpgp`` when the ``-Q`` option is
        supplied.
    """

    # TODO: This component needs the "run" method refactoring and needs a unit
    # test. Fix this once you're no longer on a strict deadline!
    REQUIRED = ["name", "a3m"]
    ADDS = ["pssm"]
    REMOVES = []

    def __init__(self, hhsuite_dir, HHLIB, blast_dir):
        self.hhsuite_dir = hhsuite_dir
        self.HHLIB = HHLIB
        self.blast_dir = blast_dir

    def run(self, data, config=None, pipeline=None):
        """Generating PSSM from MSA."""
        hhconsensus = ExternalTool({
                "input": "i",
                "seqfile": "s",
                "verbose": "v",
            }, long_prefix="-")
        hhfilter = ExternalTool({
            "input": "i",
            "oa3m": "o",
            "verbose": "v"}, long_prefix="-")
        reformat = ExternalTool({"no_lower": "r"}, long_prefix="-")
        blastpgp = ExternalTool({
                "database": "d",
                "input": "i",
                "iterations": "j",
                "db_seq_alns": "b",
                "evalue_threshold": "h",
                "input_alignment": "B",
                "checkpoint": "C",
                "ascii_pssm": "Q"})
        makemat = ExternalTool({"profile_db": "P"})

        try:
            tmpdir = tempfile.mkdtemp("-pssm", "phyreengine-")
            name, a3m = self.get_vals(data)
            query_seq = self.read_query_seq(a3m, name)

            tmp_a3m = Path(tmpdir, "msa.a3m")
            tmp_seq = Path(tmpdir, "seq.fasta")
            tmp_psi = Path(tmpdir, "msa.psi")
            env = os.environ.copy()
            env["HHLIB"] = self.HHLIB

            if "-" in query_seq:
                # Generate consensus sequence
                command_line = hhconsensus(
                    executable=(self.hhsuite_dir, "hhconsensus"),
                    options={
                        "input": a3m,
                        "seqfile": tmp_seq,
                        "oa3m": tmp_a3m})
                subprocess.run(command_line, check=True, env=env)
            else:
                # If there are no gaps, just copy the query a3m and write the
                # query sequence.
                shutil.copy2(a3m, str(tmp_a3m))
                with tmp_seq.open("w") as tmp_seq_out:
                    tmp_seq_out.write(">{name}\n{query}\n".format(
                        name=name, query=query_seq))

            # Filter query a3m to desired diversity
            command_line = hhfilter(
                executable=(self.hhsuite_dir, "hhfilter"),
                options={
                    "neff": 7,
                    "input": tmp_a3m,
                    "oa3m": tmp_a3m})
            subprocess.run(command_line, check=True, env=env)

            # Reformat to PSI-BLAST format
            command_line = reformat(
                (str(Path(self.HHLIB, "scripts")), "reformat.pl"),
                positional=["a3m", "psi", tmp_a3m, tmp_psi],
                flags=["no_lower", "noss"])
            subprocess.run(command_line, check=True, env=env)

            # Generate PSSM using blastpgp
            chk_file = "profile.chk"
            mtx_file = "profile.mtx"
            pssm_file = "profile.pssm"

            dummy_db = Path(self.HHLIB, "data/do_not_delete")
            command_line = blastpgp(
                executable=(self.blast_dir, "blastpgp"),
                options={
                    "db_seq_alns": 1,
                    "iterations": 1,
                    "evalue_threshold": 0.001,
                    "database": dummy_db,
                    "input": tmp_seq,
                    "input_alignment": tmp_psi,
                    "checkpoint": chk_file,
                    "ascii_pssm": pssm_file})
            subprocess.run(command_line, check=True)

            # Build mtx file using makemat

            # First build the profile "databases" for use with makemat. These
            # are just two files with the same prefix and the suffixes ".pn"
            # and ".sn", which contain a list of checkpoint files and the
            # corresponding list of sequences. File names are resolve relative
            # to the directory containing the *.sn and *.pn files, so we will
            # symlink the checkpoint.
            tmp_sn_file = Path(tmpdir, "makemat.sn")
            tmp_pn_file = Path(tmpdir, "makemat.pn")

            with tmp_sn_file.open("w") as sn_out:
                print(str(Path(tmp_seq.name)), file=sn_out)
            with tmp_pn_file.open("w") as pn_out:
                chk_path = Path(chk_file)
                chk_link = Path(tmpdir, "makemat.chk")
                chk_link.symlink_to(chk_path.resolve())
                print(str(chk_link.name), file=pn_out)

            command_line = makemat(
                executable=(self.blast_dir, "makemat"),
                options={"profile_db": str(Path(tmpdir, "makemat"))})
            subprocess.run(command_line, check=True)
            shutil.copy2(str(Path(tmpdir, "makemat.mtx")), mtx_file)

            data["pssm"] = {
                "mtx": mtx_file,
                "chk": chk_file,
                "ascii": pssm_file}
            return data

        finally:
            shutil.rmtree(tmpdir)

    def read_query_seq(self, a3m, name):
        """Read the query sequence from an a3m file."""
        query = []

        record_flag = False
        with open(a3m, "r") as a3m_in:
            for line in a3m_in:
                if line.startswith(">"):
                    if line.startswith(">" + name):
                        record_flag = True
                    elif record_flag:
                        break
                elif record_flag:
                    query.append(line.strip())
        if not query:
            raise ValueError(
                "Could not read query sequence '{}' from '{}'.".format(
                    name, a3m))
        return "".join(query)


class FFDatabaseUnlink(Component):
    """
    Unlink (remove from the index) entries from an ffindex database.

    Similar to :py:class:`.BuildDatabase`, this component will loop over each
    element returned by the JMESPath query `select_expr`. For each element,
    it remove the item given by the ``name`` field from all ffindex files in
    the named database.

    :param str db_prefix: Prefix for database. The databases used by hhblits
        consist of multiple files, named like
        ``<prefix>_{a3m,hhm,cs219}.ff{index,data}``. This may be supplied as
        an absolute or relative path. Relative paths are evaulatued relative
        to the current working directory.

    :param str bin_dir: Optional directory containing the ``ffindex_build``
        executable.

    :param str select_expr: JMESPath expression giving the list of templates
        in the pipeline state.
    """
    REQUIRED = []
    REMOVES = []
    ADDS = []

    @classmethod
    def config(cls, params, config):
        return config.extract({
            "foldlib": ["db_prefix", "overwrite"],
            "hhsuite": ["bin_dir"],
        }).merge_params(params)

    def __init__(self, db_prefix, bin_dir=None, select_expr="templates"):
        self.db_prefix = db_prefix
        self.bin_dir = bin_dir
        self.select_expr = select_expr

    def run(self, data, config=None, pipeline=None):
        """Collect and index the files that form an hhsuite database."""
        jmes_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        templates = jmespath.search(self.select_expr, data, jmes_opts)

        with tempfile.NamedTemporaryFile("w") as file_list:
            # Write IDs to remove to temp file
            for template in templates:
                print(template["name"], file=file_list)
            file_list.flush()

            db_types = ["a3m", "hhm", "cs219"]
            db_prefix = Path(self.db_prefix)

            for file_type in db_types:
                ffindex = Path("{}_{}.ffindex".format(db_prefix, file_type))
                cmd_line = tools.ffindex_modify(
                    (self.bin_dir, "ffindex_modify"),
                    options={"file_list": file_list.name},
                    flags=["sort", "unlink"],
                    positional=[ffindex])
                self.logger.debug("Running command %s", cmd_line)
                tools.run(cmd_line, check=True)
        return data


class BuildDatabase(Component):
    """
    Build ffindex/ffdata files for an hhsuite database.

    This component will iterate over each template in the list specified by the
    JMESPath expression `select_expr` (by default, the ``templates`` list).
    Each template must have the keys ``sequence``, ``a3m`` (from e.g.
    :py:class:`~.HHBlits`), ``hhm`` (from e.g. :py:class:`~.HHMake`) and
    ``cs219`` (from e.g. :py:class:`~.CSTranslate`) defined.

    The procedure used to build the database is essentially the same as given
    in section 3.5 of the `hh-suite user manual
    <https://github.com/soedinglab/hh-suite/raw/master/hhsuite-userguide.pdf>`_.
    This component assumes that the relevant data files have already been
    generated by, for example, :py:class:`.HHBlits`, :py:class:`.HHMake`, and
    :py:class:`.CSTranslate`. It will then call ``ffindex_build`` to build the
    ``ffindex`` and ``ffdata`` files required by hhsuite. Files are sorted by
    sequence length as indicated in the user manual, and the file names given
    in the ``ffindex`` files are cleaned up to remove the file suffixes for the
    sake of consistency.

    :param str db_prefix: Prefix for database. The databases used by hhblits
        consist of multiple files, named like
        ``<prefix>_{a3m,hhm,cs219}.ff{index,data}``. This may be supplied as
        an absolute or relative path. Relative paths are evaulatued relative
        to the current working directory.

    :param str bin_dir: Optional directory containing the ``ffindex_build``
        executable.

    :param bool overwrite: If ``True``, delete existing database files.
        Otherwise, ``ffindex_build`` may be called on existing files.

    :param str select_expr: JMESPath expression giving the list of templates
        in the pipeline state.
    """

    REQUIRED = []
    REMOVES = []
    ADDS = ["database"]

    def __init__(self, db_prefix, bin_dir=None, overwrite=False,
                 select_expr="templates"):
        self.db_prefix = db_prefix
        self.bin_dir = bin_dir
        self.overwrite = overwrite
        self.select_expr = select_expr

    @classmethod
    def config(cls, params, config):
        return config.extract({
            "foldlib": ["db_prefix", "overwrite"],
            "hhsuite": ["bin_dir"],
        }).merge_params(params)

    def run(self, data, config=None, pipeline=None):
        """Collect and index the files that form an hhsuite database."""
        jmes_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        templates = jmespath.search(self.select_expr, data, jmes_opts)

        # First, sort templates by sequence length.
        templates.sort(key=lambda t: len(t["sequence"]))

        # Make database directory if it doesn't exist.
        Path(self.db_prefix).parent.mkdir(parents=True, exist_ok=True)

        # Collect a3m/hhm/cs219 files into ffindex/ffdata databases
        to_collect = ["a3m", "hhm", "cs219"]
        ff_dbs = {}
        db_prefix = Path(self.db_prefix)

        for file_type in to_collect:
            db_name = Path("{}_{}".format(str(db_prefix), file_type))
            ffindex = Path("{}.ffindex".format(str(db_name)))
            ffdata = Path("{}.ffdata".format(str(db_name)))

            if self.overwrite:
                if ffindex.exists():
                    ffindex.unlink()
                if ffdata.exists():
                    ffdata.unlink()

            with tempfile.NamedTemporaryFile("w") as index:
                # Write all files of file_type `file_type` to a temp file
                for template in templates:
                    print(template[file_type], file=index)
                index.flush()

                # Run ffindex_build using the the temp file as the list of files
                # to include in the DB.
                cmd_line = tools.ffindex_build(
                    (self.bin_dir, "ffindex_build"),
                    positional=[ffdata, ffindex],
                    flags=["sort", "append"],
                    options={"file_list": index.name})
                self.logger.debug("Running command %s", cmd_line)
                tools.run(cmd_line, check=True)
                ff_dbs[file_type] = db_name

        # Cut useless information from the indices of each file.
        for ff_db in ff_dbs.values():
            self._trim_index_names(ff_db)

        data["database"] = str(db_prefix)
        return data


    def _trim_index_names(self, db):
        """Replace names of components in an ffindex with their stem.

        When we generate an index using ``ffindex_build``, the names of each
        element in that database are simply taken from the filenames used as
        input. Since we are trying to build a valid database for hhblits and
        hhsearch, each name in the index must be consistent across database. To
        do that, we use the name of each template as the identifier.
        """

        index = "{}.ffindex".format(db)
        with fileinput.input(index, inplace=True) as fh:
            for line in fh:
                fields = line.split("\t")
                fields[0] = PurePath(fields[0]).stem
                print("\t".join(fields), end="")
