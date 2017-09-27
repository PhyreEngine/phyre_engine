"""Parsers for hhsuite files."""
from enum import Enum
import collections
import io
import logging
import re
from collections import namedtuple
import Bio.AlignIO
from phyre_engine.tools.util import Stream


log = lambda: logging.getLogger(__name__)

def _parse_range(range_str):
    """Parse a range string (``x-y``) into a ``range`` object."""
    start, end = range_str.split("-")
    return range(int(start), int(end))

def _parse_num_seqs(num_seqs_str):
    """
    Parse the No_of_seqs header field into a 2-tuple containing filtered and
    total number of sequences.
    """
    num_filt, num_in = num_seqs_str.split(" out of ")
    return (int(num_filt), int(num_in))

class Report:
    """Parse an hhsearch report file into a list of hits.

    This class will only parse data on a per-hit basis; it will not parse the
    actual alignments from the report file. To get the alignments, you should
    generate a tabular file using the ``-atab`` option of hhblits/hhsearch and
    parse the alignment data using the :class:`Tabular` parser.

    :param str file: Path of the report file to be parsed.
    :ivar hits: List of hits.
    :vartype hits: :py:class:`.Hit`
    :ivar summary: Header information parsed from the report file.
    :vartype summary: :py:class:`.Summary`
    """

    #: Represents a single hit parsed from an hhsuite report file. For detailed
    #: information, see section 5.1 of the hh-suite user guide.
    #:
    #: :param int rank: Rank in the hit list of this hit.
    #: :param str name: Name of this hit.
    #: :param float prob: True positive probability.
    #: :param float evalue: Expected number of false positives with a score
    #:     greater than this hit.
    #: :param float pvalue: Probability in a *pairwise* comparison of finding a
    #:     hit with at least this score.
    #: :param float score: Raw score calculated by HMM-HMM alignment.
    #: :param float ss: Secondary structure score.
    #: :param int cols: Number of aligned match states in the alignment.
    #: :param range(int) query_range: Range of aligned match states in the
    #:     query.
    #: :param range(int) template_range: Range of aligned match states in the
    #:     template.
    #: :param int template_matches: Number of match states in the template HMM.
    #: :param int identities: percentage of aligned residues that are identical.
    #: :param float similarity: Arithmetic mean of substitution scores between
    #:     aligned residues.
    #: :param float sum_probs: Sum of posterior probabilities of all aligned
    #:     pairs.
    #: :param float template_neff: Effective number of sequences in the template
    #:     HMM.
    Hit = namedtuple(
        "Hit", (
            "rank name "
            "prob evalue pvalue score ss cols query_range template_range "
            "template_matches identities similarity sum_probs template_neff"
        ))

    #: Summary of an hhsuite report, stored at the top of the report file.
    #: :param str query: Name of the query sequence/MSA.
    #: :param int match_cols: Number of match columns in the query HMM.
    #: :param float neff: Effective number of sequences.
    #: :param int num_searched: Number of HMMs searched.
    #: :param str date: Date the report file was generated.
    #: :param str command: Command line used to generate the report.
    #: :param tuple(int, int) num_seqs: Number of filtered and input sequences.
    Summary = namedtuple(
        "Summary",
        "query match_cols neff num_searched date command num_seqs")

    # List of regexes, field names and functions to convert to the correct data
    # type. These match against lines in the alignment section at the bottom of
    # the report file.
    _PAIRWISE_FIELDS = [
        (r"Identities=([^%]+)%", "identities", float),
        (r"Similarity=(\S+)", "similarity", float),
        (r"Sum_probs=(\S+)", "sum_probs", float),
        (r"Template_Neff=(\S+)", "template_neff", float),
    ]

    # The scores, in order, appearing in the hit section. These appear *after*
    # the hit rank and (truncated) sequence name. This is a list of 2-tuples
    # containing the name of the field and the function used to parse the string
    # into an object of the correct type.
    _HIT_FIELDS = [
        ("prob", float),
        ("evalue", float),
        ("pvalue", float),
        ("score", float),
        ("ss", float),
        ("cols", int),
        ("query_range", _parse_range),
        ("template_range", _parse_range),
        ("template_matches", int)
    ]

    # Header fields, mapping to a 2-tuple containing the name used in the
    # summary and a conversion function.
    _HEADER_FIELDS = {
        "Query": ('query', str),
        "Match_columns": ('match_cols', int),
        "Neff": ('neff', float),
        "Searched_HMMs": ('num_searched', int),
        "Date": ('date', str),
        "Command": ('command', str),
        "No_of_seqs": ('num_seqs', _parse_num_seqs)
    }

    #Tokens describing the current state of the parser. The state reflects the
    #section of the report that is currently being parsed.
    class State(Enum):
        HEADER   = 1
        HITS     = 2
        PAIRWISE = 3
        DONE     = 4


    def __init__(self, file):
        """Parse a report file."""
        self._state = Report.State.HEADER
        self._current_index = 0
        self._summary = {}
        self._hits = []
        self._parse_file(file)
        self._hit_objs = None

    @property
    def hits(self):
        """Return a list of :py:class:`.Hit` objects."""
        if self._hit_objs is None:
            self._hit_objs = [self.Hit(**hit) for hit in self._hits]
        return self._hit_objs

    @property
    def summary(self):
        """Get summary information about the report."""
        return self.Summary(**self._summary)

    def _parse_file(self, file):
        #Parse the file line by line, maintaining the currents state in the
        #_state instance variable. Start by parsing the header.
        with open(file, "r") as in_fh:
            for line in in_fh:
                line = line.rstrip()
                if self._state == Report.State.HEADER:
                    self._parse_header_line(line)
                elif self._state == Report.State.HITS:
                    self._parse_hit_line(line)
                elif self._state == Report.State.PAIRWISE:
                    self._parse_pairwise_line(line)

    def _parse_header_line(self, line):
        #Scan for state-terminating lines
        if line.startswith(" No Hit"):
            self._state = Report.State.HITS
            return

        #Skip empty lines
        if not line.strip():
            return

        #Parse the initial few lines of the header. This is of the form
        #  ^Field       Value\s*$
        field, value = line.split(maxsplit=1)
        if field in self._HEADER_FIELDS:
            name, converter = self._HEADER_FIELDS[field]
            self._summary[name] = converter(value)
        else:
            log().warn(
                "Unknown field '%s' (value: '%s') in report",
                field, value)

    def _parse_hit_line(self, line):
        # Reproducing a full-width excerpt is worth violating my arbitrary rules
        # on code width:
        # pylint: disable=max-line-length

        #The hitlist looks like this (with column numbers added):
        #
        #  0        1         2         3         4         5         6         7         8         9
        # 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        # No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
        #   1 UP20|FOWREZABA|22|311 Krueppel 100.0  3E-103  2E-108  683.9   0.0  277    1-283    35-311 (311)
        #   2 UP20|MUPDABEBA|4|283 Kruppel-l 100.0  9E-102  6E-107  662.8   0.0  283    1-283     1-283 (283)
        #
        #   1 c2yheD_                        100.0 1.8E-36 4.4E-41  239.3  16.5  253   17-269   100-405 (634)
        #   2 c2cfuA_                        100.0 1.6E-35 3.9E-40  233.6  16.2  253   17-269   100-400 (627)
        # 139 c3izcZ_                          8.9      86  0.0022   11.9  -0.3    6  198-203    13-18  (73)

        #There are some tricks to parsing this output. The thing to note is that
        #sequence names *can* contain spaces, and that column parsing is *not*
        #reliable because long sequences or massive numbers of sequences will throw
        #off the column spacing. We should note that fields are always separated a
        #space, except for the final column which may butt right up against the
        #previous column but is always contained within parentheses.
        #
        #To parse this then, we first consume the rank, which may be of arbitrary
        #length. We then consume the 30 char name. We can then split by a
        #combination of spaces and parentheses and squeeze all empty fields.

        #A blank line terminates this section and begins the pairwise section
        if not line.strip():
            self._state = Report.State.PAIRWISE
            return

        info = {}
        #Get and remove the rank
        info["rank"], line = line.split(maxsplit=1)
        info["rank"] = int(info["rank"])

        #Discard the (possibly truncated) name of the sequence
        line = line[31:]

        # Extract and convert the remaining fields. The regex is used to split
        # to easily discard parentheses around the template_matches field.
        fields = [f for f in re.split("[ ()]+", line) if f]
        for (field, converter), value in zip(self._HIT_FIELDS, fields):
            info[field] = converter(value)

        #Finally, add this hit to the list
        self._hits.append(info)

    def _parse_pairwise_line(self, line):
        #Skip blank lines or "Done!"
        if not line.strip() or line.startswith("Done!"):
            return

        current = self._current_index
        if line.startswith("No "):
            #Lines looking like /^No \d+\s*$/ indicate the start of a new hit
            _, rank = line.split()
            self._current_index = int(rank) - 1
        elif line.startswith(">"):
            #Parse the name of the hit. This can be longer than the hit name in the
            #summary, because it is not truncated:
            #>A long description
            self._hits[current]["name"] = line[1:]
        elif line.startswith("Probab"):
            #Parse the details line:
            #Probab=100.00  E-value=9.5e-116  Score=793.48  Aligned_cols=222  Identities=66%  Similarity=1.106  Sum_probs=219.2  Template_Neff=3.452
            for regex, field, converter in self._PAIRWISE_FIELDS:
                match = re.search(regex, line)
                if match:
                    self._hits[current][field] = converter(match.group(1))


class Tabular:
    """Parse tabulated files from hhsearch (produced using the ``-atab``
    option.

    :param str file: Path of the file to parse.
    :ivar hits: List of hits contained within the file to be parsed.
    :vartype hits: :py:class:`.Hit`
    """

    #: Represents a single hit.
    #: :param str name: Name of this hit. Parsed from the tabular report.
    #: :param list[ResiduePair] alignment: List of residue pairs.
    Hit = namedtuple("Hit", "name alignment")

    #: Alignment between two residues.
    #:
    #: Some of these scores may be ``None``.
    #: :param int i: Index of the query residue.
    #: :param int j: Index of the template residue.
    #: :param float SS: Secondary structure score.
    #: :param float score: Score of this residue pair.
    #: :param float probab: Posterior probability of residue pair.
    #: :param str dssp: DSSP-assigned state of template residue.
    ResiduePair = namedtuple("ResiduePair", "i j score SS probab dssp")

    # Type conversion functions of the various possible properties.
    _TYPES = {
        "i": int, "j": int,
        "score": float, "probab": float, "dssp": str, "SS": float
    }

    def __init__(self, file):
        """Parse a tabulated file produced by an hhsuite tool."""
        self.hits = []
        with open(file, "r") as fh:
            self._parse_file(fh)

    def _records(self, fh):
        """Generator yielding records from a tabulated file.

        Records are separated by lines beginning with ``>``. Each record that we
        yield is a list of lines.
        """
        record = []
        for line in fh:
            if line.startswith(">") and record:
                yield record
                record = []
            record.append(line.rstrip("\n"))
        yield record

    def _parse_file(self, fh):
        #Records start with a line beginning with ">" containing the full name
        #of the hit.
        #
        #The following line (or lines -- we will parse as may as #necessary)
        #will describe which elements of the alignment are missing.  #For
        #example, the line may be "missing dssp", indicating that there isno
        #"dssp" column. We will ignore these, because we will parse the next
        #line, which explicitly tells us which columns are given.
        #
        #The header line gives the list of fields that we will parse. For
        #example:
        #    i     j  score     SS  probab  dssp
        #
        #The remaining lines give the values for each field.

        for lines in self._records(fh):
            missing = []
            name = None
            header = None
            pairs = []
            for line in lines:
                if line.startswith("missing"):
                    missing.append(line.replace("missing ", ""))
                elif name is None:
                    # Strip leading ">"
                    name = line[1:]
                elif header is None:
                    header = line.split()
                else:
                    # Use header to generate a dict of field:value pairs
                    pair_dict = {}
                    for field, value in zip(header, line.split()):
                        pair_dict[field] = self._TYPES[field](value)
                    # Add missing values with a value of None
                    for field in missing:
                        pair_dict[field] = None
                    # Convert to a Hit ResiduePair object
                    pairs.append(self.ResiduePair(**pair_dict))
            self.hits.append(self.Hit(name, pairs))

class Fasta:
    """
    Parse FASTA files from hhsearch (produced using the ``-Ofas`` option.
    Alignments will be stored in the `hits` attribute. Each element of `hits`
    will be a dictionary, containing the keys ``query`` and ``template``. Those
    will contain dictionaries with elements named for each sequence in the FASTA
    file and values corresponding to the sequences

    >>> from phyre_engine.tools.hhsuite.parser import Fasta
    >>> parser = Fasta("/path/to/fasta/file", "QueryName")
    >>> parser.hits
    [
        {
            "query": {
                "ss_pred": "CCHH...",
                "ss_conf": "8899...",
                "Consensus": "XXAI...",
                "sequence": "AGAI..."
            },
            "template": {
                "ss_pred": "CCHH...",
                "ss_conf": "7789...",
                "ss_dssp": "CCHT...",
                "Consensus": "XXAI...",
                "sequence": "GGAI..."
            }
        }
    ]

    :param str file: Path of the file to parse.
    :param str query_name: Name of the query sequence.
    :ivar alns: List of alignments, given as a tuple of query/template pairs.
    :vartype alns: list[dict(str, dict(str, str))]
    """
    def __init__(self, file, query_name):
        """Parse a tabulated file produced by an hhsuite tool."""
        self.hits = []
        self.query_name = query_name
        with Stream(file, "r") as fh:
            self._parse_file(fh)

    def _records(self, fh):
        """Generator yielding records from a tabulated file.

        Records are separated by blank lines.
        """
        record = []
        for line in fh:
            if not line.strip() and record:
                yield "\n".join(record)
                record = []
            if not line.startswith("No "):
                record.append(line.rstrip("\n"))
        yield "\n".join(record)

    def _parse_file(self, fh):
        for i, alignments in enumerate(self._records(fh)):
            aln_dict = {"query": {}, "template": {}}
            try:
                # The alignment is split into sections. Everything before the
                # query sequence goes in the query, everything after in the
                # template.
                section = "query"

                parsed_alns = Bio.AlignIO.read(io.StringIO(alignments), "fasta")
                for i, aln in enumerate(parsed_alns):
                    key_name = aln.id
                    if aln.id == self.query_name or i == len(parsed_alns) - 1:
                        # If the name of the sequence is the query name, then
                        # this is the query sequence. If the alignment is the
                        # last in the list, this is the template sequence.
                        key_name = "sequence"
                    aln_dict[section][key_name] = str(aln.seq)

                    # Everything before the query is query-specific
                    if aln.id == self.query_name:
                        section = "template"
            except (IndexError, ValueError):
                log().warning("Error parsing alignment %d", i, exc_info=True)
            finally:
                self.hits.append(aln_dict)
