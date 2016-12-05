"""Parsers for hhsuite files."""
from enum import Enum
import warnings
import re

class Report:
    """Parse an hhsearch report file into a list of hits.

    This class will only parse data on a per-hit basis; it will not parse the
    actual alignments from the report file. To get the alignments, you should
    generate a tabular file using the ``-atab`` option of hhblits/hhsearch and
    parse the alignment data using the :class:`Tabular` parser.

    :param str file: Path of the report file to be parsed.
    :ivar hits: List of :class:`Hit` objects.
    :ivar summary: Dictionary containing the summary lines of the report.
    """

    #Tokens describing the current state of the parser. The state reflects the
    #section of the report that is currently being parsed.
    class State(Enum):
        HEADER   = 1
        HITS     = 2
        PAIRWISE = 3
        DONE     = 4


    def __init__(self, file):
        """Parse a report file."""
        self._parse_file(file)


    def _parse_file(self, file):
        #Parse the file line by line, maintaining the currents state in the
        #_state instance variable. Start by parsing the header:
        self.summary = {}
        self._state = Report.State.HEADER
        self.hits = []

        #Index of the hit being parsed by _parse_pairwise_line
        self._current_index = 0

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
        #so we can just grab the results using (several) regexes.

        #The "No_of_seqs" field is slightly different, because it is of the
        #follwing form:
        #  ^No_of_seqs 137 out of 2499$
        #We parse this into the following hashref:
        #  {filtered => 137, in => 2499}
        fields = {
            "Query":          'query',
            "Match_columns":  'match_cols',
            "Neff":           'neff',
            "Searched_HMMs":  'num_searched',
            "Date":           'date',
            "Command":        'command'
        }
        ints = {"num_searched", "match_cols"}
        floats = {"neff"}
        field, value = line.split(maxsplit=1)

        if field == "No_of_seqs":
            num_filt, num_in = value.split(" out of ")
            self.summary["num_seqs"] = {}
            self.summary["num_seqs"]["filtered"] = int(num_filt)
            self.summary["num_seqs"]["in"]       = int(num_in)
        elif field in fields:
            short_name = fields[field]
            if short_name in ints:
                value = int(value)
            elif short_name in floats:
                value = float(value)
            self.summary[short_name] = value
        else:
            err_msg = "Unknown field {} in hhhsuite report."
            warnings.warn(err_msg.format(field))


    def _parse_hit_line(self, line):
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

        hit = Hit()

        #Get and remove the rank
        hit.info["rank"], line = line.split(maxsplit=1)

        #Discard the (possibly truncated) name of the sequence
        line = line[31:]

        #Split on spaces and parentheses, discard empty fields and fill the
        #hit dict
        (
            hit.info["prob"],
            hit.info["evalue"],
            hit.info["pvalue"],
            hit.info["score"],
            hit.info["ss"],
            hit.info["cols"],
            hit.info["query_range"],
            hit.info["template_range"],
            hit.info["template_matches"],
        ) = [x for x in re.split("[ ()]+", line) if x]

        #Turn query_range and template_range into ranges
        hit.info["template_range"] = range(
                *[int(x) for x in hit.info["template_range"].split("-")])
        hit.info["query_range"] = range(
                *[int(x) for x in hit.info["query_range"].split("-")])

        #Convert the remaining values into numbers.
        for field in ["rank", "cols", "template_matches"]:
            hit.info[field] = int(hit.info[field])
        for field in ["prob", "evalue", "pvalue", "score", "ss"]:
            hit.info[field] = float(hit.info[field])

        #Finally, add this hit to the list
        self.hits.append(hit)


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
            self.hits[current].info["name"] = line[1:]
        elif line.startswith("Probab"):
            #Parse the details line:
            #Probab=100.00  E-value=9.5e-116  Score=793.48  Aligned_cols=222  Identities=66%  Similarity=1.106  Sum_probs=219.2  Template_Neff=3.452

            #List of regexes, field names and functions to convert to the
            #correct data type.
            matches = [
                ("Identities=([^%]+)%", "identities",    float),
                ("Similarity=(\S+)",    "similarity",    float),
                ("Sum_probs=(\S+)",     "sum_probs",     float),
                ("Template_Neff=(\S+)", "template_neff", float),
            ]

            for regex, field, converter in matches:
                match = re.search(regex, line)
                if match:
                    self.hits[current].info[field] = converter(match.group(1))


class Hit:
    """A class representing a single hit.

    A hit may have several per-hit attributes describing, for example, the score
    between a query and the hit or simply giving the name of the hit.
    Additionally, a hit may contain the alignment between a query and this hit.

    :param info: Optional keys to be added to the `info` attribute.
    :ivar info: A dictionary giving per-hit properties.
    :ivar aln: A list of alignment pairs.
    """

    def __init__(self, **info):
        """Initialise a new Hit method with empty `info` and `aln`
        attributes.
        """

        self.info = info
        self.aln = []


class Tabular:
    """Parse tabulated files from hhsearch (produced using the ``-atab``
    option.

    :param str file: Path of the file to parse.
    :ivar hits: List of hits contained within the file to be parsed.
    """


    def __init__(self, file):
        """Parse a tabulated file produced by an hhsuite tool."""

        self.hits = []
        self._current_hit = None
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

    def _parse_num(self, string):
        #Parse a number from a string. If possible, this will return an int.
        #Otherwise, a float. If the string is not numeric, it will return the
        #string.
        try:
            return int(string)
        except ValueError:
            pass

        try:
            return float(string)
        except ValueError:
            pass

        return string

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
            lines = [l for l in lines if not l.startswith("missing")]

            name   = lines[0][1:]
            header = lines[1].split()
            hit = Hit(name=name)
            for l in lines[2:]:
                values = [self._parse_num(x) for x in l.split()]
                hit.aln.append(dict(zip(header, values)))
            self.hits.append(hit)

