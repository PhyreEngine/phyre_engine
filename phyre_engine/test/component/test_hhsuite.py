"""Test components in the :py:mod:`phyre_engine.component.hhsuite` module."""
import copy
import io
import textwrap
import unittest
import unittest.mock
import phyre_engine.component.hhsuite as hhsuite
import phyre_engine.tools.template as template
import collections

class TestAlignmentToFasta(unittest.TestCase):
    """Test the AlignmentToFasta component."""

    def test_template_seq(self):
        """Build template_obj part of the alignment from i,j pairs."""
        template_obj = unittest.mock.MagicMock(spec=template.Template)
        template_obj.canonical_seq = "ABCDEFGH"
        query_seq = "JKLMNOP"
        aln = [(1, 1), (3, 5)]

        template_seq = hhsuite.AlignmentToFasta.build_template_aln(
            query_seq, aln, template_obj)
        self.assertEqual(template_seq, "A-E----")

        # Test formatting of the fasta using different fields.
        pipe_state = {"sequence": query_seq, "name": "queryname"}
        hit = {"name": "templatename"}

        aln2fasta = hhsuite.AlignmentToFasta()
        self.assertEqual(
            aln2fasta.fasta(pipe_state, hit, template_seq),
            ">Query\nJKLMNOP\n>Template\nA-E----\n")

        aln2fasta = hhsuite.AlignmentToFasta(q_name="name", t_name="name")
        self.assertEqual(
            aln2fasta.fasta(pipe_state, hit, template_seq),
            ">queryname\nJKLMNOP\n>templatename\nA-E----\n")

class TestFastaParser(unittest.TestCase):
    """Test the FastaParser component."""

    FASTA = textwrap.dedent("""\
    No 1
    >Query
    AB-CDE
    >Template
    -A-A--

    No 2
    >Query
    AB
    >foo
    -A
    >Template
    GG
    """)

    PIPELINE = {
        "sequence": "XXABCDEFXX",
        "name": "Query",
        "templates": [
            {"query_range": range(3, 7)},
            {"query_range": range(3, 4)},
        ]
    }

    Pair = collections.namedtuple("Pair", "i j probab")
    ALIGNMENTS = [
        [Pair(4, 1, 0.5), Pair(5, 2, 0.4)],
        [Pair(5, 4, 0.3)],
    ]

    def setUp(self):
        """Create a copy of the pipeline."""
        self.pipeline = copy.deepcopy(self.PIPELINE)
        self.pipeline["pairwise_fasta"] = io.StringIO(self.FASTA)

    def test_parser(self):
        """Parse a FASTA alignment file into the fasta_alignment field."""
        parser = hhsuite.FastaParser()
        results = parser.run(self.pipeline)
        self.assertEqual(
            results["templates"][0]["sequence_alignments"], {
                "sequence": "---A-A-----",
                "query": "XXAB-CDEFXX"
        })

        self.assertEqual(
            results["templates"][1]["sequence_alignments"], {
                "foo": "---A------",
                "sequence": "--GG------",
                "query": "XXABCDEFXX"
        })

    def test_ignore(self):
        """Ignored sequences should be removed."""
        parser = hhsuite.FastaParser(ignore={"foo"})
        results = parser.run(self.pipeline)
        self.assertNotIn("foo", results["templates"][1]["sequence_alignments"])

    def test_confidences(self):
        """Test generation of confidence string from alignment pairs."""

        # Add alignments to pipeline
        for hit, aln in zip(self.pipeline["templates"], self.ALIGNMENTS):
            hit["alignment"] = aln

        parser = hhsuite.FastaParser()
        results = parser.run(self.pipeline)
        self.assertEqual(
            results["templates"][0]["sequence_alignments"]["confidence"],
            "---5-4-----")
        self.assertEqual(
            results["templates"][1]["sequence_alignments"]["confidence"],
            "----3-----")

class TestA3MSSParser(unittest.TestCase):
    """Test the A3MSSParser."""

    SAMPLE_A3M = textwrap.dedent("""\
    >ss_dssp DSSP Secondary structure
    CCCTHH
    >ss_pred Secondary structure from PSIPRED
    CCCCHH
    >ss_conf SS Confidence
    888889
    >Query
    AAAAGA
    """)

    def test_parser(self):
        """Test A3M parser."""
        parser = hhsuite.A3MSSParser()

        mock_open = unittest.mock.mock_open(read_data=self.SAMPLE_A3M)
        # Bit of a hack to work around a bug: mock_open doesn't provide an
        # iterator on its returned MagicMock.
        mock_open.return_value.__iter__ = lambda self: iter(self.readline, '')
        open_function = ".".join([parser.__module__, "open"])

        with unittest.mock.patch(open_function, mock_open):
            results = parser.run({"a3m": "query.a3m"})
            self.assertEqual(
                results["secondary_structure_sequence"], {
                    "ss_dssp": "CCCTHH",
                    "ss_pred": "CCCCHH",
                    "ss_conf": "888889"
                })
