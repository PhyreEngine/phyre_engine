"""Test components in the :py:mod:`phyre_engine.component.hhsuite` module."""
import copy
import io
import textwrap
import unittest
import unittest.mock
import phyre_engine.component.hhsuite as hhsuite
import phyre_engine.tools.template as template

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
    AAAAA
    >Template
    -A-A-

    No 2
    >Query
    AA
    >Template
    GG
    """)

    PIPELINE = {
        "sequence": "XXAAAAAAXX",
        "name": "Query",
        "templates": [
            {"query_range": range(3, 7), "name": "foo"},
            {"query_range": range(3, 4), "name": "bar"},
        ]
    }

    def setUp(self):
        """Create a copy of the pipeline."""
        self.pipeline = copy.deepcopy(self.PIPELINE)
        self.pipeline["pairwise_fasta"] = io.StringIO(self.FASTA)

    def test_parser(self):
        """Parse a FASTA alignment file into the fasta_alignment field."""
        parser = hhsuite.FastaParser()
        results = parser.run(self.pipeline)
        self.assertEqual(
            results["templates"][0]["fasta_alignment"],
            ">Query\nXXAAAAAAXX\n>foo\n---A-A----\n")

        self.assertEqual(
            results["templates"][1]["fasta_alignment"],
            ">Query\nXXAAAAAAXX\n>bar\n--GG------\n")
