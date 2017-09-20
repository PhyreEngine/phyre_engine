"""Test components in the :py:mod:`phyre_engine.component.hhsuite` module."""
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
