import io
import unittest
import tempfile
import textwrap
import Bio.SeqRecord
import Bio.Seq
import phyre_engine.component.input as input_cpts

class TestReadSingleSequence(unittest.TestCase):
    """Test ReadSingleSequence component."""

    def test_valid(self):
        """Create a valid FASTA file and read it back."""
        seq_lines = [
            'AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL',
            'SGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLY',
            'THMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVE',
            'AIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAK',
            'GRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG',
            'HAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDR',
            'LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVW',
            'PAAVRESVPSLL']

        fasta = ">An arbitrary identifer\n{}".format("\n".join(seq_lines))

        valid_fasta_file = io.StringIO(fasta)

        fasta_input = input_cpts.ReadSingleSequence()
        seq = fasta_input.run({"input": valid_fasta_file})
        self.assertEqual(seq["sequence"], "".join(seq_lines))

    def test_file_exists(self):
        """Check that an IOError is raised if the input file does not exist."""
        with self.assertRaises(IOError):
            input_cpts.ReadSingleSequence().run(
                {"input": "bad_file_name"})

    def test_single_seq_required(self):
        """Require a single sequence only."""

        #Create a multi-sequence FASTA file
        fasta = textwrap.dedent("""\
                >ID 1
                AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL
                >ID 2
                SGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLY
                """)


        msa_fasta_file = io.StringIO(fasta)

        exception_type = input_cpts.ReadSingleSequence.TooManySequencesError
        with self.assertRaises(exception_type):
            input_cpts.ReadSingleSequence().run({"input": msa_fasta_file})

class TestReadMultipleSequences(unittest.TestCase):
    """Test ReadMultipleSequences component."""

    _FASTA = textwrap.dedent("""\
    >FOO X Y Z
    AAAGGG
    >BAR X Y Z
    HHHEEE
    """)

    def test_multi_seq(self):
        """Read multiple sequences."""
        reader = input_cpts.ReadMultipleSequences()
        fas_out = io.StringIO(self._FASTA)
        results = reader.run({"input": fas_out})

        self.assertEqual(results["templates"][0]["sequence"], "AAAGGG")
        self.assertEqual(results["templates"][0]["name"], "FOO")
        self.assertEqual(results["templates"][0]["id"], "FOO")
        self.assertEqual(results["templates"][0]["description"], "FOO X Y Z")

        self.assertEqual(results["templates"][1]["sequence"], "HHHEEE")
        self.assertEqual(results["templates"][1]["name"], "BAR")
        self.assertEqual(results["templates"][1]["id"], "BAR")
        self.assertEqual(results["templates"][1]["description"], "BAR X Y Z")
