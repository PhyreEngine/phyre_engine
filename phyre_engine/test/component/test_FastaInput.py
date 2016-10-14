import unittest
import tempfile
import textwrap
from phyre_engine.component.FastaInput import FastaInput;

class TestFastaInput(unittest.TestCase):

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

        valid_fasta_file = tempfile.NamedTemporaryFile("w")
        valid_fasta_file.write(fasta)
        valid_fasta_file.flush()


        fasta_input = FastaInput()
        seq = fasta_input.run(input=valid_fasta_file.name)
        self.assertEqual(str(seq), "".join(seq_lines))

    def test_file_exists(self):
        """Check that an IOError is raised if the input file does not exist."""
        with self.assertRaises(IOError):
            FastaInput().run(input="bad_file_name")

    def test_single_seq_required(self):
        """Require a single sequence only."""

        #Create a multi-sequence FASTA file
        fasta = textwrap.dedent("""\
                >ID 1
                AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL
                >ID 2
                SGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLY
                """)


        msa_fasta_file = tempfile.NamedTemporaryFile("w")
        msa_fasta_file.write(fasta)
        msa_fasta_file.flush()

        with self.assertRaises(FastaInput.TooManySequencesError):
            FastaInput().run(input=msa_fasta_file.name)
