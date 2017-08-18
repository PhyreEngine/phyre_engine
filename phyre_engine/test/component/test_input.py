import unittest
import tempfile
import textwrap
import Bio.SeqRecord
import Bio.Seq
from phyre_engine.component.input import (FastaInput, MultipleFastaInput,
                                          ConvertSeqRecord)

class TestFastaInput(unittest.TestCase):
    """Test FastaInput component."""

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
        seq = fasta_input.run({"input": valid_fasta_file.name})['seq_record']
        self.assertEqual(seq.seq, "".join(seq_lines))

    def test_file_exists(self):
        """Check that an IOError is raised if the input file does not exist."""
        with self.assertRaises(IOError):
            FastaInput().run({"input": "bad_file_name"})

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
            FastaInput().run({"input":msa_fasta_file.name})

class TestMultipleFastaInput(unittest.TestCase):
    """Test MultipleFastaInput component."""

    _FASTA = textwrap.dedent("""\
    >FOO
    AAAGGG
    >BAR
    HHHEEE
    """)

    def test_multi_seq(self):
        """Read multiple sequences."""
        reader = MultipleFastaInput()
        with tempfile.NamedTemporaryFile("w") as fas_out:
            print(self._FASTA, file=fas_out)
            fas_out.flush()
            results = reader.run({"input": fas_out.name})

            self.assertEqual(
                results["templates"][0]["seq_record"].seq,
                "AAAGGG")
            self.assertEqual(
                results["templates"][0]["seq_record"].name,
                "FOO")

            self.assertEqual(
                results["templates"][1]["seq_record"].seq,
                "HHHEEE")
            self.assertEqual(
                results["templates"][1]["seq_record"].name,
                "BAR")

class TestConvertSeqRecord(unittest.TestCase):
    """Test ConvertSeqRecord component."""

    _SEQ_RECORD = Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("AAAGGG"),
        "ID", "NAME", "DESCRIPTION")

    def test_conversion(self):
        """Test ConvertSeqRecord on a single record."""
        converter = ConvertSeqRecord()
        results = converter.run({"seq_record": self._SEQ_RECORD})
        del results["seq_record"]
        self.assertDictEqual(
            results, {
                "sequence": "AAAGGG",
                "id": "ID",
                "name": "NAME",
                "description": "DESCRIPTION"
            })
