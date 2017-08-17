import unittest
from phyre_engine.component.validate import SeqValidator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class TestSeqValidator(unittest.TestCase):
    """Test sequence validator (SeqValidator) component."""

    def test_valid(self):
        valid_seq_chrs = (
            'AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL'
            'SGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLY'
            'THMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVE'
            'AIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAK'
            'GRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG'
            'HAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDR'
            'LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVW'
            'PAAVRESVPSLL')

        valid_seq_obj   = SeqRecord(Seq(valid_seq_chrs, IUPAC.protein))
        invalid_seq_obj = SeqRecord(Seq(valid_seq_chrs + "X", IUPAC.protein))

        validator = SeqValidator()
        validator.run({"sequence": valid_seq_obj}) #No exception thrown

        with self.assertRaises(SeqValidator.InvalidSeqError) as context:
            validator.run({"sequence": invalid_seq_obj})
        self.assertSetEqual(context.exception.invalid, {"X"})
