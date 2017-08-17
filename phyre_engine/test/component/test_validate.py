import unittest
from phyre_engine.component.validate import SeqValidator

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

        validator = SeqValidator()
        validator.run({"sequence": valid_seq_chrs})  #No exception thrown

        with self.assertRaises(SeqValidator.InvalidSeqError) as context:
            validator.run({"sequence": valid_seq_chrs + "X"})
        self.assertSetEqual(context.exception.invalid, {"X"})
