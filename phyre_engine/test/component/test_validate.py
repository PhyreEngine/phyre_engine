import unittest
from phyre_engine.component.validate import SeqValidator, SeqLenFilter

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

class TestSeqLenFilter(unittest.TestCase):
    """Test sequence length filter (SeqLenFilter)."""

    def test_too_short(self):
        """Exception raised when sequence is too short."""
        filt = SeqLenFilter(5, 10)
        with self.assertRaises(SeqLenFilter.SeqLenError):
            filt.run({"sequence": "A" * 4})

    def test_too_long(self):
        """Exception raised when sequence is too long."""
        filt = SeqLenFilter(5, 10)
        with self.assertRaises(SeqLenFilter.SeqLenError):
            filt.run({"sequence": "A" * 11})

    def test_just_right(self):
        """No exception raised when sequence is correct length."""
        filt = SeqLenFilter(5, 10)
        filt.run({"sequence": "A" * 5})
        filt.run({"sequence": "A" * 7})
        filt.run({"sequence": "A" * 10})

    def test_open_max(self):
        """Test sequence length filter with no maximum."""
        filt = SeqLenFilter(5, None)
        filt.run({"sequence": "A" * 1000})

    def test_open_min(self):
        """Test sequence length filter with no minimum."""
        filt = SeqLenFilter(None, 10)
        filt.run({"sequence": ""})
