import unittest
import textwrap
from data.Sequence import Sequence;

class TestSequence(unittest.TestCase):
    """Unit test Sequence class."""

    def test_valid_seq(self):
        """Instantiate a valid sequence."""
        seq_str = "AYIAKQRQISFVKSHFSRQLEERLGLIEV"
        seq = Sequence(seq_str)
        self.assertEqual(str(seq), seq_str)

    def test_invalid_seq(self):
        """Fail to build an invalid sequence."""

        #Note the "X" residue here:
        seq_str = "XAYIAKQRQISFVKSHFSRQLEERLGLIEV"
        with self.assertRaises(Sequence.UnknownAminoAcidError):
            Sequence(seq_str)

    def test_extra_residue_types(self):
        """Change list of allowed extra residue types."""

        seq_str = "XBBXZ"
        try:
            Sequence(seq_str, allowed_extra=["X", "B", "Z"])
        except:
            self.fail("Exception raised for valid sequence")
