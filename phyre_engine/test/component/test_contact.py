"""Tests for :py:mod:`phyre_engine.component.contact` module."""
import io
import os
import tempfile
import textwrap
import unittest
import unittest.mock

import phyre_engine.test
import phyre_engine.component.contact as contact

class TestMetaPsicov(unittest.TestCase):
    """Test MetaPsicov component."""

    TEST_SEQ = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHS"
    "LAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGI"
    "KATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAP"
    "DYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMP"
    "QTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSL")

    def test_parse(self):
        """Test parsing of stage3 file."""
        test_file = io.StringIO(textwrap.dedent("""\
        1 2 0 8 0.9951
        3 10 0 8 0.02
        4 1 0 8 0.35
        """))
        # Convert to tuple just so we can easily compare equality
        self.assertEqual(
            [tuple(x) for x in contact.MetaPsicov.parse_results(test_file)],
            [(1, 2, 0.9951), (3, 10, 0.02), (4, 1, 0.35)])

    @phyre_engine.test.requireFields(
        ["bin_dir", "uniref90", "uniref100", "hhblitsdb", "pdb70"],
        ["tools", "metapsicov"])
    def test_metapsicov(self):
        """Run MetaPsicov."""
        orig_dir = os.getcwd()
        with tempfile.TemporaryDirectory("-metapsicov",
                                         "phyreengine-") as tmpdir:
            try:
                os.chdir(tmpdir)
                config = phyre_engine.test.config["tools"]["metapsicov"]
                metapsicov = contact.MetaPsicov(
                    config["uniref90"], config["uniref100"],
                    config["hhblitsdb"], config["pdb70"],
                    bin_dir=config["bin_dir"],
                    cpu=15)
                results = metapsicov.run({"sequence": self.TEST_SEQ})
                self.assertGreater(len(results["contacts"]), 1)
            finally:
                os.chdir(orig_dir)
