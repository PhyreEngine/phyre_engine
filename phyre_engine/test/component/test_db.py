import json
import os
import tempfile
import unittest
import phyre_engine.component.db as db
from pathlib import Path

class TestRCSBClusterDownload(unittest.TestCase):
    """Test downloading cluster file from the RCSB"""

    @unittest.skipUnless("NET_TESTS" in os.environ, "Not doing network tests")
    def test_retrieve_and_parse(self):
        """Test download and parse of RCSB clusters."""

        this_dir = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as test_dir:
                os.chdir(test_dir)

                rcsb_clus = db.RCSBClusterDownload(100, "test_file")
                results = rcsb_clus.run({})

                self.assertIn(
                        "cluster_file",
                        results,"Added 'cluster_file' key")
                self.assertTrue(
                        os.path.isfile("test_file"),
                        "Created 'test_file'")
                self.assertGreater(
                        os.path.getsize("test_file"),
                        0, "Download file contains data")
        finally:
            os.chdir(this_dir)

    def test_requires_valid_thresh(self):
        with self.assertRaises(ValueError, msg="Requires valid threshold"):
            db.RCSBClusterDownload(12)

class TestClusterParser(unittest.TestCase):
    """Test parsing a cluster file."""

    def test_parse(self):
        """Try and parse a dummy file."""

        file_contents = (
            "3J34_5 3J34_6 3J34_7 3J34_A 3J34_B 3J34_C 3J34_D \n"
            "1CE6_B 4X6C_B 4X6C_D 4X6E_B \n"
            "1LSG_A 5JEN_B 5JEN_D\n"
            "2KXK_A 2LGB_A \n"
            "4LCD_E \n"
        )

        with tempfile.NamedTemporaryFile(mode="w") as clus_file:
            clus_file.write(file_contents)
            clus_file.flush()

            clus_parser = db.ClusterParser()
            results = clus_parser.run({"cluster_file": clus_file.name})

            self.assertIn("clusters", results, "Added 'clusters' key")

            clusters = results["clusters"]
            self.assertGreater(len(clusters), 0, "Got more than one cluster")
            self.assertListEqual(
                    clusters[0],
                    "3J34_5 3J34_6 3J34_7 3J34_A 3J34_B 3J34_C 3J34_D".split(),
                    "Cluster 0 correct")
            self.assertListEqual(
                    clusters[4],
                    ["4LCD_E"],
                    "Cluster 4 correct")

class TestSimpleRepresentativePicker(unittest.TestCase):
    """Test the (trivially) simple representative picker.."""

    def test_parse(self):
        """Pick representatives from clusters."""
        clusters = [
            "3J34_5 3J34_6 3J34_7 3J34_A 3J34_B 3J34_C 3J34_D".split(),
            "1CE6_B 4X6C_B 4X6C_D 4X6E_B".split(),
            "1LSG_A 5JEN_B 5JEN_D".split(),
            "2KXK_A 2LGB_A".split(),
            "4LCD_E".split()
        ]

        srp = db.SimpleRepresentativePicker()
        representatives = srp.run({"clusters": clusters})["templates"]
        self.assertListEqual(
                representatives, [
                    {"PDB": "3J34", "chain": "5"},
                    {"PDB": "1CE6", "chain": "B"},
                    {"PDB": "1LSG", "chain": "A"},
                    {"PDB": "2KXK", "chain": "A"},
                    {"PDB": "4LCD", "chain": "E"}
                ], "Correctly picked representatives")

class TestChainPDBBuilder(unittest.TestCase):
    """Test ChainPDBBuilder class"""

    @unittest.skipUnless("MMCIF" in os.environ, "MMCIF env var not set")
    def test_build(self):
        """Try and extract a chain"""
        with tempfile.TemporaryDirectory() as out_dir:
            builder = db.ChainPDBBuilder(os.environ["MMCIF"], out_dir)
            results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

            self.assertTrue(Path(out_dir, "2a").exists(), "Created 2a subdir")
            self.assertTrue(
                Path(out_dir, "2a/12as_A.pdb").exists(),
                "Created 12as_A"
            )
            self.assertEqual(
                results["templates"][0]["sequence"].seq, (
                "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
                "AVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDE"
                "DRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSE"
                "EFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGG"
                "KLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMG"
                "IRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTML"
                "LLQLPHIGQVQAGVWPAAVRESVPSLL"),
                "Got correct sequence")
