import unittest
import os
import tempfile

import phyre_engine.component.db.clusters as clusters
import phyre_engine.test

class TestRCSBClusterDownload(unittest.TestCase):
    """Test downloading cluster file from the RCSB"""

    @unittest.skipUnless(
        (phyre_engine.test.config is not None
            and phyre_engine.test.config.get("net_tests", False)),
        "Not doing network tests (set 'net_tests' configuration parameter).")
    def test_retrieve_and_parse(self):
        """Test download and parse of RCSB clusters."""

        this_dir = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as test_dir:
                os.chdir(test_dir)

                rcsb_clus = clusters.RCSBClusterDownload(100, "test_file")
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
            clusters.RCSBClusterDownload(12)

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

            clus_parser = clusters.ClusterParser()
            results = clus_parser.run({"cluster_file": clus_file.name})

            self.assertIn("clusters", results, "Added 'clusters' key")

            seq_clusters = results["clusters"]
            self.assertGreater(
                len(seq_clusters), 0,
                "Got more than one cluster")
            self.assertListEqual(
                    seq_clusters[0],
                    "3J34_5 3J34_6 3J34_7 3J34_A 3J34_B 3J34_C 3J34_D".split(),
                    "Cluster 0 correct")
            self.assertListEqual(
                    seq_clusters[4],
                    ["4LCD_E"],
                    "Cluster 4 correct")

class TestSimpleRepresentativePicker(unittest.TestCase):
    """Test the (trivially) simple representative picker.."""

    def test_parse(self):
        """Pick representatives from clusters."""
        seq_clusters = [
            "3J34_5 3J34_6 3J34_7 3J34_A 3J34_B 3J34_C 3J34_D".split(),
            "1CE6_B 4X6C_B 4X6C_D 4X6E_B".split(),
            "1LSG_A 5JEN_B 5JEN_D".split(),
            "2KXK_A 2LGB_A".split(),
            "4LCD_E".split()
        ]

        srp = clusters.SimpleRepresentativePicker()
        representatives = srp.run({"clusters": seq_clusters})["templates"]
        self.assertListEqual(
                representatives, [
                    {"PDB": "3J34", "chain": "5"},
                    {"PDB": "1CE6", "chain": "B"},
                    {"PDB": "1LSG", "chain": "A"},
                    {"PDB": "2KXK", "chain": "A"},
                    {"PDB": "4LCD", "chain": "E"}
                ], "Correctly picked representatives")
