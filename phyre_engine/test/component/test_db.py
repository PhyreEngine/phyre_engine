import os
import tempfile
import unittest
import phyre_engine.component.db as db

class TestRCSBClusterDownload(unittest.TestCase):
    """Test downloading cluster file from the RCSB"""

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
