"""Test components in the util package."""
import os
import pathlib
import tempfile
import unittest
import phyre_engine.component.util as util

class TestChangeDir(unittest.TestCase):
    """Test the ChangeDir component."""

    def test_chdir(self):
        """Change directory to a temporary directory."""
        orig_dir = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as tempdir:
                chdir_cmpt = util.ChangeDir()
                chdir_cmpt.run({"working_dir": tempdir})
                self.assertEqual(
                    pathlib.Path(os.getcwd()).resolve(),
                    pathlib.Path(tempdir).resolve())
        finally:
            os.chdir(orig_dir)
