"""Test components in the util package."""
import os
import pathlib
import tempfile
import unittest
import phyre_engine.component.util as util
import shutil
from pathlib import Path

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

class TestMakeDir(unittest.TestCase):
    """Test the MakeDir component."""

    def setUp(self):
        """Create tmpdir for testing."""
        self.tmpdir = tempfile.mkdtemp("-test", "mkdir-")

    def tearDown(self):
        """Remoev tmpdir."""
        shutil.rmtree(self.tmpdir)

    def test_mkdir_parents(self):
        """Create a nested subdirectory."""
        mkdir = util.MakeDir(exist_ok=True, parents=True)
        target = Path(self.tmpdir, "a/b/c")
        mkdir.run({"working_dir": target})
        self.assertTrue(target.exists(), "Created {}".format(target))

    def test_mkdir_no_exist(self):
        """Cannot create a directory that already exists when exist_ok=False."""
        mkdir = util.MakeDir(exist_ok=False, parents=True)

        target = Path(self.tmpdir, "a/b/c")
        mkdir.run({"working_dir": target})
        self.assertTrue(target.exists(), "Created {}".format(target))

        with self.assertRaises(FileExistsError):
            mkdir.run({"working_dir": target})

    def test_mkdir_no_parents(self):
        """Creating a nested subdirectory should fail with parents=False."""
        mkdir = util.MakeDir(exist_ok=True, parents=False)

        # Nested creation fails
        with self.assertRaises(FileNotFoundError):
            target = Path(self.tmpdir, "a/b/c")
            mkdir.run({"working_dir": target})

        # Non-nested works
        target = Path(self.tmpdir, "a")
        mkdir.run({"working_dir": target})
        self.assertTrue(target.exists(), "Created {}".format(target))
