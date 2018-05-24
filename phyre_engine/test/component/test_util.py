"""Test components in the util package."""
import os
import pathlib
import tempfile
import unittest
import unittest.mock
import phyre_engine.component.util as util
import shutil
from pathlib import Path


class TestSaveRootDir(unittest.TestCase):
    """Test SaveRootDir component."""

    def test_save_root(self):
        """Test saving root directory."""
        save_cpt = util.SaveRootDir()
        result = save_cpt.run({})
        self.assertEqual(result["root_directory"], os.getcwd())


class TestRestoreRootDir(unittest.TestCase):
    """Test RestoreRootDir component."""

    @unittest.mock.patch("os.chdir")
    def test_restore_root(self, chdir_mock):
        """Test restoring root directory."""
        restore_cpt = util.RestoreRootDir()
        restore_cpt.run({"root_directory": "/foo/bar"})
        chdir_mock.assert_called_once_with("/foo/bar")


class TestChangeDir(unittest.TestCase):
    """Test the ChangeDir component."""

    def test_chdir(self):
        """Change directory to a temporary directory."""
        orig_dir = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as tempdir:
                chdir_cmpt = util.ChangeDir()
                result = chdir_cmpt.run({"working_dir": tempdir})
                self.assertEqual(
                    pathlib.Path(os.getcwd()).resolve(),
                    pathlib.Path(tempdir).resolve())
                self.assertEqual(
                    result["directory_stack"],
                    [orig_dir])
        finally:
            os.chdir(orig_dir)


class TestPopDir(unittest.TestCase):
    """Test PopDir component."""

    def test_popdir(self):
        """Change directory to a previous directory."""
        with unittest.mock.patch("os.chdir") as chdir_mock:
            pop_cmpt = util.PopDir()
            state = {"directory_stack": ["a", "b", "c"]}
            result = pop_cmpt.run(state)
            result = pop_cmpt.run(state)
            result = pop_cmpt.run(state)
            result = pop_cmpt.run(state)

            # Moved back through the directory stack
            self.assertEqual(chdir_mock.call_count, 3)
            chdir_mock.assert_has_calls([
                unittest.mock.call("c"),
                unittest.mock.call("b"),
                unittest.mock.call("a"),
            ])
            # Cleared stack when empty
            self.assertNotIn("directory_stack", "state")


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


@unittest.mock.patch("pathlib.Path.symlink_to", spec=Path)
class TestSymlink(unittest.TestCase):
    """Test the Symlink component."""

    def test_symlink(self, symlink_mock):
        """Symbolic link is created with correct name pointing to target."""
        orig_name = "awkward_name.pdb"
        final_name = "01-1ABC_X.final.pdb"

        cmpt = util.Symlink(
            name="{rank:02d}-{PDB}_{chain}.final.pdb",
            target="model")
        result = cmpt.run({
            "PDB": "1ABC",
            "chain": "X",
            "rank": 1,
            "model": orig_name})

        symlink_mock.assert_called_once_with(orig_name)
        self.assertEqual(result["model"], final_name)
