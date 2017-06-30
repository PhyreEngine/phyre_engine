import os
import unittest
import tempfile
import subprocess
import shutil
from pathlib import Path
from phyre_engine.component.update.hhsuite import GitUpdate, Compile
import sys


@unittest.skipUnless("LONG_TESTS" in os.environ, "Long tests disabled")
class HHSuiteUpdateTest(unittest.TestCase):
    """Test update and compilation of hh-suite."""

    REPO_URL = "https://github.com/soedinglab/hh-suite.git"

    @classmethod
    def setUpClass(cls):
        try:
            cls._src_dir = tempfile.mkdtemp("hhsuite-compile", "test")
            cls._dest_dir = tempfile.mkdtemp("hhsuite-install", "test")
            subprocess.run(
                ["git", "clone", cls.REPO_URL, cls._src_dir],
                check=True)

            # Get the hash of the HEAD commit so we know that the update brings
            # us back there.
            cls._head_commit = cls.head_commit_hash()

            # Move back one commit so we have to actually do something
            subprocess.run(
                ["git",
                 "--git-dir", str(Path(cls._src_dir, ".git")),
                 "--work-tree", cls._src_dir,
                 "reset", "HEAD~"],
                check=True)

            # Make a change to HHPaths.pm so we can test stashing
            with Path(cls._src_dir, "scripts/HHPaths.pm").open("a") as hhpath:
                print("# Comment from phyre engine tests", file=hhpath)
        except Exception as err:
            cls.tearDownClass()
            raise err

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._src_dir, True)
        shutil.rmtree(cls._dest_dir, True)

    def test_update_and_compile(self):
        updater = GitUpdate(self._src_dir)
        compiler = Compile(self._src_dir, self._dest_dir)

        self.assertNotEqual(self._head_commit, self.head_commit_hash(),
                            "Git reset HEAD~ worked")

        updater.run(None)
        self.assertEqual(self._head_commit, self.head_commit_hash(),
                         "Update moved HEAD to correct commit")

        with Path(self._src_dir, "scripts/HHPaths.pm").open("r") as hhpath:
            final_line = hhpath.readlines()[-1].strip()
            self.assertEqual(final_line, "# Comment from phyre engine tests",
                             "Stashed line remains")

        compiler.run(None)
        self.assertTrue(Path(self._dest_dir, "bin/hhblits").exists(),
                        "Generated binaries in correct location.")

    @classmethod
    def head_commit_hash(cls):
        result = subprocess.run(
            ["git",
             "--git-dir", str(Path(cls._src_dir, ".git")),
             "--work-tree", cls._src_dir,
             "rev-list", "--max-count=1", "HEAD"],
            stdout=subprocess.PIPE, check=True)
        return result.stdout.decode(sys.stdout.encoding).strip()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
