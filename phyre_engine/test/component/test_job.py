"""Test the phyre_engine.component.job module."""
import os
from pathlib import Path
import tempfile
import unittest
import phyre_engine.component.job as job

class TestDirectory(unittest.TestCase):
    """Test the Directory component."""

    def setUp(self):
        """Create a temporary base directory."""
        self.base_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        """Remove temporary base directory."""
        self.base_dir.cleanup()

    def test_create_dir(self):
        """This should create a new directory."""
        dir_component = job.Directory(self.base_dir.name, "foo-", "-bar", False)
        results = dir_component.run({})
        self.assertIn("job_directory", results, "Created job_directory")
        job_dir = results["job_directory"]

        # Check prefix and suffix were applied
        self.assertTrue(Path(job_dir).name.startswith("foo-"), "prefix")
        self.assertTrue(Path(job_dir).name.endswith("-bar"), "suffix")

    def test_chdir(self):
        """Supplying the chdir parameter should change directory."""
        try:
            cwd = os.getcwd()
            dir_component = job.Directory(self.base_dir.name, chdir=True)
            results = dir_component.run({})
            self.assertEqual(
                os.getcwd(), results["job_directory"],
                "Changed directory")
        finally:
            os.chdir(cwd)

    def test_perms(self):
        """Set permissions of the job directory."""
        dir_component = job.Directory(self.base_dir.name, perms=0o755)
        results = dir_component.run({})
        self.assertEqual(
            os.stat(results["job_directory"]).st_mode & 0o777, 0o755,
            "Set permissions")
