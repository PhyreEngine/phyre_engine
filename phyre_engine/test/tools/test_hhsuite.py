import unittest
import subprocess

from phyre_engine.tools.hhsuite import HHBlits

class TestHHSuite(unittest.TestCase):
    """Test common methods for hhsuite tools."""

    def test_cmd_line(self):
        """Make sure we can mix long and short options."""
        self.assertListEqual(
                HHBlits(database="test").command_line,
                HHBlits(d="test").command_line)

    def test_run(self):
        """Try running hhblits. This should fail."""
        with self.assertRaises(subprocess.CalledProcessError):
            HHBlits().run()
