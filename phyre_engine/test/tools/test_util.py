import unittest
import phyre_engine.tools.util as util
import os
import tempfile
import pathlib

class TestTemporaryEnvironment(unittest.TestCase):

    def test_alteration(self):
        path = os.environ["PATH"]
        with util.TemporaryEnvironment(PATH="fake_path"):
            self.assertEqual(
                os.environ["PATH"], "fake_path",
                "Environment was altered")
        self.assertEqual(
            os.environ["PATH"], path,
            "Environment was reverted")

class TestStream(unittest.TestCase):
    """Test phyre_engine.util.Stream"""
    _TEST_DATA = "test data\n"

    def setUp(self):
        """Create a temporary file containing 'test data\n'."""
        self.temp_file = tempfile.NamedTemporaryFile("w+")
        print(self._TEST_DATA, file=self.temp_file, end="")
        self.temp_file.flush()
        self.temp_file.seek(0)

    def tearDown(self):
        """Remove temporary file."""
        self.temp_file.close()

    def _verify_content(self, stream, msg=None):
        """Check that the content of the read file is correct."""
        self.assertListEqual(stream.readlines(), [self._TEST_DATA], msg)

    def test_filename(self):
        """Open a file name."""
        with util.Stream(self.temp_file.name, "r") as tmp_in:
            self._verify_content(tmp_in, "Read from file name")

    def test_stream(self):
        """Read from a stream."""
        with util.Stream(self.temp_file, "r") as tmp_in:
            self._verify_content(tmp_in, "Read from stream")

    def test_path(self):
        """Read from a pathlib.Path object."""
        path = pathlib.Path(self.temp_file.name)
        with util.Stream(path, "r") as tmp_in:
            self._verify_content(tmp_in, "Read from path")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
