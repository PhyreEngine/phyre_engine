import unittest
from phyre_engine.tools.util import TemporaryEnvironment
import os

class TestTemporaryEnvironment(unittest.TestCase):

    def test_alteration(self):
        path = os.environ["PATH"]
        with TemporaryEnvironment(PATH="fake_path"):
            self.assertEqual(
                os.environ["PATH"], "fake_path",
                "Environment was altered")
        self.assertEqual(
            os.environ["PATH"], path,
            "Environment was reverted")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
