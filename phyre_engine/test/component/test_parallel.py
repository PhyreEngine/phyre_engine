import os
import unittest
import phyre_engine
from phyre_engine.component.parallel import ParallelComponent
from phyre_engine.component import Component

top_dir = os.path.dirname(os.path.dirname(phyre_engine.__file__))

class MockComponent(Component):
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def run(self, data):
        data["out"] = [x**2 for x in data["in"]]
        return data

class TestParallelComponent(unittest.TestCase):

    @unittest.skipUnless("QSUB_SCRATCH" in os.environ, "QSUB_SCRATCH undefined")
    def test_valid(self):
        """Attempt to run a simple parallel task."""
        scratch_dir = os.environ["QSUB_SCRATCH"]
        task = MockComponent()
        runner = ParallelComponent(task, 5, scratch_dir,
            "in", "out", path_dirs = [top_dir])

        data = {"in": list(range(0, 20))}
        results = runner.run(data)
        self.assertIn("out", results, "'out' key added")
        self.assertListEqual(
           [x**2 for x in data["in"]],
           results["out"],
           "Results are as expected")


if __name__ == "__main__":
    unittest.main()
