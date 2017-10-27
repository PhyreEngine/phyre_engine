"""Tests for :py:mod:`phyre_engine.component.sort`."""
import unittest
import copy
from phyre_engine.component.sort import Sort, Shuffle

# Sample data to sort
PIPE_STATE = {
    "parent": {
        "templates": [
            {"score": (0.3, 0.1)},
            {"score": (0.1, 0.1)},
            {"score": (0.4, 0.2)},
            {"score": (0.1, 0.3)},
            {"score": (0.5, 0.5)},
        ],
    },
    "numbers": [3, 1, 4, 1, 5, 9, 2, 6],
}

class TestSort(unittest.TestCase):
    """Test the :py:class:`phyre_engine.component.Sort` class."""

    def test_default_params(self):
        """Default parameters sort simple list into ascending order."""
        cmpt = Sort(field="numbers")
        results = cmpt.run(copy.deepcopy(PIPE_STATE))
        self.assertListEqual(
            results["numbers"],
            sorted(PIPE_STATE["numbers"], reverse=False))

    def test_default_params_nested(self):
        """Sort a nested field with the default parameters."""
        pipe_state = copy.deepcopy(PIPE_STATE)
        pipe_state["parent"]["numbers"] = pipe_state["numbers"]
        del pipe_state["numbers"]

        cmpt = Sort(field=["parent", "numbers"])
        results = cmpt.run(pipe_state)
        self.assertListEqual(
            results["parent"]["numbers"],
            sorted(PIPE_STATE["numbers"], reverse=False))

    def test_sort_nested_key(self):
        """Sort a nested list by a nested key, breaking ties."""
        cmpt = Sort(field=["parent", "templates"],
                    keys=[["ascending", "score", 0],
                          ["descending", "score", 1]])
        results = cmpt.run(copy.deepcopy(PIPE_STATE))
        self.assertListEqual(
            results["parent"]["templates"], [
                {"score": (0.1, 0.3)},
                {"score": (0.1, 0.1)},
                {"score": (0.3, 0.1)},
                {"score": (0.4, 0.2)},
                {"score": (0.5, 0.5)},
            ])


class TestShuffle(unittest.TestCase):
    """Test the :py:class:`phyre_engine.component.Shuffle` class."""
    # Note that we set the random seed for each Shuffle component to 1.

    def test_default_params(self):
        """Default parameters sort simple list into ascending order."""
        cmpt = Shuffle(field="numbers", seed=1)
        results = cmpt.run(copy.deepcopy(PIPE_STATE))
        self.assertListEqual(
            results["numbers"],
            [1, 2, 1, 9, 6, 3, 5, 4])

    def test_shuffle_nested_key(self):
        """Shuffle a nested list by a nested keys."""
        cmpt = Shuffle(field=["parent", "templates"], seed=1)
        results = cmpt.run(copy.deepcopy(PIPE_STATE))
        self.assertListEqual(
            results["parent"]["templates"], [
                {"score": (0.4, 0.2)},
                {"score": (0.1, 0.3)},
                {"score": (0.5, 0.5)},
                {"score": (0.3, 0.1)},
                {"score": (0.1, 0.1)},
            ])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
