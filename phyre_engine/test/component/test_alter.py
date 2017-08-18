"""Test the :py:mod:`phyre_engine.component.alter` module."""
import copy
import unittest
import phyre_engine.component.alter as alter

class TestAlteration(unittest.TestCase):
    """Test each component in the 'alter' module."""

    _PIPELINE = {
        "foo": "bar",
        "baz": [{"frob": 1}, 2, 3]
    }

    def setUp(self):
        """Make a copy of the sample pipeline."""
        self.pipe_state = copy.deepcopy(self._PIPELINE)

class TestRemove(TestAlteration):
    """Test Remove component."""

    def test_remove(self):
        """Use Remove component to remove a field from the pipeline state."""
        remove = alter.Remove("foo")
        self.assertDictEqual(
            remove.run(self.pipe_state),
            {"baz": [{"frob": 1}, 2, 3]})

class TestCopy(TestAlteration):
    """Test Copy component."""
    _FROM = "baz"
    _TO = "copy"
    _EXPECTED_STATE = {
        "foo": "bar",
        "baz": [{"frob": 1}, 2, 3],
        "copy": [{"frob": 1}, 2, 3]
    }

    def test_copy_shallow(self):
        """Use Copy component to shallow copy a field."""
        copy = alter.Copy(self._FROM, self._TO, deep=False)
        result = copy.run(self.pipe_state)
        self.assertDictEqual(result, self._EXPECTED_STATE)

        # Altering shallow copy should alter original state
        result["copy"][0]["frob"] = 4
        self.assertEqual(result["baz"][0]["frob"], 4)
        self.assertEqual(self.pipe_state["baz"][0]["frob"], 4)

    def test_copy_deep(self):
        """Use Copy component to deep copy a field."""
        copy = alter.Copy(self._FROM, self._TO, deep=True)
        result = copy.run(self.pipe_state)
        self.assertDictEqual(result, self._EXPECTED_STATE)

        # Altering deep copy should not alter original state
        result["copy"][0]["frob"] = 4
        self.assertEqual(result["baz"][0]["frob"], 1)
        self.assertEqual(self.pipe_state["baz"][0]["frob"], 1)
