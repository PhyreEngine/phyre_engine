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

class TestMove(TestAlteration):
    """Test Move component."""

    def test_move(self):
        """Use Move component to rename a field."""
        move = alter.Move("baz", "qux")
        result = move.run(self.pipe_state)
        self.assertDictEqual(
            result,
            {"foo": "bar", "qux": [{"frob": 1}, 2, 3]})

class TestSet(TestAlteration):
    """Test the Set component."""

    def test_set(self):
        """Set a component in the pipeline state."""
        set_cpt = alter.Set("xyz", "123")
        result = set_cpt.run(self.pipe_state)
        self.assertEqual(result["xyz"], "123")

    def test_format(self):
        """Format a value with Set."""
        set_cpt = alter.Set("xyz", "{foo}", reformat=True)
        result = set_cpt.run(self.pipe_state)
        self.assertEqual(result["xyz"], "bar")
