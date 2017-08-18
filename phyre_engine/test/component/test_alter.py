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
