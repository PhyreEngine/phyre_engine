"""Tests components in the phyre_engine.component.filter module."""
import unittest
import phyre_engine.component.filter as filters

class TestFilters(unittest.TestCase):
    """Test Whitelist and Blacklist filter components."""
    _STATE = {"foo": "bar", "baz": "qux"}

    def test_whitelist(self):
        """Whitelist keeps the specified keys."""
        whitelist = filters.Whitelist("foo")
        new_state = whitelist.run(self._STATE.copy())
        self.assertDictEqual(new_state, {"foo": "bar"}, "Filtered baz")
