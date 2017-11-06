"""Tests components in the phyre_engine.component.filter module."""
import unittest
import phyre_engine.component.filter as filters
import copy
import datetime

class TestListFilters(unittest.TestCase):
    """Test Whitelist and Blacklist filter components."""
    _STATE = {"foo": "bar", "baz": "qux"}

    def test_whitelist(self):
        """Whitelist keeps the specified keys."""
        whitelist = filters.Whitelist("foo")
        new_state = whitelist.run(self._STATE.copy())
        self.assertDictEqual(new_state, {"foo": "bar"}, "Filtered baz")

    def test_blacklist(self):
        """Blacklist removes the specified keys."""
        blacklist = filters.Blacklist("foo")
        new_state = blacklist.run(self._STATE.copy())
        self.assertDictEqual(new_state, {"baz": "qux"}, "Filtered foo")

class TestFilter(unittest.TestCase):
    """Test the Filter component."""

    PIPELINE_STATE = {"TM": 0.5, "templates": [
        {"name": "foo", "TM": 0.3},
        {"name": "bar", "TM": 0.6},
    ]}

    def setUp(self):
        """Create `self.pipeline_state` that can be arbitrarily mutated."""
        self.pipeline_state = copy.deepcopy(self.PIPELINE_STATE)

    def test_filter(self):
        """Filter elements with ``TM < self.pipeline_state["TM"]``."""
        filt = filters.Filter("templates", "templates[?TM < root('TM')]")
        results = filt.run(self.pipeline_state)
        self.assertEqual(
            results["templates"],
            self.PIPELINE_STATE["templates"][0:1])

    def test_nested_filter(self):
        """Filter a nested list."""
        self.pipeline_state["top"] = {
            "templates": self.pipeline_state["templates"]}
        del self.pipeline_state["templates"]
        filt = filters.Filter("top.templates",
                              "top.templates[?TM < root('TM')]")
        results = filt.run(self.pipeline_state)
        self.assertEqual(
            results["top"]["templates"],
            self.PIPELINE_STATE["templates"][0:1])

    def test_toordinal(self):
        """Test filtering by "toordinal" function."""
        state = {
            "today": datetime.date(1990, 1, 1),
            "events": [
                {"what": "foo", "when": datetime.date(1990, 1, 1)},
                {"what": "bar", "when": datetime.date(1989, 12, 31)},
                {"what": "baz", "when": datetime.date(1990, 1, 2)},
            ]
        }
        filt = filters.Filter(
            "events", "events[?toordinal(when) > toordinal(root('today'))]")
        results = filt.run(state)
        self.assertEqual(
            results["events"],
            [{"what": "baz", "when": datetime.date(1990, 1, 2)}])
