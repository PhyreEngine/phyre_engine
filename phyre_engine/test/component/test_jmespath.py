"""Test components in :py:mod:`phyre_engine.test.component.jmespath`."""
import unittest

import phyre_engine.component.jmespath as jmes
from copy import deepcopy

SAMPLE_STATE = {
    "TM": 0.5,
    "templates": [
        {"TM": 0.3, "name": "foo"},
        {"TM": 0.6, "name": "bar"},
    ],
    "metadata": {"date": "1970-01-01"}
}

class TestUpdate(unittest.TestCase):
    """Test the Update module."""

    def test_update_list(self):
        """When value_expr gives a list, each element is updated."""
        update = jmes.Update("templates", "{qux: `123`}")
        results = update.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, {
            "TM": 0.5,
            "templates": [
                {"TM": 0.3, "name": "foo", "qux": 123},
                {"TM": 0.6, "name": "bar", "qux": 123},
            ],
            "metadata": {"date": "1970-01-01"}
        })

    def test_update_dict(self):
        """When value_expr gives a dict, it is updated."""
        update = jmes.Update("metadata", "{qux: `123`}")
        results = update.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, {
            "TM": 0.5,
            "templates": [
                {"TM": 0.3, "name": "foo"},
                {"TM": 0.6, "name": "bar"},
            ],
            "metadata": {"date": "1970-01-01", "qux": 123}
        })


class TestReplace(unittest.TestCase):
    """Test the Update module."""

    def test_update_list(self):
        """When value_expr gives a list, it is replaced in-place."""
        replace = jmes.Replace("templates", "@[?TM < `0.5`]")
        results = replace.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, {
            "TM": 0.5,
            "templates": [
                {"TM": 0.3, "name": "foo"},
            ],
            "metadata": {"date": "1970-01-01"}
        })

    def test_update_dict(self):
        """When value_expr gives a dict, it is replaced in-place."""
        replace = jmes.Replace("metadata", "{qux: `123`}")
        results = replace.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, {
            "TM": 0.5,
            "templates": [
                {"TM": 0.3, "name": "foo"},
                {"TM": 0.6, "name": "bar"},
            ],
            "metadata": {"qux": 123}
        })



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
