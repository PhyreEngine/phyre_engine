"""Test components in :py:mod:`phyre_engine.test.component.jmespath`."""
import unittest

import phyre_engine.component.jmespath as jmes
from copy import deepcopy
import jmespath
import datetime


SAMPLE_STATE = {
    "TM": 0.5,
    "templates": [
        {"TM": 0.3, "name": "foo"},
        {"TM": 0.6, "name": "bar"},
    ],
    "metadata": {"date": "1970-01-01"}
}

class TestJMESExtensions(unittest.TestCase):
    """Test custom JMESPath extensions."""

    def search(self, expr, data):
        """Call :py:func:`jmespath.search` with extended functions."""
        extensions = jmes.JMESExtensions(data)
        jmespath_opts = jmespath.Options(custom_functions=extensions)
        return jmespath.search(expr, data, jmespath_opts)

    def test_root(self):
        """Access root node via ``root`` function."""
        self.assertEqual(
            self.search("templates[].{root: root().TM}", SAMPLE_STATE),
            [{"root": 0.5}, {"root": 0.5}])

    def test_toordinal(self):
        """Test ``toordinal`` function on a date."""
        self.assertEqual(
            self.search("toordinal(d)", {"d": datetime.date(1970, 1, 1)}),
            719163)

    def test_start_stop(self):
        """Test ``start`` and ``stop`` functions on a range."""
        self.assertEqual(
            self.search("{start: start(r), stop: stop(r)}", {"r": range(10)}),
            {"start": 0, "stop": 10})

    def test_list(self):
        """Test ``list`` function can convert tuple to list."""
        self.assertEqual(
            self.search("{list: list(t)}", {"t": (0, 1, 2)}),
            {"list": [0, 1, 2]})

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
