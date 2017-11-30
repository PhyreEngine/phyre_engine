"""Tests for tools in :py:mod:`phyre_engine.tools.jmespath."""
import datetime
import unittest

import jmespath

from phyre_engine.tools.jmespath import JMESExtensions

class TestJMESExtensions(unittest.TestCase):
    """Test custom JMESPath extensions."""

    def search(self, expr, data):
        """Call :py:func:`jmespath.search` with extended functions."""
        extensions = JMESExtensions(data)
        jmespath_opts = jmespath.Options(custom_functions=extensions)
        return jmespath.search(expr, data, jmespath_opts)

    def test_root(self):
        """Access root node via ``root`` function."""
        sample_state = {"top": 0.5, "list": [{}, {}]}
        self.assertEqual(
            self.search("list[].{root: root().top}", sample_state),
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

    def test_delete(self):
        """Test ``delete`` function can remove fields."""
        self.assertEqual(
            self.search(
                """delete(@, `["a", "b"]`)""",
                {"a": 1, "b": 2, "c": 3}),
            {"c": 3})
