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

    def test_update_identity(self):
        """Updating identical elements does nothing."""
        update = jmes.Update("@", "@")
        results = update.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, SAMPLE_STATE)

    def test_type_change(self):
        """Exception is raised if replacement type differs."""
        with self.assertRaises(TypeError):
            update = jmes.Update("metadata", "root().TM")
            update.run(deepcopy(SAMPLE_STATE))
        with self.assertRaises(TypeError):
            update = jmes.Update("templates", "root().TM")
            update.run(deepcopy(SAMPLE_STATE))

    def test_invalid_selection(self):
        """Can only update dicts and lists."""
        with self.assertRaises(TypeError):
            update = jmes.Update("TM", "`[1,2,3]`")
            update.run(deepcopy(SAMPLE_STATE))


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

    def test_replace_identity(self):
        """Replacing identical elements does nothing."""
        replace = jmes.Replace("@", "@")
        results = replace.run(deepcopy(SAMPLE_STATE))
        self.assertEqual(results, SAMPLE_STATE)

    def test_type_change(self):
        """Exception is raised if replacement type differs."""
        with self.assertRaises(TypeError):
            replace = jmes.Replace("metadata", "root().TM")
            replace.run(deepcopy(SAMPLE_STATE))
        with self.assertRaises(TypeError):
            replace = jmes.Replace("templates", "root().TM")
            replace.run(deepcopy(SAMPLE_STATE))

    def test_invalid_selection(self):
        """Can only replace dicts and lists."""
        with self.assertRaises(TypeError):
            replace = jmes.Replace("TM", "`[1,2,3]`")
            replace.run(deepcopy(SAMPLE_STATE))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
