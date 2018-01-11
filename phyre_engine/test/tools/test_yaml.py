"""Test wrapped yaml loader and dumper."""
import io
import textwrap
import unittest

import yaml.constructor
import yaml.representer

import phyre_engine.tools.yaml as custom_yaml

class TestYaml(unittest.TestCase):
    """Test wrapped yaml loader and dumper."""

    def test_load_tuple(self):
        """Tuples supported as dictionary keys."""
        yaml_doc = "[1, 2]: foo"
        buf = io.StringIO(yaml_doc)
        self.assertEqual(custom_yaml.load(buf), {(1, 2): "foo"})

    def test_safe_load(self):
        """Safe loader is used by default."""
        # Attempt to load an object.
        # This should give a ConstructorError because we used the safe loader.
        yaml_doc = "!!python/object:int {}"
        buf = io.StringIO(yaml_doc)
        with self.assertRaises(yaml.constructor.ConstructorError):
            custom_yaml.load(buf)

    def test_safe_dump(self):
        """Safe dumper is used by default."""
        # Try and dump an object. This should give a RepresenterError.
        with self.assertRaises(yaml.representer.RepresenterError):
            custom_yaml.dump(self)
