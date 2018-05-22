"""Test wrapped yaml loader and dumper."""
import io
import textwrap
import unittest

import yaml.constructor
import yaml.representer

import phyre_engine.pipeline
import phyre_engine.tools.yaml as custom_yaml

class TestYaml(unittest.TestCase):
    """Test wrapped yaml loader and dumper."""

    def test_load_tuple(self):
        """Tuples supported as dictionary keys."""
        yaml_doc = "[1, 2]: foo"
        buf = io.StringIO(yaml_doc)
        self.assertEqual(custom_yaml.load(buf), {(1, 2): "foo"})

    def test_represent_pipeline_config(self):
        """PipelineConfig objects can be represented as dicts."""
        source = {"a": "b"}
        conf = phyre_engine.pipeline.PipelineConfig(source)
        buf = io.StringIO()
        custom_yaml.dump(conf, buf)
        buf.seek(0)

        result = custom_yaml.load(buf)
        self.assertEqual(result, source)

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

    def test_implicit_load_simple(self):
        """Test single-pass implicit loader."""
        buf = io.StringIO("{a: !template '{b}', b: x}")
        doc = yaml.load(buf, custom_yaml.ImplicitLoader)
        self.assertEqual(doc.resolve(doc.unresolved)["a"], "x")

    def test_implicit_load_failure(self):
        """Unresolved templates raise UnresolvedTemplateError."""
        buf = io.StringIO("{a: !template '{b}', b: !template '{c}'}")
        doc = yaml.load(buf, custom_yaml.ImplicitLoader)
        with self.assertRaises(custom_yaml.UnresolvedTemplateError):
            doc.resolve(doc.unresolved)

    def test_implicit_load_keyerror(self):
        """Unknown fields cause a KeyError."""
        buf = io.StringIO("{a: !template '{b}'}")
        doc = yaml.load(buf, custom_yaml.ImplicitLoader)
        with self.assertRaises(KeyError):
            doc.resolve(doc.unresolved)

    def test_implicit_load_two_pass(self):
        """Test two-pass implicit loader."""
        buf = io.StringIO("{a: !template '{b}', b: !template '{c}', c: x}")
        doc = yaml.load(buf, custom_yaml.ImplicitLoader)

        doc = doc.resolve(doc.unresolved, allow_unresolved=True)
        self.assertEqual(doc.unresolved["b"], "x")
        self.assertIsInstance(doc.unresolved["a"], custom_yaml.TemplateString)

        doc = doc.resolve(doc.unresolved, allow_unresolved=True)
        self.assertIsInstance(doc, dict)
        self.assertEqual(doc["a"], "x")
