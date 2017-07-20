"""Tests for phyre_engine.component.dump module."""

import io
import unittest
import json
import yaml
import phyre_engine.component.dump as dump

class TestDumpers(unittest.TestCase):
    """Test all dumpers."""

    _TEST_PIPELINE = {"foo": "bar", "baz": [1, 2, 3]}

    def setUp(self):
        """Create an in-memory stream we can read and write from."""
        self.stream = io.StringIO()

    def tearDown(self):
        """Close in-memory stream."""
        self.stream.close()


    def _verify_dump(self, parsed):
        return self.assertDictEqual(
            parsed, self._TEST_PIPELINE,
            "Parsed dict and written dict equal")

    def test_json(self):
        """Dump data in JSON format."""
        json_dumper = dump.Json(self.stream)
        json_dumper.run(self._TEST_PIPELINE)
        self.stream.seek(0)
        self._verify_dump(json.load(self.stream))

    def test_yaml(self):
        """Dump data in YAML format."""
        yaml_dumper = dump.Yaml(self.stream)
        yaml_dumper.run(self._TEST_PIPELINE)
        self.stream.seek(0)
        self._verify_dump(yaml.safe_load(self.stream))
