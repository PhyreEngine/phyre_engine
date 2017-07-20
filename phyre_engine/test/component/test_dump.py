"""Tests for phyre_engine.component.dump module."""

import io
import unittest
import json
import yaml
import phyre_engine.component.dump as dump
import Bio.SeqRecord
import Bio.Seq

# This will be serialised...
_PIPELINE_INPUT = {
    "foo": "bar",
    "baz": [1, 2, 3],
    "seq": Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("AHHEG"),
        id="foo", description="desc"),
    "range": range(0, 100)
}

# Into this:
_PIPELINE_OUTPUT = {
    "foo": "bar",
    "baz": [1, 2, 3],
    "seq": ">foo desc\nAHHEG\n",
    "range": [0, 100]
}

class DumperTestBase(unittest.TestCase):
    """Base class for testing dumpers."""

    def setUp(self):
        """Create an in-memory stream we can read and write from."""
        self.stream = io.StringIO()

    def tearDown(self):
        """Close in-memory stream."""
        self.stream.close()

    def _verify_dump(self, parsed):
        return self.assertDictEqual(
            parsed, _PIPELINE_OUTPUT,
            "Parsed dict and written dict equal")

class TestJsonDumper(DumperTestBase):
    """Test JSON dumper."""

    def test_json(self):
        """Dump data in JSON format."""
        json_dumper = dump.Json(self.stream)
        json_dumper.run(_PIPELINE_INPUT)
        self.stream.seek(0)
        self._verify_dump(json.load(self.stream))

    def test_range_exception(self):
        """An exception should be raised for a range with step != 1."""
        with self.assertRaises(TypeError):
            dump.Json(self.stream).run({"range": range(0, 10, 2)})

class TestYamlDumper(DumperTestBase):
    """Test YAML dumper."""

    def test_yaml(self):
        """Dump data in YAML format."""
        yaml_dumper = dump.Yaml(self.stream)
        yaml_dumper.run(_PIPELINE_INPUT)
        self.stream.seek(0)
        self._verify_dump(yaml.safe_load(self.stream))

    def test_range_exception(self):
        """An exception should be raised for a range with step != 1."""
        with self.assertRaises(yaml.representer.RepresenterError):
            dump.Yaml(self.stream).run({"range": range(0, 10, 2)})
