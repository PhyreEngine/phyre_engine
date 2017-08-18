"""Tests for phyre_engine.component.dump module."""

import copy
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
        self.pipe_input = copy.deepcopy(_PIPELINE_INPUT)
        self.expected = copy.deepcopy(_PIPELINE_OUTPUT)

    def tearDown(self):
        """Close in-memory stream."""
        self.stream.close()

    def _verify_dump(self, parsed):
        return self.assertDictEqual(
            parsed, self.expected,
            "Parsed dict and written dict equal")

class TestJsonDumper(DumperTestBase):
    """Test JSON dumper."""

    def _verify_dump(self):
        """Verify JSON dump."""
        self.stream.seek(0)
        return super()._verify_dump(json.load(self.stream))

    def test_json(self):
        """Dump data in JSON format."""
        json_dumper = dump.Json(self.stream)
        json_dumper.run(self.pipe_input)
        self._verify_dump()

    def test_exclude(self):
        """Test exclusions."""
        json_dumper = dump.Json(self.stream, exclude="^.{3}$")
        json_dumper.run(self.pipe_input)
        del self.expected["foo"]
        del self.expected["baz"]
        del self.expected["seq"]
        self._verify_dump()

    def test_include(self):
        """Exclude all and specifically include "seq"."""
        json_dumper = dump.Json(self.stream, exclude="^.*$", include="^seq$")
        json_dumper.run(self.pipe_input)
        del self.expected["foo"]
        del self.expected["baz"]
        del self.expected["range"]
        self._verify_dump()

    def test_range_exception(self):
        """An exception should be raised for a range with step != 1."""
        with self.assertRaises(TypeError):
            dump.Json(self.stream).run({"range": range(0, 10, 2)})

class TestYamlDumper(DumperTestBase):
    """Test YAML dumper."""

    def _verify_dump(self):
        """Verify YAML dump."""
        self.stream.seek(0)
        return super()._verify_dump(yaml.safe_load(self.stream))

    def test_yaml(self):
        """Dump data in YAML format."""
        yaml_dumper = dump.Yaml(self.stream)
        yaml_dumper.run(self.pipe_input)
        self._verify_dump()

    def test_exclude(self):
        """Test exclusions."""
        yaml_dumper = dump.Yaml(self.stream, exclude="^.{3}$")
        yaml_dumper.run(self.pipe_input)
        del self.expected["foo"]
        del self.expected["baz"]
        del self.expected["seq"]
        self._verify_dump()

    def test_include(self):
        """Exclude all and specifically include "seq"."""
        yaml_dumper = dump.Yaml(self.stream, exclude="^.*$", include="^seq$")
        yaml_dumper.run(self.pipe_input)
        del self.expected["foo"]
        del self.expected["baz"]
        del self.expected["range"]
        self._verify_dump()

    def test_range_exception(self):
        """An exception should be raised for a range with step != 1."""
        with self.assertRaises(yaml.representer.RepresenterError):
            dump.Yaml(self.stream).run({"range": range(0, 10, 2)})
