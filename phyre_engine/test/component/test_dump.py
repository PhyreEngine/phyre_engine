"""Tests for phyre_engine.component.dump module."""

import copy
import io
import unittest
import json
import yaml
import phyre_engine.component.dump as dump
import Bio.SeqRecord
import Bio.Seq
import csv

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
        # Can't just use assertDictEqual because we don't care whether lists
        # were deserialised as lists or tuples.
        self.assertSequenceEqual(parsed.keys(), self.expected.keys())

        if "foo" in parsed:
            self.assertEqual(parsed["foo"], self.expected["foo"])
        if "seq" in parsed:
            self.assertEqual(parsed["seq"], self.expected["seq"])
        if "baz" in parsed:
            self.assertSequenceEqual(parsed["baz"], self.expected["baz"])
        if "range" in parsed:
            self.assertSequenceEqual(parsed["range"], self.expected["range"])

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

class TestCsvDumper(unittest.TestCase):
    """Test CSV dumper."""

    _SAMPLE_DATA = [
        {"a": 1, "b": 2, "c": 3},
        {"a": 1, "b": 2},
        {"a": 1, "b": 2, "c": None},
        {"a": 1, "b": 2, "c": 3},
        {"a": 1, "b": 2, "c": 3},
    ]
    _SAMPLE_PIPE = {"sample": _SAMPLE_DATA}

    _EXPECTED_ROWS = [
        ["a", "b", "c"],
        ["1", "2", "3"],
        ["1", "2", "MISSING"],
        ["1", "2", "NULL"],
        ["1", "2", "3"],
        ["1", "2", "3"],
    ]

    def roundtrip(self, *args, **kwargs):
        with io.StringIO("w+") as buffer:
            csv_dumper = dump.Csv(file=buffer, *args, **kwargs)
            pipeline = copy.deepcopy(self._SAMPLE_PIPE)
            csv_dumper.run(pipeline)
            buffer.seek(0)
            reader = csv.reader(buffer)
            rows = [row for row in reader]
        return rows

    def test_csv_default_fields(self):
        """Round-trip sample data using default fields."""
        results = self.roundtrip(
            "sample",
            null_placeholder="NULL", missing_placeholder="MISSING")
        self.assertEqual(results, self._EXPECTED_ROWS)

    def test_csv_explicit_fields(self):
        """Round-trip sample data using pre-selected fields."""
        results = self.roundtrip(
            "sample", ("b", "c"),
            null_placeholder="NULL", missing_placeholder="MISSING")
        expected = [row[1:] for row in self._EXPECTED_ROWS]
        self.assertEqual(results, expected)

    def test_csv_no_header(self):
        """Round-trip output with no header."""
        results = self.roundtrip(
            "sample", ("a", "b", "c"), header=False,
            null_placeholder="NULL", missing_placeholder="MISSING")
        self.assertEqual(results, self._EXPECTED_ROWS[1:])
