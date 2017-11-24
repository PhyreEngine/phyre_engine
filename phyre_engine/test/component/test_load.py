"""Test components in the :py:mod:`phyre_engine.component.load` module."""
import copy
import io
import textwrap
import unittest
import phyre_engine.component.load as load
import pickle

_INPUT_PIPELINE = {"qux": "frob"}

# Input pipeline updated from serialised data.
_EXPECTED_PIPELINE = {
    "qux": "frob",
    "foo": "bar",
    "baz": [1, 2, 3]
}

class LoaderTestBase(unittest.TestCase):
    """Base class for loader tests."""

    def setUp(self):
        """Set up StringIO from self.SOURCE and deep copy expected result."""
        if self.MODE == "r":
            self.stream = io.StringIO()
        elif self.MODE == "rb":
            self.stream = io.BytesIO()

        self.stream.write(self.SOURCE)
        self.stream.seek(0)
        self.expected = copy.deepcopy(_EXPECTED_PIPELINE)
        self.pipeline = copy.deepcopy(_INPUT_PIPELINE)

    def _verify_load(self, loaded):
        """Check that the loaded pipeline is equal to self.expected."""
        # We don't care about whether lists were loaded as tuples or lists, so
        # we can't just use assertDictEqual
        self.assertSequenceEqual(loaded.keys(), self.expected.keys())
        if "qux" in loaded:
            self.assertEqual(loaded["qux"], self.expected["qux"])
        if "foo" in loaded:
            self.assertEqual(loaded["foo"], self.expected["foo"])
        if "baz" in loaded:
            self.assertSequenceEqual(loaded["baz"], self.expected["baz"])

class TestYaml(LoaderTestBase):
    """Test Yaml loader."""


    MODE = "r"
    SOURCE = textwrap.dedent("""\
    foo: bar
    baz: [1, 2, 3]
    """)

    def test_load(self):
        """Load simple pipeline state."""
        loader = load.Yaml(self.stream)
        result = loader.run(self.pipeline)
        self._verify_load(result)

class TestJson(LoaderTestBase):
    """Test Json loader."""

    MODE = "r"
    SOURCE = """{"foo": "bar", "baz": [1, 2, 3]}"""

    def test_load(self):
        """Load simple pipeline state."""
        loader = load.Json(self.stream)
        result = loader.run(self.pipeline)
        self._verify_load(result)

class TestPickle(LoaderTestBase):
    """Test Pickle loader."""

    MODE = "rb"
    SOURCE = pickle.dumps({"foo": "bar", "baz": [1, 2, 3]})

    def test_load(self):
        """Load simple pipeline state from pickle."""
        loader = load.Pickle(self.stream)
        result = loader.run(self.pipeline)
        self._verify_load(result)
