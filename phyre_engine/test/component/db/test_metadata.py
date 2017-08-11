import unittest
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
import io
from phyre_engine.component.db.metadata import NameTemplate, ParseSequenceName, \
    ParseField
import copy

class TestNameTemplate(unittest.TestCase):
    """Test that NameTemplate can assign names correctly."""

    FASTA = ">FOO.BAR-BAZ-QUX\nAAAEEEAAA"

    def setUp(self):
        with io.StringIO(self.FASTA) as seq_fh:
            self.data = {
                "templates": [{"sequence": Bio.SeqIO.read(seq_fh, "fasta")}]}

    @staticmethod
    def result_name(data):
        return data["templates"][0]["name"]

    def test_builtin_full(self):
        """Using full_name builtin"""
        self.assertEqual(
            self.result_name(NameTemplate("full_name").run(self.data)),
            "FOO.BAR-BAZ-QUX")

    def test_builtin_split(self):
        """Using split builtin"""
        self.assertEqual(
            self.result_name(NameTemplate("split", ".", 0).run(self.data)),
            "FOO")
        self.assertEqual(
            self.result_name(NameTemplate("split", ".", 1).run(self.data)),
            "BAR-BAZ-QUX")

    def test_custom(self):
        """Using custom function"""
        custom_fn = lambda t: "custom"
        self.assertEqual(
            self.result_name(NameTemplate(custom_fn).run(self.data)),
            "custom")

class TestParseSequenceName(unittest.TestCase):
    """Test ParseSequenceName component."""

    _TEMPLATES = [
        {"sequence": SeqRecord(None, name="1ABC_X extra data")},
        {"sequence": SeqRecord(None, name="1ABC_Y ignore me")},
    ]
    _STATE = {"templates": _TEMPLATES}

    def _test_attributes(self, templates, atts):
        """Test all attributes match for each template."""
        for i, (template, att_dict) in enumerate(zip(templates, atts)):
            for key, value in att_dict.items():
                with self.subTest("Testing attribute", attribute=key, i=i):
                    self.assertIn(key, template, "Key is present")
                    self.assertEqual(template[key], value, "Values match")

    def test_full_name(self):
        """Test the full name regex given as an example."""
        regex = "^(?P<name>.*)$"
        parser = ParseSequenceName(regex)
        results = parser.run(copy.deepcopy(self._STATE))

        self._test_attributes(results["templates"], [
            {"name": "1ABC_X extra data"},
            {"name": "1ABC_Y ignore me"}
        ])

    def test_pdb_id(self):
        """Test the regex for parsing PDB and chain ID given as an example."""
        regex = "^(?P<name>(?P<PDB>\w{4})_(?P<chain>\w+))"
        parser = ParseSequenceName(regex)
        results = parser.run(copy.deepcopy(self._STATE))

        self._test_attributes(results["templates"], [
            {"name": "1ABC_X", "PDB": "1ABC", "chain": "X"},
            {"name": "1ABC_Y", "PDB": "1ABC", "chain": "Y"}
        ])

    def test_failure(self):
        """Test that an exception is raised on failure if must_match."""
        parser = ParseSequenceName("^XXX", must_match=True)
        with self.assertRaises(ValueError, msg="must_match raises exception"):
            parser.run(copy.deepcopy(self._STATE))

    def test_no_match(self):
        """If not must_match, match failures shouldn't have an effect."""
        parser = ParseSequenceName("^XXX", must_match=False)
        results = parser.run(copy.deepcopy(self._STATE))
        for i, template in enumerate(results["templates"]):
            with self.subTest("No attributes added", i=i):
                self.assertListEqual(["sequence"], list(template.keys()))

    def test_unicode_matching(self):
        """Check that the unicode_matching flag works."""
        state = copy.deepcopy(self._STATE)
        state["templates"][0]["sequence"].name = "αβ: and so on"
        parser = ParseSequenceName("^(?P<id>\w+)", unicode_matching=True)
        results = parser.run(state)
        self.assertEqual(
            results["templates"][0]["id"], "αβ",
            "Non-ASCII matching")

class TestParseField(unittest.TestCase):
    """Test ParseField component."""

    _TEMPLATES = [
        {"haystack": "1foo"},
        {"haystack": "2bar"},
    ]
    _STATE = {"templates": _TEMPLATES}

    def test_matcher(self):
        """Test matching arbitrary field."""
        parser = ParseField("haystack", r"^(?P<num>\d)(?P<char>\w{3})$")
        results = parser.run(self._STATE)
        self.assertDictEqual(
            results["templates"][0],
            {"haystack": "1foo", "num": "1", "char": "foo"})
        self.assertDictEqual(
            results["templates"][1],
            {"haystack": "2bar", "num": "2", "char": "bar"})


if __name__ == "__main__":
    unittest.main()
