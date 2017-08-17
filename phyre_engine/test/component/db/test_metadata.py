import unittest
import io
from phyre_engine.component.db.metadata import ParseField
import copy

class TestParseField(unittest.TestCase):
    """Test ParseField component."""

    _TEMPLATES = [
        {"haystack": "1foo"},
        {"haystack": "2bar"},
    ]
    _STATE = {"templates": _TEMPLATES}

    def setUp(self):
        """Set 'state' var."""
        self.state = copy.deepcopy(self._STATE)

    def test_matcher(self):
        """Test matching arbitrary field."""
        parser = ParseField("haystack", r"^(?P<num>\d)(?P<char>\w{3})$")
        results = parser.run(self.state)
        self.assertDictEqual(
            results["templates"][0],
            {"haystack": "1foo", "num": "1", "char": "foo"})
        self.assertDictEqual(
            results["templates"][1],
            {"haystack": "2bar", "num": "2", "char": "bar"})


    def test_failure(self):
        """Test that an exception is raised on failure if must_match."""
        parser = ParseField("haystack", "^XXX", must_match=True)
        with self.assertRaises(ValueError, msg="must_match raises exception"):
            parser.run(copy.deepcopy(self.state))

    def test_no_match(self):
        """If not must_match, match failures shouldn't have an effect."""
        parser = ParseField("haystack", "^XXX", must_match=False)
        results = parser.run(copy.deepcopy(self.state))
        for i, template in enumerate(results["templates"]):
            with self.subTest("No attributes added", i=i):
                self.assertListEqual(["haystack"], list(template.keys()))

    def test_unicode_matching(self):
        """Check that the unicode_matching flag works."""
        state = copy.deepcopy(self.state)
        state["templates"][0]["haystack"] = "αβ: and so on"
        parser = ParseField("haystack", "^(?P<id>\w+)", unicode_matching=True)
        results = parser.run(state)
        self.assertEqual(
            results["templates"][0]["id"], "αβ",
            "Non-ASCII matching")

if __name__ == "__main__":
    unittest.main()
