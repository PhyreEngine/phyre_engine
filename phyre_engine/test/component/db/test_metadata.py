import unittest
import io
from phyre_engine.component.db.metadata import ParseField
import copy

class TestParseField(unittest.TestCase):
    """Test ParseField component."""
    _STATE = {"haystack": "1foo"}

    def setUp(self):
        """Set 'state' var."""
        self.state = copy.deepcopy(self._STATE)

    def test_matcher(self):
        """Test matching arbitrary field."""
        parser = ParseField("haystack", r"^(?P<num>\d)(?P<char>\w{3})$")
        results = parser.run(self.state)
        self.assertDictEqual(
            results,
            {"haystack": "1foo", "num": "1", "char": "foo"})

    def test_failure(self):
        """Test that an exception is raised on failure if must_match."""
        parser = ParseField("haystack", "^XXX", must_match=True)
        with self.assertRaises(ValueError, msg="must_match raises exception"):
            parser.run(self.state)

    def test_no_match(self):
        """If not must_match, match failures shouldn't have an effect."""
        parser = ParseField("haystack", "^XXX", must_match=False)
        results = parser.run(self.state)
        self.assertDictEqual(results, self._STATE)

    def test_unicode_matching(self):
        """Check that the unicode_matching flag works."""
        self.state["haystack"] = "αβ: and so on"
        parser = ParseField("haystack", "^(?P<id>\w+)", unicode_matching=True)
        results = parser.run(self.state)
        self.assertEqual(
            results["id"], "αβ",
            "Non-ASCII matching")

if __name__ == "__main__":
    unittest.main()
