import unittest
import Bio.SeqIO
import io
from phyre_engine.component.db.metadata import NameTemplate

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

if __name__ == "__main__":
    unittest.main()
