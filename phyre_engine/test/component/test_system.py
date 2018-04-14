"""Test :py:mod:phyre_engine.component.system."""
import textwrap
import unittest

import phyre_engine.component.system as system

class TestPython(unittest.TestCase):
    """Test phyre_engine.component.syste.Python component."""

    SCRIPT = textwrap.dedent("""\
    data["result"] = data["input"] * 2
    """)

    def test_execution(self):
        """Execute arbitrary code."""
        py = system.Python(self.SCRIPT)
        result = py.run({"input": 1})
        self.assertEqual(result["result"], 2)

    def test_adds(self):
        """Test ADDS fields."""
        py = system.Python(self.SCRIPT, adds="['a', 'b', 'c']")
        self.assertEqual(py.ADDS, ["a", "b", "c"])

    def test_removes(self):
        """Test REMOVES fields."""
        py = system.Python(self.SCRIPT, removes="['a', 'b', 'c']")
        self.assertEqual(py.REMOVES, ["a", "b", "c"])

    def test_required(self):
        """Test REQUIRED fields."""
        py = system.Python(self.SCRIPT, required="['a', 'b', 'c']")
        self.assertEqual(py.REQUIRED, ["a", "b", "c"])
