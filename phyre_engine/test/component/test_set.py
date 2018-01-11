"""Test components in :py:mod:`phyre_engine.component.set`."""
import copy
import unittest

import phyre_engine.component.set

SAMPLE_STATE = {
    "A": [
        {"value": 1},
        {"value": 2},
        # Extra parameter used to test ordering of sets.
        {"value": 3, "extra": 0xDEADBEEF},
    ],
    "B": [
        {"value": 3},
        {"value": 4},
        {"value": 5},
    ],
    "C": [
        {"value": 5},
        {"value": 6},
        {"value": 7},
    ],
    "D": [
        {"a": 1, "b": 1},
        {"a": 1, "b": 2},
    ],
    "E": [
        {"a": 1, "b": 2},
        {"a": 2, "b": 2},
    ]

}
class SetTestBase(unittest.TestCase):
    """Base class for testing set components."""

    def setUp(self):
        """Set self.state to a copy of SAMPLE_STATE."""
        self.state = copy.deepcopy(SAMPLE_STATE)

    def verify_state(self, result, expected_extra):
        """
        Check that the expected result is the original result with the
        specified extra field.
        """
        for key, expected in expected_extra.items():
            # Unordered list comparison
            got = result.pop(key)
            with self.subTest("Unordered equal lists", key=key):
                self.assertEqual(len(expected), len(got))
                for elem in expected:
                    self.assertIn(elem, got)
        self.assertEqual(result, SAMPLE_STATE)


class TestDifference(SetTestBase):
    """Test Difference component."""

    def test_difference(self):
        """Take the difference of A, B and C."""
        cmpt = phyre_engine.component.set.Difference(
            ["A", "B", "C"],
            "value",
            "result")
        self.verify_state(
            cmpt.run(self.state),
            {"result": [{"value": 1}, {"value": 2}]})

    def test_tuple_key(self):
        """Component works when using multiple keys."""
        cmpt = phyre_engine.component.set.Difference(
            ["D", "E"], "@.[a, b]", "result")
        self.verify_state(
            cmpt.run(self.state),
            {"result": [{"a": 1, "b": 1}]})

    def test_single_set(self):
        """Pass data through unchanged for a single set."""
        cmpt = phyre_engine.component.set.Difference(["A"], "value", "result")
        self.verify_state(cmpt.run(self.state), {"result": SAMPLE_STATE["A"]})

    def test_ordering(self):
        """First element with identical keys is retained."""
        cmpt = phyre_engine.component.set.Difference(["D"], "a", "result")
        self.verify_state(cmpt.run(self.state),
                          {"result": SAMPLE_STATE["D"][0:1]})


class TestUnion(SetTestBase):
    """Test Union component."""

    def test_union(self):
        """Take the union of A, B and C."""
        cmpt = phyre_engine.component.set.Union(
            ["A", "B", "C"],
            "value",
            "result")
        self.verify_state(
            cmpt.run(self.state), {
                "result": [
                    {"value": 1},
                    {"value": 2},
                    {"value": 3, "extra": 0xDEADBEEF},
                    {"value": 4},
                    {"value": 5},
                    {"value": 6},
                    {"value": 7},
                ]
            }
        )

    def test_tuple_key(self):
        """Component works when using multiple keys."""
        cmpt = phyre_engine.component.set.Union(
            ["D", "E"], "@.[a, b]", "result")
        self.verify_state(
            cmpt.run(self.state), {
                "result": [
                    {"a": 1, "b": 1},
                    {"a": 1, "b": 2},
                    {"a": 2, "b": 2},
                ]
            }
        )

    def test_single_set(self):
        """Pass data through unchanged for a single set."""
        cmpt = phyre_engine.component.set.Union(["A"], "value", "result")
        self.verify_state(cmpt.run(self.state), {"result": SAMPLE_STATE["A"]})

    def test_ordering(self):
        """First element with identical keys is retained."""
        cmpt = phyre_engine.component.set.Union(["D"], "a", "result")
        self.verify_state(cmpt.run(self.state),
                          {"result": SAMPLE_STATE["D"][0:1]})


class TestIntersection(SetTestBase):
    """Test Intersection component."""

    def test_intersect(self):
        """Take the intersection of A and B."""
        cmpt = phyre_engine.component.set.Intersection(
            ["A", "B"],
            "value",
            "result")
        self.verify_state(
            cmpt.run(self.state),
            {"result": [{"value": 3, "extra": 0xDEADBEEF}]})

    def test_tuple_key(self):
        """Component works when using multiple keys."""
        cmpt = phyre_engine.component.set.Intersection(
            ["D", "E"], "@.[a, b]", "result")
        self.verify_state(cmpt.run(self.state), {"result": [{"a": 1, "b": 2}]})

    def test_single_set(self):
        """Pass data through unchanged for a single set."""
        cmpt = phyre_engine.component.set.Intersection(
            ["A"], "value", "result")
        self.verify_state(cmpt.run(self.state), {"result": SAMPLE_STATE["A"]})

    def test_ordering(self):
        """First element with identical keys is retained."""
        cmpt = phyre_engine.component.set.Intersection(["D"], "a", "result")
        self.verify_state(cmpt.run(self.state),
                          {"result": SAMPLE_STATE["D"][0:1]})
