import unittest
from phyre_engine.component import Component

class TestComponent(unittest.TestCase):
    """Test non-abstract methods in the Component abstract base class."""

    class MockComponentScalar(Component):
        REQUIRED = ["a"]
        ADDS = []
        REMOVES = []
        def run(self, data, config=None, pipeline=None): pass

    class MockComponentList(Component):
        REQUIRED = ["a", "b"]
        ADDS = []
        REMOVES = []
        def run(self, data, config=None, pipeline=None): pass

    def test_get_vals(self):
        data = {"a":123, "b":456}

        a = TestComponent.MockComponentScalar().get_vals(data)
        self.assertEqual(a, 123)

        a, b = TestComponent.MockComponentList().get_vals(data)
        self.assertEqual(a, 123)
        self.assertEqual(b, 456)
