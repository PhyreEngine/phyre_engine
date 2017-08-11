import unittest
from phyre_engine.component import Component
from phyre_engine.component.component import Map
import phyre_engine.pipeline

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

class TestMap(unittest.TestCase):
    """Test the phyre_engine.component.Map class."""

    class Double(Component):
        """Dummy component that multiplies a number by 2."""
        ADDS = []
        REQUIRED = ["num"]
        REMOVES = []

        def run(self, data, config=None, pipeline=None):
            data["num"] *= 2
            return data

    def setUp(self):
        """Create some test data to map."""
        self.state = {
            "values": [
                {"num": 1},
                {"num": 2},
                {"num": 3},
            ]
        }
        self.expected_state = {
            "values": [
                {"num": 2},
                {"num": 4},
                {"num": 6},
            ]
        }

    def test_map_pipeline(self):
        """Run Map with Pipeline as input."""
        components = [self.Double()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline)

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)

    def test_map_load(self):
        """Run Map with Dict as input."""
        components = [
            (self.Double.__module__, self.Double.__qualname__)
        ]
        map_cmpt = Map("values", {"components": components})

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)
