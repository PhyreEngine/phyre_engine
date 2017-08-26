import copy
import unittest
from phyre_engine.component import Component
from phyre_engine.component.component import Map, Conditional
import phyre_engine.pipeline

class Double(Component):
    """Dummy component that multiplies a number by 2."""
    ADDS = []
    REQUIRED = ["num"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        data["num"] *= 2
        return data

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
        components = [Double()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline)

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)

    def test_map_load(self):
        """Run Map with Dict as input."""
        components = ["{}.{}".format(Double.__module__, Double.__qualname__)]
        map_cmpt = Map("values", {"components": components})

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)

class TestConditional(unittest.TestCase):
    """Test the phyre_engine.component.Conditional class."""

    def setUp(self):
        """Create some test data to test on."""
        self.states = copy.deepcopy([
            {"num": 1, "do_double": False},
            {"num": 2, "do_double": True}
        ])

        self.expected_states = copy.deepcopy([
            {"num": 1, "do_double": False},
            {"num": 4, "do_double": True}
        ])

    def test_conditional(self):
        """Run Conditional with Pipeline as input."""
        for start_state, end_state in zip(self.states, self.expected_states):
            components = [Double()]
            pipeline = phyre_engine.pipeline.Pipeline(components)
            cond_cmpt = Conditional("do_double", pipeline)

            results = cond_cmpt.run(start_state)
            self.assertDictEqual(results, end_state)
