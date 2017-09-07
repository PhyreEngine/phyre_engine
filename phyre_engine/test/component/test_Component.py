import copy
import io
import logging
import unittest
from phyre_engine.component import Component
from phyre_engine.component.component import Map, Conditional, TryCatch
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

    class _NoneComponent(Component):
        ADDS = []
        REMOVES = []
        REQUIRED = []
        def run(self, data, config=None, pipeline=None):
            return None

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

    def test_exclude_none(self):
        """Map should discard results when the child pipeline returns None."""
        components = [self._NoneComponent()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline)

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, {"values": []})


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

class TestTryCatch(unittest.TestCase):
    """Test the TryCatch component."""

    class _RaisingComponent(Component):
        """Raises an exception when run."""
        REQUIRED = []
        ADDS = []
        REMOVES = []

        def run(self, data, config=None, pipeline=None):
            raise RuntimeError("Deliberately raising an exception.")

    class _LogCapture:
        """Simple context manager for capturing log output."""
        def __init__(self):
            self.log_buffer = io.StringIO()
            logger = logging.getLogger("test_logger")
            logger.setLevel(logging.DEBUG)
            handler = logging.StreamHandler(self.log_buffer)
            formatter = logging.Formatter("%(levelname)s")
            handler.setFormatter(formatter)
            handler.setLevel(logging.DEBUG)
            logger.addHandler(handler)
            self.logger = logger

        def __enter__(self):
            return self.logger

        def __exit__(self, _exc_type, _exc_value, _traceback):
            self.log_buffer.seek(0)
            self.log_contents = self.log_buffer.read()

        def __str__(self):
            return self.log_contents

    def test_raise(self):
        """Ensure that our dummy component correctly raises an error."""
        raise_pipe = phyre_engine.pipeline.Pipeline([self._RaisingComponent()])
        with self.assertRaises(RuntimeError):
            raise_pipe.run()

    def test_trycatch(self):
        """Squash an exception with the TryCatch component."""
        raise_pipe = phyre_engine.pipeline.Pipeline([self._RaisingComponent()])
        initial_state = {"foo": "bar"}

        # Use a dummy logger to avoid false error messages on the console
        capture = self._LogCapture()
        with capture as logger:
            trycatch = TryCatch(pipeline=raise_pipe, logger=logger)
            results = trycatch.run(copy.deepcopy(initial_state))
        self.assertDictEqual(initial_state, results)

    def test_logging(self):
        """Test log contents after running TryCatch."""
        raise_pipe = phyre_engine.pipeline.Pipeline([self._RaisingComponent()])
        initial_state = {"foo": "bar"}

        # Use a dummy logger to avoid false error messages on the console
        for log_level in ("ERROR", "WARN", "INFO"):
            with self.subTest(log_level=log_level):
                capture = self._LogCapture()
                with capture as logger:
                    trycatch = TryCatch(
                        pipeline=raise_pipe,
                        logger=logger,
                        log_level=log_level)
                    trycatch.run(copy.deepcopy(initial_state))
                self.assertRegex(str(capture), "^" + log_level)
