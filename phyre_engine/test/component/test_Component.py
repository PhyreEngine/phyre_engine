import copy
import io
import logging
import unittest
import unittest.mock
from phyre_engine.component import Component
from phyre_engine.component.component import (Map, Conditional, TryCatch,
                                              PipelineComponent, Branch)
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

class TestPipelineComponent(unittest.TestCase):
    """Test PipelineComponent."""

    class SubPipeline(PipelineComponent):
        """Simple pipeline component that returns nothing."""
        ADDS = []
        REMOVES = []
        REQUIRED = []
        def run(self, data, config=None, pipeline=None):
            pass

    _RUNTIME_CONFIG = {"foo": "bar"}
    _STATIC_CONFIG = {"baz": "qux"}

    def test_update_config(self):
        """Update runtime config with static config."""
        pipeline = phyre_engine.pipeline.Pipeline(
            [], config=self._STATIC_CONFIG)
        sub_pipe = self.SubPipeline(pipeline)

        self.assertDictEqual(
            sub_pipe.config(self._RUNTIME_CONFIG),
            {"foo": "bar", "baz": "qux"})

    def test_discard_config(self):
        """Discard runtime config, replacing with static config."""
        pipeline = phyre_engine.pipeline.Pipeline(
            [], config=self._STATIC_CONFIG)
        sub_pipe = self.SubPipeline(pipeline, discard_config=True)

        self.assertDictEqual(
            sub_pipe.config(self._RUNTIME_CONFIG),
            self._STATIC_CONFIG)


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

    class _ListComponent(Component):
        ADDS = []
        REMOVES = []
        REQUIRED = []
        def run(self, data, config=None, pipeline=None):
            return [1, 2, 3]

    def test_map_pipeline(self):
        """Run Map with Pipeline as input."""
        components = [Double()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline)

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)

    def test_map_load(self):
        """Run Map with Dict as input."""
        components = [Double.qualname]
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

    def test_extend_list(self):
        """Child pipelines may return multiple items in a list."""
        components = [self._ListComponent()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline)

        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, {"values": [1, 2, 3] * 3})


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
            # pass_through = False, so results is None
            self.assertIsNone(results)

        with capture as logger:
            trycatch = TryCatch(
                pipeline=raise_pipe, logger=logger, pass_through=True)
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

class TestBranch(unittest.TestCase):
    """Test Branch component."""

    class _AlteringComponent(Component):
        # Used to alter the pipeline state so we can verify that the main branch
        # does not change.
        ADDS = []
        REMOVES = []
        REQUIRED = []

        def run(self, data, config=None, pipeline=None):
            data["foo"]["bar"] = "baz"
            return data

    def test_branch_copy(self):
        """Branched pipeline does not alter data in main branch."""
        pipe = phyre_engine.pipeline.Pipeline([self._AlteringComponent()])
        branch = Branch(pipe)
        start_pipe = {"foo": {"bar": "qux"}}
        results = branch.run(start_pipe)
        self.assertEqual(results, {"foo": {"bar": "qux"}})
        self.assertIs(results, start_pipe)

    def test_branch_run(self):
        """Ensure that the sub-pipeline is actually run."""
        component = self._AlteringComponent()
        component.run = unittest.mock.MagicMock(return_value={})
        pipe = phyre_engine.pipeline.Pipeline([component])
        Branch(pipe).run({})
        component.run.assert_called_once()

    def test_branch_keep(self):
        """Keep certain fields from the branch."""
        component = self._AlteringComponent()
        pipe = phyre_engine.pipeline.Pipeline([component])
        start_pipe = {"foo": {"bar": "qux"}}
        results = Branch(pipe, keep=('foo',)).run(start_pipe)
        self.assertIn("foo", results)
        self.assertEqual(results["foo"], {"bar": "baz"})
