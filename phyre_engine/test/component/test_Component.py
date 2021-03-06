import copy
import io
import logging
import unittest
import unittest.mock
import phyre_engine.pipeline
from phyre_engine.component import Component
from phyre_engine.component.component import (Map, Conditional, TryCatch,
                                              PipelineComponent, Branch,
                                              ConfigLoader)
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

    def test_config(self):
        pipeline_config = phyre_engine.pipeline.PipelineConfig(
            {"conf": {"a": 1, "b": 2, "c": {"x": 1}}})
        component_params = {"a": 2, "c": {"x": 2, "y": 3}}

        # CONFIG_SECTION not defined
        self.assertEqual(
            Double.config(component_params, pipeline_config),
            component_params)

        with unittest.mock.patch.object(Double, "CONFIG_SECTION", new="conf"):
            # CONFIG_SECTION defined, so "a" is overridden.
            self.assertEqual(
                Double.config(component_params, pipeline_config),
                {"a": 2, "b": 2, "c": {"x": 2, "y": 3}})

class TestPipelineComponent(unittest.TestCase):
    """Test PipelineComponent."""

    class SubPipeline(PipelineComponent):
        """Simple pipeline component that returns nothing."""
        ADDS = []
        REMOVES = []
        REQUIRED = []
        def run(self, data, config=None, pipeline=None):
            pass

    _STATIC_CONFIG = {"foo": "bar"}
    _EMPTY_PIPELINE = {"components": [], "config": _STATIC_CONFIG}

    def test_prefer_child_config(self):
        """Update runtime config with static config."""
        mode = self.SubPipeline.ConfigurationPreference.PREFER_CHILD
        sub_pipe_cpt = self.SubPipeline(self._EMPTY_PIPELINE, config_mode=mode)

        # Non-conflicting fields added
        sub_pipe = sub_pipe_cpt.pipeline({"baz": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "bar", "baz": "qux"})

        # Static config wins conflicts
        sub_pipe = sub_pipe_cpt.pipeline({"foo": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "bar"})

    def test_prefer_parent_config(self):
        """Update static config with parent configuration."""
        mode = self.SubPipeline.ConfigurationPreference.PREFER_PARENT
        sub_pipe_cpt = self.SubPipeline(self._EMPTY_PIPELINE, config_mode=mode)

        # Non-conflicting fields added
        sub_pipe = sub_pipe_cpt.pipeline({"baz": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "bar", "baz": "qux"})

        # Runtime config wins conflicts
        sub_pipe = sub_pipe_cpt.pipeline({"foo": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "qux"})

    def test_discard_parent_config(self):
        """Runtime configuration is discarded."""
        mode = self.SubPipeline.ConfigurationPreference.DISCARD_PARENT
        sub_pipe_cpt = self.SubPipeline(self._EMPTY_PIPELINE, config_mode=mode)

        # Runtime config completely ignored
        sub_pipe = sub_pipe_cpt.pipeline({"baz": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "bar"})

    def test_no_lazy_load(self):
        """Lazily-loaded pipelines always use static config."""
        mode = self.SubPipeline.ConfigurationPreference.PREFER_PARENT
        sub_pipe_cpt = self.SubPipeline(self._EMPTY_PIPELINE, config_mode=mode,
                                        lazy_load=False)

        # Runtime config ignored
        sub_pipe = sub_pipe_cpt.pipeline({"baz": "qux"})
        self.assertDictEqual(sub_pipe.config, {"foo": "bar"})

    def test_load_components(self):
        """Load components from `components` parameter."""
        sub_pipe_cpt = self.SubPipeline(components=[Double.qualname])
        self.assertIsInstance(sub_pipe_cpt.pipeline({}).components[0], Double)

    def test_accept_pipeline_or_components(self):
        """Only a single one of `pipeline` or `components` allowed."""
        with self.assertRaises(ValueError):
            self.SubPipeline(pipeline=self._EMPTY_PIPELINE, components=[])
        with self.assertRaises(ValueError):
            self.SubPipeline()


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

    class _CopyComponent(Component):
        ADDS = ["result"]
        REMOVES = []
        REQUIRED = ["copy"]
        def run(self, data, config=None, pipeline=None):
            data["result"] = data["copy"]
            return data

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

    def test_discard(self):
        """Passing "discard" flag discards pipeline state."""
        components = [Double()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline, discard=True)

        self.expected_state["values"] = []
        results = map_cmpt.run(self.state)
        self.assertDictEqual(results, self.expected_state)

    def test_copy_list(self):
        """Copy items from the root into child elements."""
        components = [self._CopyComponent()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        map_cmpt = Map("values", pipeline, copy=["copy"])

        # CopyComponent will set result = copy for each element.
        state = {"copy": 1, "values": [{}, {}]}
        results = map_cmpt.run(state)
        self.assertDictEqual(
            results,
            {"copy": 1, "values": [{"result": 1}, {"result": 1}]})

    def test_copy_list(self):
        """Copy items from the root into child elements."""
        components = [self._CopyComponent()]
        pipeline = phyre_engine.pipeline.Pipeline(components)
        # Rename "foo" to "copy"
        map_cmpt = Map("values", pipeline, copy={"foo": "copy"})

        # CopyComponent will set result = copy for each element.
        state = {"foo": 1, "values": [{}, {}]}
        results = map_cmpt.run(state)
        self.assertDictEqual(
            results,
            {"foo": 1, "values": [{"result": 1}, {"result": 1}]})


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

    def test_branch_deep_copy(self):
        """Branched pipeline does not alter data in main branch."""
        pipe = phyre_engine.pipeline.Pipeline([self._AlteringComponent()])
        branch = Branch(pipe, shallow=False)
        start_pipe = {"foo": {"bar": "qux"}}
        results = branch.run(start_pipe)
        self.assertEqual(results, {"foo": {"bar": "qux"}})
        self.assertIs(results, start_pipe)

    def test_branch_shallow_copy(self):
        """Shallow copy does allow component to alter deep state."""
        pipe = phyre_engine.pipeline.Pipeline([self._AlteringComponent()])
        branch = Branch(pipe, shallow=True)
        start_pipe = {"foo": {"bar": "qux"}}
        results = branch.run(start_pipe)
        self.assertEqual(results, {"foo": {"bar": "baz"}})
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


class TestConfigLoader(unittest.TestCase):
    """Test ConfigLoader component."""

    def test_config_transfer(self):
        """Check that values are transferred from state to config."""
        conf_loader = ConfigLoader(
            mapping={"slice": "jmespath.value_expr"},
            components=[])
        self.assertEqual(
            conf_loader.generate_config({"slice": "foo"}, None),
            {"jmespath": {"value_expr": "foo"}})
