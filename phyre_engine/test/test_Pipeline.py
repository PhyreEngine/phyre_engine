import unittest
from unittest.mock import MagicMock
import unittest.mock
from phyre_engine.component import Component
from phyre_engine import Pipeline
import tempfile
import json
import pickle
from phyre_engine.pipeline import Checkpoint
import pathlib

class MockNonNestedComponent(Component):
    """Non-nested component."""
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        return data

class TestPipeline(unittest.TestCase):
    """Check that we can assemble, validate and execute pipelines."""

    class MockComponentStart(Component):
        """Seed a couple of input elements."""
        REQUIRED = []
        ADDS = ["MCA_1", "MCA_2"]
        REMOVES = []

        def run(self, data, config=None, pipeline=None):
            data["MCA_1"] = "MCA_1"
            data["MCA_2"] = "MCA_2"
            return data

    class MockComponentMid(Component):
        """Removes a required key, so will fail if called twice."""
        REQUIRED = ["MCA_1"]
        ADDS = ["MCB_1"]
        REMOVES = ["MCA_1"]

        # Read a config section for this one
        CONFIG_SECTION = "mid"

        def __init__(self, *args, **kwargs):
            """Store args so they can be checked later."""
            self.args = args
            self.kwargs = kwargs

        def run(self, data, config=None, pipeline=None):
            del data["MCA_1"]
            data["MCB_1"] = "MCB_1"
            return data

    class MockComponentBad(Component):
        """Badly-behaving mock component that does not add a key."""
        REQUIRED = []
        ADDS = ["MCA_1"]
        REMOVES = []

        def run(self, data, config=None, pipeline=None):
            return data

    def test_valid(self):
        """Test a valid pipeline. All required attributes are present."""

        #No exception raised
        pipe = Pipeline([
            TestPipeline.MockComponentStart(),
            TestPipeline.MockComponentMid()
        ])
        pipe.validate()

        #Run the pipeline
        out = pipe.run()
        expected = {"MCA_2": "MCA_2", "MCB_1": "MCB_1"}
        self.assertDictEqual(expected, out)

    def test_invalid(self):
        """Test an invalid pipeline."""

        #The second MockComponentB requires MCA_1 but it was deleted.
        with(self.assertRaises(Pipeline.ValidationError)) as cm:
            Pipeline([
                TestPipeline.MockComponentStart(),
                TestPipeline.MockComponentMid(),
                TestPipeline.MockComponentMid()
            ]).validate()
        self.assertSetEqual(set(cm.exception.missing), set(["MCA_1"]))

    def test_init_values(self):
        """Test passing initialiser values to a pipeline."""

        #Should raise because MockComponentB requires MCA_1
        with(self.assertRaises(Pipeline.ValidationError)):
            Pipeline([TestPipeline.MockComponentMid()]).validate()

        #Should not raise because we pass an initial value
        Pipeline([TestPipeline.MockComponentMid()], {"MCA_1":10}).validate()


    def test_bad_pipeline(self):
        """Badly-behaving pipes should be caught at runtime."""

        #No exception raised because MockComponentBad promises to add "MCA_1"
        pipe = Pipeline([
            TestPipeline.MockComponentBad(),
            TestPipeline.MockComponentMid()
        ])
        pipe.validate()

        #Run the pipeline
        with(self.assertRaises(Pipeline.ValidationError)) as cm:
            pipe.run()
        self.assertSetEqual(set(cm.exception.missing), set(["MCA_1"]))
        self.assertDictEqual({}, cm.exception.data)

    def test_load_from_dict(self):
        """Load a pipeline from a dictionary of components."""

        dict_pipe = {
            "start": {"abc": 123, "xyz": 789},
            "components": [
                TestPipeline.MockComponentStart.qualname,
                {TestPipeline.MockComponentMid.qualname: [
                    "foo", "bar", {"baz": "qux"}]},
                MockNonNestedComponent.qualname
            ]
        }
        pipe = Pipeline.load(dict_pipe)
        self.assertDictEqual(
            pipe.start,
            dict_pipe["start"],
            "'Start' argument correctly passed to pipeline.")

        self.assertIsInstance(
            pipe.components[0],
            TestPipeline.MockComponentStart)
        self.assertIsInstance(
            pipe.components[1],
            TestPipeline.MockComponentMid)
        self.assertIsInstance(
            pipe.components[2],
            MockNonNestedComponent)
        self.assertTupleEqual(pipe.components[1].args, ("foo", "bar"))
        self.assertDictEqual(pipe.components[1].kwargs, {"baz": "qux"})

    @unittest.mock.patch("importlib.import_module")
    def test_default_namespace(self, import_mock):
        """Load a component from the default namespace."""
        dict_pipe = {"components": [".x.y.Foo"]}
        pipe = Pipeline.load(dict_pipe)
        import_mock.assert_called_once_with("phyre_engine.component.x.y")
        import_mock.return_value.Foo.assert_called_once_with()

    @unittest.mock.patch("importlib.import_module")
    def test_new_namespace(self, import_mock):
        """Load a component from a new namespace."""
        dict_pipe = {"namespace": "a.b", "components": [".x.y.Foo"]}
        pipe = Pipeline.load(dict_pipe)
        import_mock.assert_called_once_with("a.b.x.y")
        import_mock.return_value.Foo.assert_called_once_with()

    def test_checkpointing(self):
        """Test that we can save and load checkpoints."""

        start_cpt = TestPipeline.MockComponentStart()
        mid_cpt = TestPipeline.MockComponentMid()
        # Introduce a "bug" that we can try and fix.
        def crash(_data, _config, _pipeline):
            raise RuntimeError("Bad things")
        mid_cpt.run = crash

        with tempfile.TemporaryDirectory() as tempdir:
            check_file = pathlib.Path(tempdir, "checkpoint.chk")

            pipe = Pipeline([
                start_cpt,
                mid_cpt,
            ], checkpoint=check_file)

            try:
                pipe.run()
            except RuntimeError:
                pass

            with check_file.open("rb") as check_in:
                checkpoint = pickle.load(check_in)

            self.assertEqual(
                checkpoint.current_component, 1,
                "Checkpoint points to correct component")
            self.assertDictEqual(
                {"MCA_1": "MCA_1", "MCA_2": "MCA_2"},
                checkpoint.state,
                "Checkpoint contains pipeline state")

            # Alter state to verify component 0 wasn't re-run
            checkpoint = Checkpoint(1, {"MCA_1": "foo", "MCA_2": "bar"})
            with check_file.open("wb") as check_out:
                checkpoint = pickle.dump(checkpoint, check_out)

            # "Fix" pipeline
            mid_cpt = TestPipeline.MockComponentMid()
            pipe = Pipeline([
                start_cpt,
                mid_cpt,
            ], checkpoint=check_file)
            result = pipe.run()
            with check_file.open("rb") as check_in:
                checkpoint = pickle.load(check_in)

            self.assertDictContainsSubset(
                checkpoint.state, result,
                "Pipeline didn't alter previously-run data")

    def test_global_config(self):
        """Check we correctly pass config values to constructors."""
        dict_pipe = {
            "config": {"mid": {"foo": 1, "bar": 2}},
            "components": [
                {TestPipeline.MockComponentMid.qualname: [{"bar": 3}]}
            ]
        }
        pipe = Pipeline.load(dict_pipe)
        self.assertDictEqual(
            pipe.components[0].kwargs,
            {"foo": 1, "bar": 3},
            "Read arguments from config but allowed overriding")

    def test_statusfile(self):
        """Test status file after completing a pipeline."""

        with tempfile.NamedTemporaryFile("r") as status_file:
            pipeline = Pipeline([
                TestPipeline.MockComponentStart(),
            ], statusfile=status_file.name)

            status_components = [
                "INITIALISING",
                type(pipeline.components[0]).__name__,
                "FINISHED",
            ]

            pipeline.run()

            # Pipeline runs to completion, state is FINISHED
            status_file.seek(0)
            status_info = json.load(status_file)
            self.assertEqual(
                status_info,
                {"components": status_components, "index": 2})

            # Monkey patch component so that it raises. We can then check that
            # the index is correct.
            run_mock = MagicMock(side_effect=Exception(""))
            pipeline.components[0].run = run_mock

            # Disable warning about this overly-broad exception, because we
            # deliberately raised it.
            # pylint: disable=broad-except
            try:
                pipeline.run()
            except Exception:
                # deliberately raised
                pass

            # Pipeline fails on component 1, state is that component.
            status_file.seek(0)
            status_info = json.load(status_file)
            self.assertEqual(
                status_info,
                {"components": status_components, "index": 1})
