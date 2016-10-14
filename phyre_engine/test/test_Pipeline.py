import unittest
from phyre_engine.component import Component
from phyre_engine import Pipeline

class TestPipeline(unittest.TestCase):
    """Check that we can assemble, validate and execute pipelines."""

    class MockComponentStart(Component):
        """Seed a couple of input elements."""
        REQUIRED = []
        ADDS = ["MCA_1", "MCA_2"]
        REMOVES = []

        def run(self, data):
            data["MCA_1"] = "MCA_1"
            data["MCA_2"] = "MCA_2"
            return data

    class MockComponentMid(Component):
        """Removes a required key, so will fail if called twice."""
        REQUIRED = ["MCA_1"]
        ADDS = ["MCB_1"]
        REMOVES = ["MCA_1"]

        def run(self, data):
            del data["MCA_1"]
            data["MCB_1"] = "MCB_1"
            return data

    class MockComponentBad(Component):
        """Badly-behaving mock component that does not add a key."""
        REQUIRED = []
        ADDS = ["MCA_1"]
        REMOVES = []

        def run(self, data):
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
