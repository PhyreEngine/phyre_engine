import unittest
import phyre_engine.component.update.build as build
import tempfile
import yaml

class TestBuild(unittest.TestCase):

    def test_nested_format(self):
        """Test recursive string formatting of nested object."""

        input_field = {
            "foo": {
                "dirs":[
                    "a/{test}",
                    "b/{test}",
                    {"c": "{test}"}
                ],
                "static": 1337,
                "deep": {"dictionary": "{test}"}
            }
        }
        expected_output = {
            "foo": {
                "dirs":[
                    "a/replacement",
                    "b/replacement",
                    {"c": "replacement"}
                ],
                "static": 1337,
                "deep": {"dictionary": "replacement"}
            }
        }
        pipe_state = {
            "update_required": [{"name": "foo", "test": "replacement"}]
        }

        with tempfile.NamedTemporaryFile() as tmpfile:
            conf_updater = build.UpdateConfigFile(tmpfile.name, input_field)
            conf_updater.run(pipe_state)
            tmpfile.seek(0)
            serialised = yaml.load(tmpfile, yaml.SafeLoader)

        self.assertDictEqual(serialised, expected_output, "Nested replacement")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
