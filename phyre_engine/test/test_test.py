"""This module tests our custom test decorators."""
import unittest
import phyre_engine.test

# Missing docstrings are for our dummy functions, and invalid-name warnings are
# for methods overriden from TestCase.
# pylint: disable=missing-docstring, invalid-name

# Sample data used for testing
_FIELDS = ["foo", "bar"]
_SECTIONS = ["section1", "section2"]
_SAMPLE_CONFIG = {
    "section1": {
        "section2": {
            "foo": True,
            "bar": "baz"
        }
    }
}

class TestRequireConfig(unittest.TestCase):
    """
    The :py:func:`phyre_engine.test.requireConfig` decorator should skip tests
    if no config file was supplied.
    """

    def test_missing(self):
        """@requireConfig should cause tests to skip if no config is set."""

        # config is not set, so calling this should raise a SkipTest exception.
        @phyre_engine.test.requireConfig(override_config=None)
        def dummy_fail():
            return 1
        self.assertRaises(unittest.SkipTest, dummy_fail, "Test skipped")

    def test_present(self):
        """@requireConfig should not skip if config is set."""
        try:
            @phyre_engine.test.requireConfig(override_config=True)
            def dummy_pass():
                return 1
            self.assertEqual(dummy_pass(), 1, "Test not skipped")
        except unittest.SkipTest:
            self.fail("Test was skipped erroneously")

class TestRequireFields(unittest.TestCase):
    """
    The :py:func:`phyre_engine.test.requireFields` decorator should skip tests
    if the configuration is missing or the required fields are missing.
    """

    def _dummy(self):
        pass

    def test_missing_config(self):
        """@requireFields should skip if config is not set."""

        # This should fail because the configuration is None
        @phyre_engine.test.requireFields(_FIELDS, _SECTIONS,
                                         override_config=None)
        def dummy_fail_none():
            return 1
        self.assertRaises(
            unittest.SkipTest, dummy_fail_none,
            "Test skipped with config == None")

    def test_good_config(self):
        """Set test config to something useful and require those fields."""
        try:
            @phyre_engine.test.requireFields(_FIELDS, _SECTIONS,
                                             override_config=_SAMPLE_CONFIG)
            def dummy_pass():
                return 1
            self.assertEqual(
                dummy_pass(), 1,
                "Function call succeeds with correct config.")
        except unittest.SkipTest:
            self.fail("Test was skipped erroneously")

    def test_missing_field(self):
        """Fields not present in the config should cause a skip."""

        @phyre_engine.test.requireFields(_FIELDS + ["qux"], _SECTIONS)
        def dummy_fail_missing():
            return 1
        self.assertRaises(
            unittest.SkipTest, dummy_fail_missing,
            "Test skipped when missing fields")

    def test_missing_section(self):
        """Sections not present in the config should cause a skip."""

        @phyre_engine.test.requireFields(_FIELDS, _SECTIONS + ["qux"])
        def dummy_fail_missing():
            return 1
        self.assertRaises(
            unittest.SkipTest, dummy_fail_missing,
            "Test skipped when missing sections")

    def test_single_field(self):
        """Single fields should be treated like a list of length 1."""

        self.assertRaises(
            unittest.SkipTest,
            phyre_engine.test.requireFields(
                "bad", _SECTIONS, override_config=_SAMPLE_CONFIG
            )(self._dummy),
            "Skip when a single bad field is passed")
        try:
            phyre_engine.test.requireFields(
                "foo", _SECTIONS, override_config=_SAMPLE_CONFIG
            )(self._dummy)
        except unittest.SkipTest:
            self.fail("Should not have skipped.")

    def test_validator_dict(self):
        """We should be able to pass a dict of field names and validators."""

        # For brevity's sake
        deco = phyre_engine.test.requireFields
        bad_fields = {"foo": bool, "bar": lambda v: len(v) > 5}
        good_fields = {"foo": bool, "bar": lambda v: len(v) < 5}

        self.assertRaises(
            unittest.SkipTest,
            deco(
                bad_fields, _SECTIONS, override_config=_SAMPLE_CONFIG
            )(self._dummy),
            "Skip when a validator fails")

        try:
            phyre_engine.test.requireFields(
                good_fields, _SECTIONS, override_config=_SAMPLE_CONFIG
            )(self._dummy)
        except unittest.SkipTest:
            self.fail("Should not have skipped.")

class TestWithMissingConfig(unittest.TestCase):
    """
    This test case uses our test decorators in the wild, wrapped around test
    methods. This class sets the config to be ``None`` and ensure that
    decorators skip tests that rely on the config.
    """

    @phyre_engine.test.requireConfig(override_config=None)
    def test_require_config(self):
        """This should be skipped because config is None."""
        self.fail("This should have been skipped.")

    @phyre_engine.test.requireFields(["foo"], override_config=None)
    def test_require_fields(self):
        """This should be skipped because config is None."""
        self.fail("This should have been skipped.")

class TestWithConfig(unittest.TestCase):
    """
    Set a test configuration and make sure that our custom test decorators work
    with it when attached to actual test cases.
    """

    @phyre_engine.test.requireConfig(override_config=_SAMPLE_CONFIG)
    def test_require_config(self):
        """This should be run because test config was not None."""
        pass

    @phyre_engine.test.requireFields(_FIELDS, _SECTIONS,
                                     override_config=_SAMPLE_CONFIG)
    def test_require_fields_valid(self):
        """This should be run because test config is valid."""
        pass

    @phyre_engine.test.requireFields(_FIELDS + ["qux"], _SECTIONS,
                                     override_config=_SAMPLE_CONFIG)
    def test_require_fields_invalid(self):
        """This should not run because field "qux" is not in config."""
        self.fail("This should have been skipped.")

    @phyre_engine.test.requireFields(_FIELDS, _SECTIONS + ["qux"],
                                     override_config=_SAMPLE_CONFIG)
    def test_require_sections_invalid(self):
        """This should not run because section "qux" is not in config."""
        self.fail("This should have been skipped.")

class TestClassDecorators(unittest.TestCase):

    @phyre_engine.test.requireConfig(override_config=None)
    class TestRequireConfigClassSkipped(unittest.TestCase):
        """@requireConfig should skip this class because config is None."""

        def runTest(self):
            """We should have been skipped."""
            self.fail()

    @phyre_engine.test.requireFields(_FIELDS, _SECTIONS, override_config=None)
    class TestRequireFieldsClassSkipped(unittest.TestCase):
        """@requireFields should skip this class because config is None."""

        def runTest(self):
            """We should have been skipped."""
            self.fail()

    @phyre_engine.test.requireConfig(override_config={})
    class TestRequireConfigClassNotSkipped(unittest.TestCase):
        """This class should NOT be skipped because the config is not None."""

        def runTest(self):
            """Class not skipped."""
            pass

    @phyre_engine.test.requireFields(_FIELDS, _SECTIONS,
                                     override_config=_SAMPLE_CONFIG)
    class TestRequireFieldsClassNotSkipped(unittest.TestCase):
        """@requireFields should NOT skip this class."""

        def runTest(self):
            """Class not skipped."""
            pass

    @phyre_engine.test.requireFields(_FIELDS + ["qux"], _SECTIONS,
                                     override_config=_SAMPLE_CONFIG)
    class TestRequireFieldsClassSkippedWithInvalidFields(unittest.TestCase):
        """@requireFields should skip this class because of invalid fields."""

        def runTest(self):
            """We should have been skipped."""
            self.fail()

    @phyre_engine.test.requireFields(_FIELDS, _SECTIONS + ["qux"],
                                     override_config=_SAMPLE_CONFIG)
    class TestRequireFieldsClassSkippedWithInvalidSections(unittest.TestCase):
        """@requireFields should skip this class because of invalid sections."""

        def runTest(self):
            """We should have been skipped."""
            self.fail()

    # Mapping of test classes to a boolean indicating whether they should be
    # skipped (True) or not (False).
    _TESTS = {
        TestRequireConfigClassSkipped: True,
        TestRequireFieldsClassSkipped: True,
        TestRequireConfigClassNotSkipped: False,
        TestRequireFieldsClassNotSkipped: False,
        TestRequireFieldsClassSkippedWithInvalidFields: True,
        TestRequireFieldsClassSkippedWithInvalidSections: True,
    }

    def test_class_decorators(self):
        """Test class decorators."""
        result = unittest.TestResult()
        suite = unittest.TestSuite()
        for test in self._TESTS:
            suite.addTest(test())
        suite.run(result)

        self.assertEqual(result.testsRun, len(self._TESTS))
        for test, reason in result.skipped:
            with self.subTest("Skipped test", test=test, reason=reason):
                self.assertTrue(self._TESTS[type(test)])
