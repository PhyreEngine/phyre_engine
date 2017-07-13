"""
Module containing utilities to make testing slightly easier.

PhyreEngine is designed to implement bioinformatics pipelines, and so it
contains several components that wrap external tools. These components require
some configuration, which cannot be taken for granted, and so some of the
PhyreEngine tests rely on the variable :py:data:`phyre_engine.test.config`,
which may contain a test configuration.

As well as defining :py:data:`phyre_engine.test.config`, this module contains
decorators aiming to make it easier to skip tests if the required configuration
parameters are not set.
"""

# The global variable "config" is not a constant, and so shouldn't be uppercase.
# Similarly, the decorators defined here are in camelCase because that is the
# (annoying) convention used by the unittest module.
# pylint: disable=invalid-name

import unittest
import functools

#: Can be set by the test runner
config = None

def requireConfig(func):
    """Decorator used to skip a test if the test configuration is not set."""
    if config is None:
        return unittest.skip("Test requires configuration to be set")(func)
    return func

def requireFields(fields, parents=None):
    """
    Decorator used to skip a test if the given fields are not set in the test
    configuration.

    For example, the following usage requires the fields ``foo`` and ``bar`` in
    the section ``section2``, which must be nested under ``section1``:

    >>> from phyre_engine.test import requireFields
    >>> @requireFields(["foo", "bar"], ["section1", "section2"])
    ... def test_widget():
    ...     ...

    This would correspond to the following configuration file (in YAML format):

    .. code-block:: YAML

        section1:
            section2:
                foo: true
                bar: "Hello world"

    :param list[str] fields: List of fields that must be defined.
    :param list[str] parents: Section of the configuration in which fields must
        be present.
    """
    if config is None:
        return unittest.skip("Test requires configuration to be set")

    parents = parents if parents is not None else []

    # pylint: disable=unsubscriptable-object, unsupported-membership-test
    config_section = config
    for i, parent in enumerate(parents):
        if parent not in config_section:
            return unittest.skip(
                "Section {sec} not found in {previous}".format(
                    sec=parent,
                    previous=parents[0:i]))
        config_section = config_section[parent]
    for field in fields:
        if field not in config_section:
            return unittest.skip("Field {field} not found".format(field=field))

    # pylint: disable=missing-docstring
    def decorator(obj):
        if not isinstance(obj, type):
            @functools.wraps(obj)
            def wrapper(*args, **kwargs):
                return obj(*args, **kwargs)
            return wrapper
        return obj
    return decorator
