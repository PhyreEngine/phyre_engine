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

# Use the default config
DEFAULT_CONFIG = object()

def requireConfig(override_config=DEFAULT_CONFIG):
    """Decorator used to skip a test if the test configuration is not set."""
    if override_config is DEFAULT_CONFIG:
        override_config = config

    def decorator(obj):
        if override_config is None:
            return unittest.skip("Test requires configuration to be set")(obj)
        return obj
    return decorator

def requireFields(fields, parents=None, override_config=DEFAULT_CONFIG):
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

    Fields may be specified in any of the following ways:

    Single field
        A single string ``foo`` is interpreted identically to ``[foo]``; that
        is, it must be present and be ``True`` when coerced to a ``bool``.

    List of fields
        A list of fields ``[x, y, z]`` is interpreted identically to the dict
        ``{x: bool, y: bool, z: bool}``. That is, each field must be present
        and ``True`` when coerced to a ``bool``.

    Dictionary of fields and validators
        In the most general case, a dictionary of field names and validation
        functions (or any callable) may be passed. Each field must be present
        and the result of the validation function must be ``True``.

    :param dict[str, callable] fields: List of fields that must be defined.
    :param list[str] parents: Section of the configuration in which fields must
        be present.
    """
    if override_config is DEFAULT_CONFIG:
        override_config = config

    if override_config is None:
        return unittest.skip("Test requires configuration to be set")

    # Convert single fields to a dict
    if isinstance(fields, str):
        fields = {fields: bool}
    elif isinstance(fields, list):
        fields = {f: bool for f in fields}

    parents = parents if parents is not None else []

    # pylint: disable=unsubscriptable-object, unsupported-membership-test
    config_section = override_config
    for i, parent in enumerate(parents):
        if parent not in config_section:
            return unittest.skip(
                "Section {sec} not found in {previous}".format(
                    sec=parent,
                    previous=parents[0:i]))
        config_section = config_section[parent]
    for field, validator in fields.items():
        if field not in config_section:
            return unittest.skip("Field {field} not found".format(field=field))
        if not validator(config_section[field]):
            return unittest.skip("Validation of {field} failed".format(
                field=field))

    # pylint: disable=missing-docstring
    def decorator(obj):
        if not isinstance(obj, type):
            @functools.wraps(obj)
            def wrapper(*args, **kwargs):
                return obj(*args, **kwargs)
            return wrapper
        return obj
    return decorator
