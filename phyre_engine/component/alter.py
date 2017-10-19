"""
This module contains components for altering the pipeline state without regard
to the content by removing, moving or copying fields.
"""
import copy
from phyre_engine.component import Component

class _AlterationComponent(Component):
    """Base class for components in this module."""
    ADDS = []
    REMOVES = []
    REQUIRED = []

class Remove(_AlterationComponent):
    """
    Remove a field from the pipeline state.

    :param str field: Name of the field to be removed.
    """

    def __init__(self, field):
        self.field = field

    def run(self, data, config=None, pipeline=None):
        """Remove a field from the pipeline state."""
        del data[self.field]
        return data

class Copy(_AlterationComponent):
    """
    Make a copy of a field. A shallow copy is performed by default.

    :param str src: Name of the source field.
    :param str dst: Name of the destination field.
    :param bool deep: Perform a deep copy.
    """

    def __init__(self, src, dst, deep=False):
        self.src = src
        self.dst = dst
        self.copy_fn = copy.deepcopy if deep else copy.copy

    def run(self, data, config=None, pipeline=None):
        """Remove a field from the pipeline state."""
        data[self.dst] = self.copy_fn(data[self.src])
        return data

class Move(_AlterationComponent):
    """
    Rename a field.

    :param str from_field: Name of the field to be renamed.
    :param str to_field: Name of the new field.
    """

    def __init__(self, from_field, to_field):
        self.from_field = from_field
        self.to_field = to_field

    def run(self, data, config=None, pipeline=None):
        """Rename a field."""
        data[self.to_field] = data[self.from_field]
        del data[self.from_field]
        return data

class Set(_AlterationComponent):
    """
    Set a field in the pipeline state.

    :param str field: Name of the field to set.
    :param str value: Value to set
    """

    def __init__(self, field, value):
        self.field = field
        self.value = value

    def run(self, data, config=None, pipeline=None):
        """Set field in pipeline state."""
        data[self.field] = self.value
        return data
