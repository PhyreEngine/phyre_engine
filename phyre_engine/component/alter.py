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
