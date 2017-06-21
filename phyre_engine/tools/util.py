"""
Module containing various utility classes and functions that may be used in many
places.
"""
import os

class TemporaryEnvironment:
    """
    Context manager used for temporarily setting an environment variable.

    >>> with AlteredEnvironment(PATH=/tmp/path):
    >>>     # Do work with altered environment
    >>> # Environment is back to normal
    """

    def __init__(self, **kwargs):
        """
        :param kwargs: Environment variables to set.
        """
        self.env = kwargs
        self._original_env = None

    def __enter__(self):
        self._original_env = os.environ.copy()
        os.environ.update(self.env)

    def __exit__(self, *args):
        os.environ.clear()
        os.environ.update(self._original_env)
