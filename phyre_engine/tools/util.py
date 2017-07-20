"""
Module containing various utility classes and functions that may be used in many
places.
"""
import os
import pathlib

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

class Stream:
    """
    Context manager used to open streams starting from a filename,
    :py:class:`pathlib.Path` object, or existing stream. In the case of
    pre-existing streams, the stream is *not* closed when the context manager is
    exited.

    :param file_like: Either a filename, stream or :py:class:`pathlib.Path`.
    :param str mode: Stream mode. This is mandatory, but is ignored if the
        stream is already open.

    >>> from phyre_engine.util import Stream
    >>> with Stream("foo.txt", "w") as stream:
    ...     print("Hello world", file=stream)

    >>> with open("foo.txt", "w") as file_handle:
    ...    with Stream(file_handle, "w") as stream:
    ...        print("Hello again", file=stream)

    >>> from pathlib import Path
    >>> with Stream(Path("foo.txt"), "w") as stream:
    ...    print("Once more", file=stream)
    """

    def __init__(self, file_like, mode):
        self.file_like = file_like
        self.mode = mode
        # Set when a stream needs to be opened. Closed on exit.
        self._stream = None

    def __enter__(self):
        if isinstance(self.file_like, str):
            self._stream = open(self.file_like, self.mode)
            return self._stream
        elif isinstance(self.file_like, pathlib.Path):
            self._stream = self.file_like.open(self.mode)
            return self._stream
        return self.file_like

    def __exit__(self, *_args):
        if self._stream is not None:
            self._stream.close()
