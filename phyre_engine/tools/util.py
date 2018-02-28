"""
Module containing various utility classes and functions that may be used in many
places.
"""
import collections.abc
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


class NamedTuple:
    """
    Implementation of a "named tuple", similar to
    :py:func:`collections.namedtuple`.

    A common pattern (at least in PhyreEngine) is to use a named tuple to
    represent the input or output of a class. For example, class ``ThingDoer``
    might represent its output as a ``ThingDoer.Result``:

    .. code-block:: python

        from collections import namedtuple
        class ThingDoer:
            Result = namedtuple("Result", "a b c")

            def do_thing(self):
                return self.Result(1, 2, 3)

    This way, all the little struct-like classes that tend to accumulate are
    kept local to the class that generates them, at least until it is necessary
    to lift them up to the module level. The named tuples are just a convenient
    way to implement this.

    Unfortunately, the named tuples cause problems with serialisation: some
    serialisers treat them like tuples, and some treat them like objects. More
    importantly, named tuples defined as class variables cannot be pickled. The
    following will raise a ``PicklingError``:

    .. code-block:: python

        import pickle
        result = ThingDoer().do_thing()
        pickle.dumps(result)

    This class is a reimplementation of named tuples that can be serialised
    effectively. To create a named tuple, subclass this and define a
    class-level ``FIELDS`` attribute:

    >>> from phyre_engine.tools.util import NamedTuple
    >>> class Result(NamedTuple):
    ...     FIELDS = "a b c"
    >>> Result(1, 2, 3)
    <Result a=1, b=2, c=3>
    >>> Result(a=4, b=5, c=6)
    <Result a=4, b=5, c=6>
    >>> x, y, z = Result(1, 2, 3)

    The ``FIELDS`` attribute may be defined either as a list of names or as a
    single space-separated string. Some magic (using ``__init_subclass__``)
    ensures that ``FIELDS`` is always a list when the class is examined.

    The class property ``_fields`` operates the same as it does for
    :py:func:`collections.namedtuple`, as does the instance method ``_asdict``.
    """


    _EXTRA_POS_ARGS_MSG = "Expected {} positional arguments but got {}."
    _FIELDS_REQUIRED_MSG = "Class {} does not have the class field 'FIELDS'."
    _EXTRA_ATTR_MSG = "'{}' object has no attribute '{}'."
    _MISSING_ATTRS_MSG = "Missing attributes '{}'"


    def __init__(self, *args, **kwargs):
        name = type(self).__name__
        fields = set(self.FIELDS)
        attrs = set()

        if len(args) > len(self.FIELDS):
            raise AttributeError(self._EXTRA_POS_ARGS_MSG.format(
                len(self.FIELDS), len(args)))

        for attr, value in zip(self.FIELDS, args):
            setattr(self, attr, value)
            attrs.add(attr)

        for attr, value in kwargs.items():
            if attr not in fields:
                raise AttributeError(self._EXTRA_ATTR_MSG.format(name, attr))
            setattr(self, attr, value)
            attrs.add(attr)

        missing_attrs = fields - attrs
        if missing_attrs:
            raise AttributeError(self._MISSING_ATTRS_MSG.format(missing_attrs))

    def __init_subclass__(cls):
        """
        Called whenever a subclass of this base class is generated.

        This method ensures that the ``FIELDS`` attribute actually exists, and
        converts it to a list of strings if a single space-separated string
        was supplied.
        """
        name = cls.__name__
        fields = getattr(cls, "FIELDS", None)
        if fields is None:
            raise NotImplementedError(cls._FIELDS_REQUIRED_MSG.format(name))

        if isinstance(fields, str):
            cls.FIELDS = fields.split()
        cls._fields = cls.FIELDS

    def _asdict(self):
        return {f: getattr(self, f) for f in self.FIELDS}

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [getattr(self, field) for field in self.FIELDS[index]]
        return getattr(self, self.FIELDS[index])

    def __repr__(self):
        return "<{} {}>".format(
            type(self).__name__,
            ", ".join([f + "=" + str(getattr(self, f)) for f in self.FIELDS]))

    def __iter__(self):
        return (getattr(self, f) for f in self.FIELDS)


def apply_dotted_key(dictionary, dotted_key, value):
    """
    Set a key deep in a dictionary. The key is split on each dot (``.``), and
    each level is assumed to be a nested map. For example, ``a.b.c`` will set
    the key ``{"a": {"b": {"c": value}}}``.
    """
    keys = dotted_key.split(".")
    dict_section = dictionary
    for key in keys[:-1]:
        if key not in dict_section:
            dict_section[key] = {}
        dict_section = dict_section[key]
    dict_section[keys[-1]] = value

def deep_merge(src, dst):
    """
    Merge two dictionaries, overwriting elements of `dst` with the
    corresponding elements of `src`. Dictionary elements are in turn merged.

    >>> from phyre_engine.tools.util import deep_merge
    >>> deep_merge({"a": {"x": 1}, "b": 2}, {"a": {"y": 2}, "b": 1})
    {"a": {"x": 1, "y": 2}, "b": 2}
    """
    for key, value in src.items():
        if isinstance(value, collections.abc.Mapping):
            node = dst.setdefault(key, {})
            deep_merge(value, node)
        else:
            dst[key] = value
    return dst
