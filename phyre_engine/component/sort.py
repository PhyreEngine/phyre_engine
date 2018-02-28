"""
Components used for sorting fields in the pipeline state.

The components in this module can be used to sort a list within the pipeline
state. The list to be sorted should be specified to each component as the
`list` parameter. Components that accept a sorting key use that key relative
to the item being sorted

The field on which items are sorted is selected by supplying a `key` field.
Keys are specified as JMESPath paths, allowing you to drill down into
nested items. The key path must be stored in the ``key`` field of a
dictionary. Each dictionary also accepts a ``reverse`` field and a
``allow_none`` field, which respectively reverse the sort order and allow `None`
values when sorting. The primary sort happens according to the first key, then
the second, and so on:

>>> from phyre_engine.component.sort import Sort
>>> pipe_state = {
...     "templates": [
...         {"scores": (0.3, 0.1)},
...         {"scores": (0.1, 0.1)},
...         {"scores": (0.4, 0.2)},
...         {"scores": (0.1, 0.3)},
...         {"scores": (0.5, 0.5)},
...     ]
... }
>>> sort_component = Sort(field="templates", keys=[
...     {"key": "scores[0]"},
...     {"key": "scores[1]", "reverse": True},
... ])
>>> sort_component.run(pipe_state)
>>> {
...     "templates": [
...         {"scores": (0.1, 0.3)}, # Ties broken by ["scores"][1]
...         {"scores": (0.1, 0.1)},
...         {"scores": (0.3, 0.1)},
...         {"scores": (0.4, 0.2)},
...         {"scores": (0.5, 0.5)},
...     ]
... }

If the `keys` parameter is not defined, the list will be sorted into ascending
order using the default comparison operators.
"""
from phyre_engine.component.component import Component
from phyre_engine.tools.jmespath import JMESExtensions
import collections
import random
import jmespath

class Sort(Component):
    """
    Sort the pipeline. See :py:mod:`phyre_engine.component.sort` for more
    details.

    :param field: Key of the list to sort. Can either be specified as a single
        string to sort a list at the top of the pipeline state, or as a list of
        identifiers to drill down.

    :param key: List of keys to sort by. Each key must be a list, the first
        element of which must be "ascending" or "descending".
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, field, keys=None):
        self.field = field

        if keys is None:
            self.keys = [{"key": "@"}]
        else:
            self.keys = keys

    def run(self, data, config=None, pipeline=None):
        """Sort pipeline state."""
        jmespath_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        to_sort = jmespath.search(self.field, data, jmespath_opts)

        # Sort according to each key, running from last to first to take
        # advantage of Python's stable sorting.
        for sort_key in reversed(self.keys):
            reverse = sort_key.get("reverse", False)
            allow_none = sort_key.get("allow_none", False)
            to_sort = jmes_sort(
                to_sort, sort_key["key"], root=data,
                reverse=reverse, allow_none=allow_none)
        jmespath.search(self.field, data, jmespath_opts)[:] = to_sort
        return data

class Shuffle(Component):
    """
    Randomly shuffle a field of the pipeline state.

    :param str field: JMESPath expression giving the field to shuffle.
    :param int seed: Random seed to use for shuffling.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, field, seed=None):
        self.field = field
        self.random = random.Random(seed)

    def run(self, data, config=None, pipeline=None):
        """Shuffle pipeline state."""
        jmespath_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        to_shuffle = jmespath.search(self.field, data, jmespath_opts)
        self.random.shuffle(to_shuffle)
        return data

def jmes_sort(data, jmespath_key, root=None, reverse=False, allow_none=False):
    """
    Sort `data` according to JMESPath key.

    Returns a list sorted according to the result of the JMESPath expression
    `jmespath_key`, which is evaluated relative to the item being sorted.

    For example, to sort a list by the value of each item, use the JMESPath
    expression "``@``", meaning "the entire item".

    >>> from phyre_engine.component.sort import jmes_sort
    >>> data = [3, 1, 4, 1, 5, 9]
    >>> jmes_sort(data, "@")
    [1, 1, 3, 4, 5, 9]

    More complicated paths can be used to search according to a key buried deep
    in the hierarchy:

    >>> from phyre_engine.component.sort import jmes_sort
    >>> data = [{"x": {"y": 1}}, {"x": {"y": -1}}, {"x": {"y": None}}]
    >>> jmes_sort(sort(data, "x.y", allow_none=True)
    [{'x': {'y': -1}}, {'x': {'y': 1}}, {'x': {'y': None}}]

    :param list to_sort: List to sort.
    :param str jmespath_key: JMESPath expression relative to each item,
        returning the value by which to sort the item.
    :param dict root: Pipeline state, for use with the ``root()`` extension.
    :param bool reverse: Reverse the sort order.
    :param bool allow_none: If `True`, allow `None` values in the list to be
        sorted. `None` is sorted to the end of the list (or the beginning) if
        `reverse=True`.
    :return: Sorted list.
    :rtype: list
    """
    def key_fn(datum):
        """Getter closure using `jmespath_key`."""
        jmespath_opts = jmespath.Options(custom_functions=JMESExtensions(root))
        field_value = jmespath.search(jmespath_key, datum, jmespath_opts)
        if allow_none:
            return (field_value is None, field_value)
        return field_value

    return sorted(data, key=key_fn, reverse=reverse)
