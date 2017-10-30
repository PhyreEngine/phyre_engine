"""
Components used for sorting fields in the pipeline state.

The components in this module can be used to sort a list within the pipeline
state. The list to be sorted should be specified to each component as the
`list` parameter. Components that accept a sorting key use that key relative
to the item being sorted

The field on which items are sorted is selected by supplying a `keys` field.
Keys are specified as a list of identifiers, allowing you to drill down into
nested items. The list of identifiers must be stored in the ``key`` field of a
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
...     {"keys": ["scores", 0]},
...     {"keys": ["scores", 1], "reverse": True},
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
import collections
import random

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
            self.keys = [{"keys": []}]
        else:
            self.keys = keys

        # Convert scalar to list so we can use the same logic in "run"
        if is_scalar(self.field):
            self.field = [self.field]

    def run(self, data, config=None, pipeline=None):
        """Sort pipeline state."""
        to_sort = get_nested(data, self.field)

        # Normalise the keys into a list of lists. Then sort, running from the
        # last to first, taking advantage of Python's stable sorting.
        for sort_key in reversed(self.keys):
            reverse = sort_key.get("reverse", False)

            # We can ignore this warning because we are calling the closure
            # immediately, not storing it for later.
            # pylint: disable=cell-var-from-loop
            if sort_key.get("allow_none", False):
                getter = lambda i: (i is None, get_nested(i, sort_key["keys"]))
            else:
                getter = lambda i: get_nested(i, sort_key["keys"])

            to_sort.sort(key=getter, reverse=reverse)
        set_nested(data, self.field, to_sort)
        return data

class Shuffle(Component):
    """
    Randomly shuffle a field of the pipeline state.

    :param field: Name of the field to shuffle. If this is a scalar value, it
        must point to a field at the top of the pipeline state. To drill down
        into nested objects, pass a list of identifiers.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, field, seed=None):
        self.field = field
        if is_scalar(self.field):
            self.field = [self.field]
        self.random = random.Random(seed)

    def run(self, data, config=None, pipeline=None):
        """Shuffle pipeline state."""
        to_shuffle = get_nested(data, self.field)
        self.random.shuffle(to_shuffle)
        return data


def get_nested(item, key):
    """
    Get element `key` from `item`. This handles the logic for drilling
    down into a list or dictionary according to a list of keys. Keys must
    be supplied as a list of identifiers.

    :param item: Container from which `key` will be extracted.
    :param key: String or list of identifiers.
    """
    field_value = item
    for k in key:
        field_value = field_value[k]
    return field_value

def set_nested(item, key, value):
    """Analogous to :py:meth:`.get`, but sets the key to `value`."""
    field_value = item
    for k in key[:-1]:
        field_value = field_value[k]
    field_value[key[-1]] = value

def is_scalar(item):
    """Return `True` if `item` is a string or a non-iterable type."""
    return isinstance(item, str) or not isinstance(item, collections.Iterable)
