"""
Components used for sorting fields in the pipeline state.

The components in this module can be used to sort a list within the pipeline
state. The list to be sorted should be specified to each component as the
`list` parameter. Components that accept a sorting key use that key relative
to the item being sorted

The field on which items are sorted is selected by supplying a `keys` field.
Keys are specified as a list of identifiers, allowing you to drill down into
nested items. The first element of each key must be either ``ascending`` or
``descending``, indicating the sort order. The primary sort happens according to
the first key, then the second, and so on:

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
...     ["ascending", "scores", 0],
...     ["descending", "scores", 1],
    ])
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

The default list of `keys` is `[["ascending"]]`, which will sort the list
chosen by the `fields` parameter into ascending order.
"""
from phyre_engine.component.component import Component
import collections

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

    SORT_DIRECTIONS = ("ascending", "descending")

    def __init__(self, field, keys=(("ascending",),)):
        self.field = field
        self.keys = keys

        self.validate_keys()
        # Convert scalar to list so we can use the same logic in "run"
        if is_scalar(self.field):
            self.field = [self.field]

    def run(self, data, config=None, pipeline=None):
        """Sort pipeline state."""
        to_sort = get_nested(data, self.field)

        # Normalise the keys into a list of lists. Then sort, running from the
        # last to first, taking advantage of Python's 
        for sort_key in reversed(self.keys):
            reverse = True if sort_key[0] == "descending" else False
            # We can ignore this warning because we are calling the closure
            # immediately, not storing it for later.
            # pylint: disable=cell-var-from-loop
            to_sort.sort(key=lambda i: get_nested(i, sort_key[1:]),
                         reverse=reverse)
        set_nested(data, self.field, to_sort)
        return data

    def validate_keys(self):
        """
        Check that the list of keys are a list of lists, the first element of
        which is either ``ascending`` or ``descending``.
        """
        if not self.keys:
            raise ValueError("List of keys cannot be empty.")

        for i, key in enumerate(self.keys):
            if key[0] not in self.SORT_DIRECTIONS:
                raise ValueError(
                    "Invalid sort direction {} for key {}".format(key[0], i))

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
