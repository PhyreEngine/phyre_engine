"""
Components in this module mutate the pipeline state according to the results
of a `JMESPath <http://jmespath.org/>`_ expression. JMESPath is an XPath-like
language for selecting data from a data model that can be represented in JSON.

Two components are made available here, :py:class:`.Update` and
:py:class:`.Replace`. Both components take two expressions: one to select to
be modified and one that is evaluated in the context of the results of the
first expression and specifies how to alter each element.

We make several extensions available to each JMESPath expression. These are
described in the :py:class:`phyre_engine.tools.jmespath` class.

"""
from phyre_engine.component.component import Component
from phyre_engine.tools.jmespath import JMESExtensions
import collections.abc
import jmespath

class JMESPathBase(Component):
    """Base class for JMESPath components."""

    REQUIRED = []
    ADDS = []
    REMOVES = []

    CONFIG_SECTION = "jmespath"

    _TYPE_ERROR_MSG = (
        "Cannot replace '{select_expr}' (type: {select_expr_type}) "
        "with the result of '{value_expr}' (type: {value_expr_type}).")

    _INVALID_SELECTION_MSG = (
        "Invalid selection '{select_expr}' (type: {select_expr_type}). "
        "Expected either a dictionary or a list.")

    def __init__(self, select_expr, value_expr):
        self.select_expr = select_expr
        self.value_expr = value_expr

    def _invalid_selection(self, to_replace):
        raise TypeError(self._INVALID_SELECTION_MSG.format(
            select_expr=self.select_expr,
            select_expr_type=type(to_replace)))

    def _type_error(self, to_replace, replace_with):
        raise TypeError(self._TYPE_ERROR_MSG.format(
            select_expr=self.select_expr, select_expr_type=type(to_replace),
            value_expr=self.value_expr, value_expr_type=type(replace_with)))


class Update(JMESPathBase):
    """
    Alter an existing value in the pipeline state, retaining existing values.

    This component will operate on the dictionary or list of dictionaries
    returned by `select_expr`. If `select_expr` returns a dictionary, then the
    result of `select_expr` has its `update` method called with the results of
    `value_expr`. If `select_expr` returns a list, each element of the list is
    operated on in turn.

    For example, consider the following pipeline state:

    .. code-block:: python

        {"TM": 0.5, "n": 2}

    We could add an extra element to this by first selecting the root of the
    pipeline state with `select_expr` and then returning a dictionary from
    `value_expr`:

    >>> from phyre_engine.component.jmespath import Update
    >>> state = {"TM": 0.5, "n": 2}
    >>> Update("@", "{method: `Test`}").run(state)
    {"TM": 0.5, "n": 2, "method": "Test"}

    Here, the ``@`` symbol just means "current node", which is by default the
    root of the pipeline state. More usefully, this component could be used to
    copy a top-level parameter down a nested part of the pipeline state:

    >>> from phyre_engine.component.jmespath import Update
    >>> state = {
    ...     "structure": "native.pdb",
    ...     "templates": [
    ...         {"model": "model-1.pdb"},
    ...         {"model": "model-2.pdb"}
    ...     ]
    ... }
    >>> Update("templates", "{native: root().structure}").run(state)
    {
        "structure": "native.pdb",
        "templates": [
            {"native": "native.pdb", "model": "model-1.pdb"},
            {"native": "native.pdb", "model": "model-2.pdb"}
        ]
    }

    Here, we selected the ``templates`` list with the `select_expr`, and made
    use of the ``root()`` extension to examine the root of the pipeline state.
    This example could be useful to prepare the pipeline state for running a
    component from the :py:mod:`phyre_engine.component.strucaln` module.
    """

    def run(self, data, config=None, pipeline=None):
        """Update results of `self.select_expr` with `self.value_expr`."""
        jmespath_opts = jmespath.Options(
            custom_functions=JMESExtensions(data))

        to_replace = jmespath.search(self.select_expr, data, jmespath_opts)
        if isinstance(to_replace, collections.abc.Mapping):
            value = jmespath.search(self.value_expr, to_replace, jmespath_opts)
            if not isinstance(value, collections.abc.Mapping):
                self._type_error(to_replace, value)
            to_replace.update(value)
        elif isinstance(to_replace, collections.abc.Sequence):
            for item in to_replace:
                value = jmespath.search(self.value_expr, item, jmespath_opts)
                item.update(value)
        else:
            self._invalid_selection(to_replace)
        return data

class Replace(JMESPathBase):
    """
    Alter an existing value in the pipeline state, clearing existing values.

    This component operates in a similar manner to :py:class:`.Update`. The
    values returned by `select_expr` are *replaced* by the results of
    `value_expr`.

    If `select_expr` returns a list, the contents of that list are replaced
    with the results of `value_expr`, which should also be a list. If the
    result of `select_expr` is a dictionary, it is cleared and replaced with
    the result of `value_expr`.

    If the replacement value is not of the same type as the element to be
    replaced, then a :py:exc:`TypeError` will be raised. This is in part due to
    technical limitations of the JMESPath library, but in general it makes
    little sense to alter the pipeline state in such a way.

    >>> from copy import deepcopy
    >>> from phyre_engine.component.jmespath import Replace
    >>> state = {
    ...     "structure": "native.pdb",
    ...     "templates": [
    ...         {"model": "model-1.pdb"},
    ...         {"model": "model-2.pdb"}
    ...     ]
    ... }
    >>> Replace("templates", "@[].model").run(deepcopy(state))
    {
        "structure": "native.pdb",
        "templates": ["model-1.pdb", "model-2.pdb"]
    }
    >>> Replace("@", "{structure: templates[0].model}").run(deepcopy(state))
    {"structure": "model-1.pdb"}
    """

    def run(self, data, config=None, pipeline=None):
        """Replace results of `self.select_expr` with `self.value_expr`."""
        opts = jmespath.Options(custom_functions=JMESExtensions(data))

        to_replace = jmespath.search(self.select_expr, data, opts)
        replace_with = jmespath.search(self.value_expr, to_replace, opts)
        if isinstance(to_replace, collections.abc.Mapping):
            if not isinstance(replace_with, collections.abc.Mapping):
                self._type_error(to_replace, replace_with)
            # Edge case: If the two elements are the same, then we cannot
            # call "clear", because it will erase both. In that case, do
            # nothing.
            if to_replace is not replace_with:
                to_replace.clear()
                to_replace.update(replace_with)
        elif isinstance(to_replace, collections.abc.Sequence):
            if not isinstance(replace_with, collections.abc.Sequence):
                self._type_error(to_replace, replace_with)
            self.logger.info(
                ("Replacing result of '%s' (a list of length %d) "
                 "with a %d-element list"),
                self.select_expr, len(to_replace), len(replace_with))
            to_replace[:] = replace_with
        else:
            self._invalid_selection(to_replace)
        return data
