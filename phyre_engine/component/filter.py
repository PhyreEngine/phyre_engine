"""Components for filtering keys from the pipeline state."""
from phyre_engine.component.component import Component
import jmespath

class Whitelist(Component):
    """
    Keep only the specified keys in the pipeline state.

    :param list[str] *args: List of keys to keep.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, *args):
        self.whitelist = set(args)

    def run(self, data, config=None, pipeline=None):
        """Apply a whitelist to the pipeline state."""
        new_state = {}
        for key, value in data.items():
            if key in self.whitelist:
                new_state[key] = value
        return new_state

class Blacklist(Component):
    """
    Remove the specified keys from the pipeline state.

    :param list[str] *args: List of keys to remove.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, *args):
        self.blacklist = set(args)

    def run(self, data, config=None, pipeline=None):
        """Apply a whitelist to the pipeline state."""
        new_state = {}
        for key, value in data.items():
            if key not in self.blacklist:
                new_state[key] = value
        return new_state


class Filter(Component):
    """
    Filter elements of a list by comparing fields in the pipeline state.

    This component uses `JMESPath <http://jmespath.org/>`_ expressions to select
    the list to be modified and to specify the filtered list.

    The extension function ``root(expr)`` is provided: it takes a JMESPath
    expression as its only argument and evaluates that expression relative to
    the root of the pipeline state.

    For example, to filter a list of templates in order to keep only those with
    a TM-score lower than the field ``TM`` in the root of the pipeline state:

    >>> from phyre_engine.component.filter import Filter
    >>> state = {"TM": 0.5, "templates": [{"TM": 0.3}, {"TM": 0.6}]}
    >>> filter = Filter("templates", "templates[?TM < root('TM')]")
    >>> filter.run(state)
    {"TM": 0.5, "templates": [{"TM": 0.3}]}

    This component performs an in-place list assignment. The contents of the
    list returned by the first JMESPath expression are replaced by the results
    of the list returned by the second expression. This allows the component to
    mutate the pipeline state in arbitrarily powerful ways, but it is
    recommended that you reserve this component for fairly simple uses. Create
    a component in Python if you find yourself manipulating the pipeline state
    in complex ways.

    :param str list_expr: JMESPath expression returning the list to be
        modified.

    :param str replace_expr: JMESPath expression that will be used to
        replace the result of `list_expr`.

    .. seealso::

        :py:class:`.JMESExtensions`
            For a full list of available extension functions.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    class JMESExtensions(jmespath.functions.Functions):
        """
        Class used to provide extensions to JMESPath.

        :param root: Root of the pipeline state.
        """

        def __init__(self, root):
            self.root = root

        @jmespath.functions.signature({"types": ["string"]})
        def _func_root(self, expr):
            """Evaluate the JMESPath `expr` relative to `self.root`."""
            return jmespath.search(expr, self.root)

        @jmespath.functions.signature({"types": []})
        def _func_toordinal(self, date):
            return date.toordinal()

    def __init__(self, list_expr, replace_expr):
        self.list_expr = list_expr
        self.replace_expr = replace_expr

    def run(self, data, config=None, pipeline=None):
        """Filter pipeline state."""
        jmespath_opts = jmespath.Options(
            custom_functions=self.JMESExtensions(data))
        to_replace = jmespath.search(self.list_expr, data, jmespath_opts)
        replace_with = jmespath.search(self.replace_expr, data, jmespath_opts)
        to_replace[:] = replace_with
        return data
