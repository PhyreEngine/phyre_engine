"""Set operations that can be applied to the pipeline state."""
import collections.abc
import jmespath

from phyre_engine.component import Component
from phyre_engine.tools.jmespath import JMESExtensions

class SetOperation(Component):
    """
    Base class for components that apply a set operation.

    Sets are chosen according to `sets`, a list of JMESPath expressions. The
    value used to compare sets is chosen by evaluating the JMESPath expression
    `key` relative to each item in the set. If an item appears in multiple
    sets, the version of that item that appears in the leftmost set in the
    `sets` list is preferred. If identical identifiers are found in a single
    set, the first identifier is the one that is kept, preserving the sort
    order.

    The results of the operation are stored in the `destination` field.

    :param list[str] sets: List of JMESPath expressions, each of which must
        return a list. The lists are the sets on which the set operation is
        applied.

    :param str key: JMESPath expression evaluated relative to each element of
        each set. This is the key used to determine when two items are
        equivalent. If the JMESPath expression returns a list, it is coerced
        into a tuple.

    :param str destination: Name of the field in which the results of the
        operation are stored.
    """
    @property
    def ADDS(self):
        return [self.destination]

    REMOVES = []
    REQUIRED = []

    def __init__(self, sets, key, destination):
        self.jmespath_sets = sets
        self.jmespath_key = key
        self.destination = destination

    def _sets(self, root):
        """
        Return a list of sets, each of which contains the results of
        evaluating `self.jmespath_key` on each item. This also returns a map
        of those identifiers to the corresponding objects.
        """

        key_map = {}
        key_sets = []
        opts = jmespath.Options(custom_functions=JMESExtensions(root))
        # Reverse the list of sets so the assignment to key_map occurs in
        # reverse order. That is, we prefer to keep elements from the first
        # set over the last set.
        for set_expr in reversed(self.jmespath_sets):
            item_set = set()
            item_list = jmespath.search(set_expr, root, opts)

            # Get list of identifiers by evaluating self.jmespath_key on each
            # item. We explicitly turn any sequences into tuples so they can
            # be hashed.
            identifiers = []
            for item in item_list:
                identifier = jmespath.search(self.jmespath_key, item, opts)
                if isinstance(identifier, collections.abc.Sequence):
                    identifier = tuple(identifier)
                identifiers.append(identifier)

            # Build map of identifiers to items, and list of sets of IDs.
            for ident, item in zip(identifiers, item_list):
                if ident not in item_set:
                    key_map[ident] = item
                    item_set.add(ident)
            key_sets.append(item_set)
        # Compensate for initial "reverse"
        key_sets.reverse()
        return key_sets, key_map

class Difference(SetOperation):
    """
    Calculate difference between sets of items.

    The difference is evaluated from left to right. That is, if `sets` is
    ``["a", "b", "c"]``, then the result is :math:`a - b - c`.

    .. seealso::

        :py:class:`.SetOperation`
            For the constructor parameters.
    """

    def run(self, data, config=None, pipeline=None):
        """Calculate set difference."""
        key_sets, key_map = self._sets(data)
        result = key_sets[0]
        for operand in key_sets[1:]:
            result = result - operand
        data[self.destination] = [key_map[k] for k in result]
        return data


class Union(SetOperation):
    """
    Take the union of each set of items.

    .. seealso::

        :py:class:`.SetOperation`
            For the constructor parameters.
    """

    def run(self, data, config=None, pipeline=None):
        """Calculate set union."""
        key_sets, key_map = self._sets(data)
        result = set().union(*key_sets)
        data[self.destination] = [key_map[k] for k in result]
        return data


class Intersection(SetOperation):
    """
    Take the intersection of each set of items.

    .. seealso::

        :py:class:`.SetOperation`
            For the constructor parameters.
    """

    def run(self, data, config=None, pipeline=None):
        """Calculate set union."""
        key_sets, key_map = self._sets(data)
        result = key_sets[0].intersection(*key_sets[1:])
        data[self.destination] = [key_map[k] for k in result]
        return data
