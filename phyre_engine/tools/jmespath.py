"""
Tools used by PhyreEngine when dealing with `JMESPath <http://jmespath.org/>`_
expressions.
"""
import jmespath


class JMESExtensions(jmespath.functions.Functions):
    """
    Class used to provide extensions to JMESPath.

    These can be called as functions. For example, you can retrieve the start
    and stop of a range object using the ``start`` and ``stop`` functions:

    .. code:: none

        {range_start: start(range_object), range_stop: stop(range_object)}


    ``root()``
        Returns the root of the pipeline state.

    ``toordinal(date)``
        Calls :py:meth:`datetime.date.toordinal` on the given
        :py:class:`datetime.date` object.

    ``start(range)`` and ``stop(range)``
        Return the start or stop of the supplied :py:class:`range` object.

    ``list(sequence)``
        Call the :py:func:`list` method on the given sequence. This is useful
        for converting non-list sequences such as tuples into lists. Without
        calling this function, the underlying JMESPath library will see a tuple
        as an opaque Python object.



    :param root: Root of the pipeline state.
    """

    def __init__(self, root):
        self.root = root

    @jmespath.functions.signature()
    def _func_root(self):
        return self.root

    @jmespath.functions.signature({"types": []})
    def _func_toordinal(self, date):
        return date.toordinal()

    @jmespath.functions.signature({"types": []})
    def _func_start(self, value):
        return value.start

    @jmespath.functions.signature({"types": []})
    def _func_stop(self, value):
        return value.stop

    @jmespath.functions.signature({"types": []})
    def _func_list(self, value):
        return list(value)

    @jmespath.functions.signature({"types": ["object"]}, {"types": ["array"]})
    def _func_delete(self, obj, fields):
        for field in fields:
            del obj[field]
        return obj
