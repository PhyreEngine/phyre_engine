"""
Tools used by PhyreEngine when dealing with `JMESPath <http://jmespath.org/>`_
expressions.
"""
import datetime
import re

import jmespath


class JMESExtensions(jmespath.functions.Functions):
    """
    Class used to provide extensions to JMESPath.

    These can be called as functions. For example, you can retrieve the start
    and stop of a range object using the ``start`` and ``stop`` functions:

    .. code-block:: none

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

    ``delete(object, fields)``
        Delete each field in the ``fields`` list from ``object``.

    ``date(string, format)``
        Parse ``string`` into a :py:class:`datetime.date` object. The
        ``format`` string is interpreted as for
        :py:meth:`~datetime.datetime.strptime`.

    ``datetime(string, format)``
        Similar to ``date()``, but returning a :py:class:`~datetime.datetime`
        object.

    ``except(object, fields)``
        Return a copy of ``object`` with all of ``fields`` removed.

    ``select(object, fields)``
        The opposite of ``except``: keep only the keys in ``fields``.

    ``regex_keys(object, fields)``
        Return a list of the keys in ``object`` that match any of the regexes
        in ``fields``.

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
            if field in obj:
                del obj[field]
        return obj

    @jmespath.functions.signature({"types": ["string"]}, {"types": ["string"]})
    def _func_date(self, string, format_str):
        return datetime.datetime.strptime(string, format_str).date()

    @jmespath.functions.signature({"types": ["string"]}, {"types": ["string"]})
    def _func_datetime(self, string, format_str):
        return datetime.datetime.strptime(string, format_str)

    @jmespath.functions.signature({"types": ["object"]}, {"types": ["array"]})
    def _func_except(self, obj, fields):
        copy = obj.copy()
        for field in fields:
            if field in copy:
                del copy[field]
        return copy

    @jmespath.functions.signature({"types": ["object"]}, {"types": ["array"]})
    def _func_regex_keys(self, obj, fields):
        regexes = [re.compile(regex) for regex in fields]
        matching_fields = set()

        for regex in regexes:
            for key in obj.keys():
                if regex.search(key):
                    matching_fields.add(key)
        return list(matching_fields)

    @jmespath.functions.signature({"types": ["object"]}, {"types": ["array"]})
    def _func_select(self, obj, keys):
        result = {}

        for key in keys:
            if key in obj:
                result[key] = obj[key]
        return result
