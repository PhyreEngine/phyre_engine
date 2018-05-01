"""
Components used for setting sequence metadata when building a fold database.
Here, "metadata" is used to mean the sequence name and description. For example,
the following FASTA sequence has a ``name`` of ``1ABC_D`` and a ``description``
of ``Example sequence``:

.. code-block:: none

    > 1ABC_D Example sequence
    ...

The components in this module operate on the ``templates`` list, which is stored
in the pipeline state. Each element of the ``templates`` list is a dictionary,
and components in this module add or modify keys in this dictionary based on
other information about the template. For example, the :py:class:`.NameTemplate`
component parses information from the ``sequence`` key of each template and
sets the ``name`` key.
"""

from phyre_engine.component import Component
import re

class RegexComponent(Component):
    """
    Base class of regex-using metadata parsers. Components that read a regular
    expression from a string should inherit from this to ensure that regular
    expressions are treated consistently across components.

    :param str regex: Regular expression to search against the sequence name.
    :param bool must_match: If true, any non-matching templates cause an error.
    :param bool unicode_matching: Enable unicode matching, rather than the
        default ASCII-only matching.

    .. note::

        Regex matching is by default done using the ``ASCII`` flag, because it
        is rare to see sequences with non-ASCII characters in their metadata.
        Unicode matching may be enabled with the ``unicode_match`` parameter.
    """

    def __init__(self, regex, must_match=False, unicode_matching=False):
        self.must_match = must_match
        if not unicode_matching:
            self.flags = re.ASCII
        else:
            self.flags = 0
        self.regex = re.compile(regex, self.flags)

    def _search(self, haystack):
        """
        Search for ``self.regex`` in ``haystack``. Aborts if ``must_match``
        and no match was found.
        """
        match = self.regex.search(haystack)
        if match is None:
            if self.must_match:
                raise ValueError(
                    "Field '{}' did not match regex {!s}".format(
                        haystack, self.regex.pattern))
        return match

class ParseField(RegexComponent):
    """
    Apply a regular expression to a field of each dictionary in the
    ``templates`` list, adding fields matching all named captures. This may be
    useful, for example, for parsing a PDB ID and chain from the ``name``
    element of a template parsed by
    :py:class:`phyre_engine.component.hhsuite.ReportParser`.

    In addition to the ``field`` parameter, this component accepts the
    parameters of its base class, :py:class:`.RegexComponent`.

    The following regular expressions may be useful when applied to the ``name``
    field of a sequence parsed from a file:

    ``^(?P<name>.*)$``
        Use the entire sequence name to set the ``name`` attribute.

    ``^(?P<name>(?P<PDB>\w{4})_(?P<chain>\w+))``
        Match a PDB ID separated from a chain ID by an underscore, and store
        the PDB ID and chain ID in the ``PDB`` and ``chain`` elements. Set the
        ``name`` key to the PDB ID and chain ID, separated by an underscore.

    :param str field: Name of the field to parse.
    """
    CONFIG_SECTION = "metadata"

    REMOVES = []

    @property
    def REQUIRED(self):
        return [self.field]

    @property
    def ADDS(self):
        return list(self.regex.groupindex)

    def __init__(self, field, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.field = field

    def run(self, data, config=None, pipeline=None):
        """Parse a field using a regex."""
        match = self._search(data[self.field])
        if match is not None:
            data.update(match.groupdict())
        return data
