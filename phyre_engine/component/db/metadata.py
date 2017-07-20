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

from phyre_engine.component.component import Component
import re

class NameTemplate(Component):
    """
    Add a ``name`` attribute to each template.

    Each element in the ``templates`` array must contain a ``sequence`` key.
    Names are assigned to the ``name`` element of each template dictionary based
    on some transformation of the sequence name.

    Transformation functions may be supplied directly as a function or as a
    string, in which case the corresponding built-in function is called.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, name_fn="full_name", *args, **kwargs):
        """
        :param function name_fn: Function accepting a template dictionary and
            returning a name. Uses the full sequence name if not specified.
        """

        # Map for converting method names to actual functions.
        BUILTIN_SPLITTERS = {
            "full_name": self.full_name,
            "split": self.split,
        }

        if isinstance(name_fn, str):
            # Raise error if unsupported
            if name_fn not in BUILTIN_SPLITTERS:
                raise KeyError((
                    "Built-in splitter '{}' not found."
                    "Available splitters: {}").format(
                        name_fn, list(BUILTIN_SPLITTERS.keys())
                    ))

            # Save splitter function and args
            self.namer = (BUILTIN_SPLITTERS[name_fn], args, kwargs)
        elif callable(name_fn):
            self.namer = (name_fn, args, kwargs)
        else:
            raise TypeError("Argument name_fn must be a str or callable.")

    def run(self, data, config=None, pipeline=None):
        templates = self.get_vals(data)
        for template in templates:
            template["name"] = self.namer[0](
                template,
                *self.namer[1], **self.namer[2])
        return data

    @staticmethod
    def full_name(template):
        """Built-in splitter using the full name of the sequence."""
        return template["sequence"].name

    @staticmethod
    def split(template, delimiter, index):
        """Built-in splitter splitting on a delimiter."""
        return template["sequence"].name.split(delimiter)[index]

class DescribeTemplate(Component):
    """
    Alter the description of the template sequence.

    This component alters the ``.description`` field of the
    :py:class:`Bio.PDB.SeqRecord` object pointed to by the ``sequence`` key of
    each template.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, description_fn):
        """
        :param function name_fn: Function accepting a template dictinoary and
            returning a description.
        """
        self.description_fn = description_fn

    def run(self, data, config=None, pipeline=None):
        templates = self.get_vals(data)
        for template in templates:
            template["sequence"].description = self.description_fn(template)
        return data

class ParseSequenceName(Component):
    """
    Parse the sequence name (the ``name`` attribute of the ``sequence`` element)
    of each template in the ``templates`` list.

    A regular expression with *named captures* must be supplied to the
    constructor. All matching groups are assigned to the template with the
    corresponding keys.

    The following regular expressions may be useful:

    ``^(?P<name>.*)$``
        Use the entire sequence name to set the ``name`` attribute.

    ``^(?P<name>(?P<PDB>\w{4})_(?P<chain>\w+))``
        Match a PDB ID separated from a chain ID by an underscore, and store
        the PDB ID and chain ID in the ``PDB`` and ``chain`` elements. Set the
        ``name`` key to the PDB ID and chain ID, separated by an underscore.

    .. note::

        Regex matching is by default done using the ``ASCII`` flag, because it
        is rare to see sequences with non-ASCII characters in their metadata.
        Unicode matching may be enabled with the ``unicode_match`` parameter.

    :param str regex: Regular expression to search against the sequence name.
    :param bool must_match: If true, any non-matching templates cause an error.
    :param bool unicode_matching: Enable unicode matching, rather than the
        default ASCII-only matching.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = ["templates"]

    def __init__(self, regex, must_match=False, unicode_matching=False):
        self.must_match = must_match
        if not unicode_matching:
            self.flags = re.ASCII
        else:
            self.flags = 0
        self.regex = re.compile(regex, self.flags)

    def run(self, data, config=None, pipeline=None):
        """Parse sequence name using a regex."""

        templates = self.get_vals(data)
        for template in templates:
            match = self.regex.search(template["sequence"].name)
            if match is None:
                if self.must_match:
                    raise ValueError(
                        "Sequence name {} did not match regex {!s}".format(
                            template["sequence"].name, self.regex.pattern))
                continue
            template.update(match.groupdict())
        return data
