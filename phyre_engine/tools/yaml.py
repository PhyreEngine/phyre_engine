"""
This module contains utility functions wrapping parts of the :py:mod:`yaml`
module. The big changes are:

* Automatically attempt to use the C versions of the YAML parser and emitter if
  they are available.

* Lists that do not contain a child list or dictionary will be converted into
  tuples, allowing tuples to be used as hash keys.

* The :py:func:`.load` and :py:func:`.dump` functions default to using the
  safe loader and dumper.

In addition to the :py:func:`.load` and :py:func:`.dump` functions, this module
also defines the YAML loader :py:class:`.ImplicitLoader`. This loader returns a
:py:class:`.ImplicitDocument` object rather than a `dict`, and understands the
``!template`` tag, which it parses into :py:class:`.TemplateString` objects.
The :py:class:`.ImplicitDocument` object returned by the loader contains
unresolved template strings (represented with :py:class:`.TemplateString`
classes) which may be resolved using the :py:meth:`.ImplicitDocument.resolve`
method. Templates are formatted using the :py:meth:`str.format` method.

This allows for fairly powerful templating. For example, consider the following
pipeline definition:

    .. code-block:: YAML

    pipeline:
      config:
        tool_a: /long/path/to/bin/dir/a
        tool_b: /long/path/to/bin/dir/b
        tool_c: /long/path/to/bin/dir/c

This is tedious, and relocating the tools (by, for example, switching conda
environments) requires three strings to be updated. In complex pipelines you
will inevitably miss one and break the pipeline in ways that aren't obvious
until runtime. Instead, we can take advantage of the implicit loader:

    .. code-block:: YAML

    pipeline:
      config:
        bin_dir: /long/path/to/bin/dir
        tool_a: "{bin_dir}/a"
        tool_b: "{bin_dir}/b"
        tool_c: "{bin_dir}/c"

The document may then be loaded and the templates resolved:

    >>> import io
    >>> import phyre_engine.tools.yaml as yaml
    >>> yaml_doc = io.StringIO('''
    ... pipeline:
    ...   config:
    ...     bin_dir: /long/path/to/bin/dir
    ...     tool_a: "{bin_dir}/a"
    ...     tool_b: "{bin_dir}/b"
    ...     tool_c: "{bin_dir}/c"''')
    >>> pipeline = yaml.load(yaml_doc, Loader=yaml.ImplicitLoader)
    >>> pipeline = pipeline.resolve(pipeline.unresolved["pipeline"]["config"])
    >>> pipeline["config"]
    {
        "bin_dir": "/long/path/to/bin/dir",
        "tool_a": "/long/path/to/bin/dir/a",
        "tool_b": "/long/path/to/bin/dir/b",
        "tool_c": "/long/path/to/bin/dir/c",
    }
"""
import collections.abc

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeLoader as SafeLoader, CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader, SafeDumper

class UnresolvedTemplateError(Exception):
    """Raised when an unresolved value is used in a format string."""
    def __init__(self, template):
        super().__init__("Unresolved template '{}' referenced.".format(
            str(template)))


class TemplateString:
    """
    Unresolved template string.

    :param str template_str: Template string (i.e. a string to be passed to
        :py:meth:`str.format`).
    """
    def __init__(self, template_str):
        self.template_str = template_str

    def resolve(self, data):
        """
        Resolve this template using the fields in `data`.

        :param data: Fields used to format the template string. We recommend
            that this is a `dict` and that you only use named fields in your
            template string, but you can pass a list and use numbered fields if
            you really want to.
        """
        if isinstance(data, collections.Sequence):
            return self.template_str.format(*data)
        elif isinstance(data, collections.Mapping):
            return self.template_str.format(**data)
        else:
            raise TypeError("Expected dict or list, got {}".format(type(data)))

    def __format__(self, _spec):
        """
        This method shouldn't be formatted, because it is supposed to represent
        an unresolved (i.e. unformatted) string.
        """
        raise UnresolvedTemplateError(self)

class ImplicitDocument:
    """
    Document that may contain unresolved template strings.

    This document should be resolved by calling :py:meth:`.resolve` before use.
    The unresolved document may be accessed via the attribute `unresolved`. It
    is valid (and encouraged) to pass all or part of `unresolved` to the
    `resolve` method, but you should take care to avoid circular references to
    unresolved templates.

    :param document: Results of :py:meth:`yaml.load` using a vanilla loader.
    """
    def __init__(self, document):
        self.unresolved = document
        # "None" indicates that we don't know to begin with.
        self.contains_unresolved = None

    def resolve(self, data, allow_unresolved=False):
        """
        Resolve the entire document. If unresolved elements remain in the
        document after processing, an :py:class:`ImplicitDocument` is returned.
        Otherwise, an object of type matching `document` is returned.

        :param data: Dict or list (we encourage you to use a `dict` and named
            fields in your templates) passed as arguments to
            :py:meth:`str.format`.

        :param element: String, list, or scalar to resolve.

        :param bool allow_unresolved: If `False`, a
            :py:exc:`.UnresolvedTemplateError` will be raised if an unresolved
            template is referenced in a format field. If `True`, the error will
            be caught and the :py:class:`.TemplateString` object will remain
            in the document, allowing for multi-pass resolution.

        :raises UnresolvedTemplateError: If an unresolved template is referenced
            and `allow_unresolved` is `False`.

        :raises KeyError: If an unknown field is referenced.
        """
        resolved = self._resolve_part(data, self.unresolved, allow_unresolved)
        if self.contains_unresolved:
            return type(self)(resolved)
        self.contains_unresolved = False
        return resolved

    def _resolve_part(self, data, element, allow_unresolved=False):
        """
        Resolve each unresolved template with `data` using
        :py:meth:`TemplateString.resolve`.
        """
        if (isinstance(element, collections.abc.Sequence) and
                not isinstance(element, (str, bytes))):
            return [self._resolve_part(data, i, allow_unresolved)
                    for i in element]
        elif isinstance(element, collections.abc.Mapping):
            return {k: self._resolve_part(data, v, allow_unresolved)
                    for k, v in element.items()}
        elif isinstance(element, TemplateString):
            try:
                return element.resolve(data)
            except UnresolvedTemplateError as err:
                if not allow_unresolved:
                    raise err
                self.contains_unresolved = True
                return element
        else:
            return element

class ImplicitLoader(yaml.SafeLoader):
    """
    Load a YAML document, parsing ``!template`` tags as
    :py:class:`.TemplateString` objects to be resolve later.

    Passing this class to :py:func:`.load` will cause that function to return a
    :py:class:`.ImplicitDocument` object.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_constructor("!template", self.template_str)

    def construct_document(self, node):
        return ImplicitDocument(super().construct_document(node))

    @staticmethod
    def template_str(self, node):
        return TemplateString(self.construct_scalar(node))


def _construct_yaml_tuple(self, node):
    # Used to convert sequences from lists to tuples. Only applies to lists
    # without any nested structures.
    seq = self.construct_sequence(node)
    if any(isinstance(e, (list, tuple, dict)) for e in seq):
        return seq
    return tuple(seq)


def dump(data, stream=None, Dumper=SafeDumper, **kwargs):
    """See :py:func:`yaml.dump`."""
    return yaml.dump(data, stream=stream, Dumper=Dumper, **kwargs)


def load(stream, Loader=SafeLoader):
    """See :py:func:`yaml.load`."""
    return yaml.load(stream, Loader)

SafeLoader.add_constructor('tag:yaml.org,2002:seq', _construct_yaml_tuple)
