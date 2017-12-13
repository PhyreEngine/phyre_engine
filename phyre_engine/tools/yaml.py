"""
This module contains utility functions wrapping parts of the :py:mod:`yaml`
module. The big changes are:

* Automatically attempt to use the C versions of the YAML parser and emitter if
  they are available.

* Lists that do not contain a child list or dictionary will be converted into
  tuples, allowing tuples to be used as hash keys.

* The :py:func:`.load` and :py:func:`.dump` functions default to using the
  safe loader and dumper.
"""

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeLoader as SafeLoader, CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader, SafeDumper


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
