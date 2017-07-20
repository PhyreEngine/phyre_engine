"""
This module contains components for dumping the pipeline state in various
formats. By default, each dumper will write to standard output. The ``output``
parameter may be set to either a file name or output stream to determine where
the output will be written.
"""
import json.encoder
import sys
from phyre_engine.component.component import Component
from phyre_engine.tools.util import Stream
import Bio.SeqRecord

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeDumper

class JsonStateEncoder(json.JSONEncoder):
    """
    Custom JSON encoder for pipeline state.

    This decoder includes JSON serialisers for the following objects types:

    :py:class:`Bio.SeqRecord` objects:
        Sequences objects are serialised as a FASTA string.

    :py:class:`range` objects:
        Ranges are serialised as tuples containing a start and stop. Ranges with
        a step other than 1 will result on a :py:exc:`TypeError` being raised.
    """

    def default(self, o):
        """Called to encode a."""
        # Disable warnings about hidden method; this is the recommended usage.
        # pylint: disable=method-hidden
        if isinstance(o, Bio.SeqRecord.SeqRecord):
            return serialise_seqrecord(o)
        elif isinstance(o, range) and o.step == 1:
            return serialise_range(o)
        else:
            super().default(o)

class YamlStateDumper(SafeDumper):
    """
    Custom decoder for pipeline state. See :py:class:`.JsonStateEncoder` for the
    types encoded.
    """
    # Disable warning about too many ancestors. Don't blame me, blame pyyaml.
    # pylint: disable=too-many-ancestors

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_representer(Bio.SeqRecord.SeqRecord, self._serialise_seqrecord)
        self.add_representer(range, self._serialise_range)

    @staticmethod
    def _serialise_seqrecord(dumper, record):
        return dumper.represent_scalar(
            "tag:yaml.org,2002:str",
            serialise_seqrecord(record),
            style="|")

    @staticmethod
    def _serialise_range(dumper, range_obj):
        if range_obj.step != 1:
            raise yaml.representer.RepresenterError(
                "Cannot represent range {!s} with step {}".format(
                    range_obj, range_obj.step))

        return dumper.represent_sequence(
            "tag:yaml.org,2002:seq",
            serialise_range(range_obj))


def serialise_seqrecord(record):
    """Use :py:mod:`Bio.SeqIO` to generate a FASTA sequence."""
    return record.format("fasta")

def serialise_range(range_obj):
    """Return a tuple of ``(stop, start)`` values."""
    return (range_obj.start, range_obj.stop)

class Dumper(Component):
    """
    Base class for dumpers.

    :param output: A file name, stream or :py:class:`pathlib.Path` object.
    """
    # Disable warnings for this half-implemented class.
    # pylint: disable=abstract-method

    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, output=sys.stdout):
        self.output = output

class Json(Dumper):
    """Dump the pipeline state in JSON."""
    def run(self, data, config=None, pipeline=None):
        """Dump pipeline state."""
        with Stream(self.output, "w") as out_fh:
            json.dump(data, out_fh, cls=JsonStateEncoder)
        return data

class Yaml(Dumper):
    """Dump pipeline state in YAML."""
    def run(self, data, config=None, pipeline=None):
        with Stream(self.output, "w") as out_fh:
            yaml.dump(data, out_fh, Dumper=YamlStateDumper)
