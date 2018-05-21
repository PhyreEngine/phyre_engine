"""
This module contains components for dumping the pipeline state in various
formats. By default, each dumper will write to standard output. The ``output``
parameter may be set to either a file name or output stream to determine where
the output will be written.
"""
import csv
import json.encoder
import pickle
import re
import sys

import Bio.SeqRecord

from phyre_engine.component.component import Component
from phyre_engine.tools.util import Stream, NamedTuple
import jmespath
from phyre_engine.tools.jmespath import JMESExtensions
import datetime


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
        elif isinstance(o, datetime.date):
            return serialise_date(o)
        elif isinstance(o, NamedTuple):
            return tuple(o)
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
        self.add_representer(datetime.date, self._serialise_date)
        self.add_multi_representer(NamedTuple, self._serialise_named_tuple)

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

    @staticmethod
    def _serialise_date(dumper, date):
        return dumper.represent_scalar(
            "tag:yaml.org,2002:str",
            serialise_date(date))

    @staticmethod
    def _serialise_named_tuple(dumper, tuple_obj):
        return dumper.represent_sequence(
            "tag:yaml.org,2002:seq",
            tuple(tuple_obj))

def serialise_seqrecord(record):
    """Use :py:mod:`Bio.SeqIO` to generate a FASTA sequence."""
    return record.format("fasta")

def serialise_range(range_obj):
    """Return a tuple of ``(stop, start)`` values."""
    return (range_obj.start, range_obj.stop)

def serialise_date(date):
    """Serialise a date object into ISO8601 format."""
    return date.isoformat()

class Dumper(Component):
    """
    Base class for dumpers.

    The state may be manipulated via the JMESPath expression `select_expr`. If
    `select_expr` is not supplied, the entire pipeline state will be dumped.

    :param output: A file name, stream or :py:class:`pathlib.Path` object.
    :param str select_expr: JMESPath expression to be applied to the pipeline
        state. The result of the expression is dumped.
    """
    # Disable warnings for this half-implemented class.
    # pylint: disable=abstract-method

    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, output=sys.stdout, select_expr=None):
        self.output = output
        self.select_expr = select_expr

    def _filter(self, data):
        """Filter data to contain the required fields."""
        if self.select_expr is None:
            return data

        opts = jmespath.Options(custom_functions=JMESExtensions(data))
        return jmespath.search(self.select_expr, data, opts)

class Json(Dumper):
    """Dump the pipeline state in JSON."""
    def run(self, data, config=None, pipeline=None):
        """Dump pipeline state."""
        with Stream(self.output, "w") as out_fh:
            json.dump(self._filter(data), out_fh, cls=JsonStateEncoder)
        return data

class Yaml(Dumper):
    """Dump pipeline state in YAML."""
    def run(self, data, config=None, pipeline=None):
        with Stream(self.output, "w") as out_fh:
            yaml.dump(self._filter(data), out_fh, Dumper=YamlStateDumper)
        return data

class Pickle(Dumper):
    """Dump pipeline state as a pickle."""
    def run(self, data, config=None, pipeline=None):
        with Stream(self.output, "wb") as out_fh:
            pickle.dump(self._filter(data), out_fh)
        return data

class Csv(Component):
    """
    Dump a list in the pipeline state in CSV format.

    :param str select_expr: JMESPath expression returning a list of
        dictionaries, all values of which will be written to the output file.

    :param file: File name or file-like object to write to. Defaults to standard
        output.

    :param bool header: Whether or not to include a header field giving the
        column names.

    :param str null_placeholder: Written for fields containing `None` of
        missing fields.

    :param \\*\\*csv_args: Extra arguments passed directly to the
        :py:func:`csv.writer` function.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, select_expr, output=sys.stdout, header=True,
                 null_placeholder="NA", **csv_args):
        self.select_expr = select_expr
        self.output = output
        self.header = header
        self.null_placeholder = null_placeholder
        self.csv_args = csv_args

    def _fill_placeholders(self, record):
        """
        Fill all `None`s in a dictionary with `self.null_placeholder`.
        """
        for key in list(record.keys()):
            if record[key] is None:
                record[key] = self.null_placeholder

    def run(self, data, config=None, pipeline=None):
        """Write CSV file."""
        jmespath_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        results = jmespath.search(self.select_expr, data, jmespath_opts)

        with Stream(self.output, "w") as csv_out:
            if not results:
                print("# No results", file=csv_out)
            else:
                writer = csv.DictWriter(
                    csv_out, sorted(results[0].keys()),
                    restval=self.null_placeholder)
                if self.header:
                    writer.writeheader()
                for record in results:
                    self._fill_placeholders(record)
                    writer.writerow(record)

        return data
