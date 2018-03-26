"""
Module containing components for reading sequences from files.

PhyreEngine uses `BioPython <http://biopython.org/>`_ to parse sequence files,
so any format supported by BioPython is also supported by PhyreEngine. The
sequence is stored as a string in the IUPAC alphabet in the ``sequence`` field
of the pipeline state, or the ``sequence`` field of each element in the ``templates``
list for :py:class:`~.ReadMultipleSequences`.

By default, these components also attempt to parse metadata from the input files.
This can be disabled by setting `metadata=False`. The metadata is read into the
``id``, ``name``, and ``description`` fields, and varies depending on the format
being read.
"""

from phyre_engine.component import Component
import Bio.SeqIO
from Bio.Alphabet import IUPAC

class ReadSingleSequence(Component):
    """
    Read a single sequence from the file ``input``. The sequence is stored
    in the ``sequence`` field, and the ``id``, ``name`` and ``description``
    are parsed from the sequence file and stored in the pipeline state.

    If any of the metadata fields are missing, the value `None` is used in its
    place. This behaviour is different to the BioPython defaults, which uses
    a special string (like ``<unknown id>``) to indicate a missing field.

    :param str seq_format: Format of the sequence file. This must be a format
        understood by BioPython's :py:mod:`Bio.SeqIO` module.

    :param bool metadata: If `False`, ignore all metadata supplied along with
        the sequence. Otherwise the ``name``, ``id`` and ``description`` fields
        in the pipeline state are parsed from the input file.
    """

    #: :param str input: Path of the FASTA file from which to read.
    REQUIRED = ['input']
    #: :param str sequence: Parsed sequence.
    #: :param str id: Sequence ID.
    #: :param str name: Sequence name, if any.
    #: :param str description: Sequence description, if any.
    ADDS = ['seq_record', 'name', 'id', 'description']
    REMOVES  = []

    def __init__(self, seq_format="fasta", metadata=True):
        self.seq_format = seq_format
        self.metadata = metadata

    def run(self, data, config=None, pipeline=None):
        """Read the sequence from a FASTA file.

        Reads a single sequence from the file with the path given by input.
        If the file can not be read or contains multiple sequences an exception
        will be thrown.

        :raises IOError: Error reading the file.
        :raises FastaInput.TooManySequencesError: The FASTA file contained
            multiple sequences.
        """

        input_file = self.get_vals(data)
        try:
            seq_record = Bio.SeqIO.read(
                input_file,
                format=self.seq_format,
                alphabet=IUPAC.protein)

            data["sequence"] = str(seq_record.seq)
            if self.metadata:
                data["id"] = seq_record.id
                data["name"] = seq_record.name
                data["description"] = seq_record.description

        except ValueError as e:
            raise self.TooManySequencesError() from e
        return data


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass


class ReadMultipleSequences(Component):
    """
    Create a template for each sequence in the ``input`` file.

    Adds the ``templates`` key to the pipeline data. Each template is a
    dictionary containing the fields added by
    :py:class:`~.ReadSingleSequenceRecord`.

    >>> from phyre_engine.component.input import ReadMultipleSequences
    >>> rms = ReadMultipleSequences("fasta")
    >>> rms.run({"input": "/path/to/fasta/file"})
    {"templates": [
        {
            "sequence": "AAA", "id": "foo", "name": "foo sequence",
            "description": "A description of this sequence."},
        }
        {
            "sequence": "AGX", "id": "bar", "name": "bar sequence",
            "description": "A description of this sequence."},
        }
        {
            "sequence": "AGG", "id": "baz", "name": "baz sequence",
            "description": "A description of this sequence."
        },
    ]}

    See :py:class:`~.ReadSingleSequenceRecord` for a description of the
    component parameters.
    """

    REQUIRED = ["input"]
    ADDS = ["templates"]
    REMOVES = []

    def __init__(self, seq_format="fasta", metadata=True):
        self.seq_format = seq_format
        self.metadata = metadata

    def run(self, data, config=None, pipeline=None):
        """Read multiple sequences from a FASTA file."""
        input_file = self.get_vals(data)
        templates = []
        for record in Bio.SeqIO.parse(input_file, "fasta"):
            template = {"sequence": str(record.seq)}
            if self.metadata:
                template["id"] = record.id
                template["name"] = record.name
                template["description"] = record.description
            templates.append(template)
        data["templates"] = templates
        return data
