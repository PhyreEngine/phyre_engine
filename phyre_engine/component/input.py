"""
Module containing components for reading sequences from files.

Each of the "Input" components in this module add a ``seq_record`` field to the
pipeline state. The object contained within the ``seq_record`` field is an
instance of BioPython's :py:class:`~Bio.SeqRecord.SeqRecord` class. The
:py:class:`.ConvertSeqRecord` component may be used to copy the sequence in the
``seq_record`` field to the ``sequence`` field as a plain Python string. The
:py:class:`.ConvertSeqRecord` component will also copy any metadata from the
``seq_record`` field into the pipeline state.

.. versionadded:: 0.1a2
   Early versions stored a :py:class:`~Bio.SeqRecord.SeqRecord` object in the
   ``sequence`` field. This behaviour was changed because it was felt to be too
   unexpected: nearly every other top-level element of the pipeline state is a
   simple Python datatype, and unexpectedly having an object in the ``sequence``
   field led to some confusion, especially when attempting to serialise the
   pipeline state.
"""
from phyre_engine.component import Component
import Bio.SeqIO
from Bio.Alphabet import IUPAC

class FastaInput(Component):
    """Read a FASTA file as input and output the sequence."""
    #: :param str input: Path of the FASTA file from which to read.
    REQUIRED = ['input']
    #: :param seq_record: Parsed sequence.
    #: :type seq_record: :class:`Bio.SeqRecord.SeqRecord`
    ADDS = ['seq_record']
    REMOVES  = []

    def run(self, data, config=None, pipeline=None):
        """Read the sequence from a FASTA file.

        Reads a single sequence from the file with the path given by input.
        If the file can not be read or contains multiple sequences an exception
        will be thrown.

        :raises IOError: Error reading the file.
        :raises FastaInput.TooManySequencesError: The FASTA file contained
            multiple sequences.
        """

        input = self.get_vals(data)
        with open(input, "r") as fasta:
            try:
                data['seq_record'] = Bio.SeqIO.read(fasta, format="fasta",
                        alphabet=IUPAC.protein)
            except ValueError as e:
                raise FastaInput.TooManySequencesError() from e
        return data


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass

class MultipleFastaInput(Component):
    """
    Create a template for each sequence in a FASTA file.

    Adds the ``templates`` key to the pipeline data. Each template is a
    dictinoary containing a ``sequence'' key point to a
    :py:class:`Bio.PDB.SeqRecord` object.
    """

    REQUIRED = ["input"]
    ADDS = ["templates"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Read multiple sequences from a FASTA file."""
        input = self.get_vals(data)
        templates = []
        for record in Bio.SeqIO.parse(input, "fasta"):
            templates.append({"seq_record": record})
        data["templates"] = templates
        return data

class ConvertSeqRecord(Component):
    """
    Convert the ``seq_record`` field to a simple text string in the ``sequence``
    field. This component will also convert the ``name``, ``id`` and
    ``description`` fields of the ``seq_record`` to pipeline state attributes.

    The ``sequence`` field will consist of one-letter amino acid codes.

    >>> from io import StringIO
    >>> import Bio.SeqIO
    >>> from phyre_engine.component.input import ConvertSeqRecord
    >>> fasta = ">FOO BAR\nAAAGGG\n"
    >>> with StringIO(fasta) as fasta_in:
    ...     seq_record = Bio.SeqIO.read(fasta_in, 'fasta')
    >>> results = ConvertSeqRecord().run({'seq_record': seq_record})
    >>> results.sequence
    'AAAGGG'
    >>> results.name
    'FOO'
    >>> results.id
    'FOO'
    >>> results.description
    'FOO BAR'
    """
    REQUIRED = ["seq_record"]
    ADDS = ["sequence", "name", "id", "description"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        seq_record = self.get_vals(data)
        data["sequence"] = str(seq_record.seq)
        data["id"] = seq_record.id
        data["name"] = seq_record.name
        data["description"] = seq_record.description
        return data

