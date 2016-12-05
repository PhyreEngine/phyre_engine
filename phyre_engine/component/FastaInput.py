"""Module containing components for reading sequences from files."""
from phyre_engine.component import Component
import Bio.SeqIO
from Bio.Alphabet import IUPAC

class FastaInput(Component):
    """Read a FASTA file as input and output the sequence."""
    #: :param str input: Path of the FASTA file from which to read.
    REQUIRED = ['input']
    #: :param sequence: Parsed sequence.
    #: :type sequence: :class:`Bio.SeqRecord`
    ADDS     = ['sequence']
    REMOVES  = []

    def run(self, data):
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
                data['sequence'] = Bio.SeqIO.read(fasta, format="fasta",
                        alphabet=IUPAC.protein)
            except ValueError as e:
                raise FastaInput.TooManySequencesError() from e
        return data


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass
