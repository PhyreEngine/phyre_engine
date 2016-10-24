from phyre_engine.component import Component
import Bio.SeqIO

class FastaInput(Component):
    """Read a FASTA file as input and output the sequence."""
    REQUIRED = ['input']
    ADDS     = ['sequence']
    REMOVES  = []

    def run(self, data):
        """Read the sequence from a FASTA file.

        Reads a single sequence from the file with the path given by input.
        If the file can not be read or contains multiple sequences an exception
        will be thrown.

        Args:
            data: Key-value mapping of data. The following keys are required:
                `input`: Path of the FASTA file from which to read.

        Returns:
            A Bio.Seq object representing the sequence from the FASTA file.

        Raises:
            IOError: Error reading the file.
            FastaInput.TooManySequencesError: The FASTA file contained multiple sequences.
        """

        input = self.get_vals(data)
        with open(input, "r") as fasta:
            try:
                data['sequence'] = Bio.SeqIO.read(fasta, format="fasta")
            except ValueError as e:
                raise FastaInput.TooManySequencesError() from e
        return data


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass
