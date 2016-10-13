from phyre_engine.component import Component
from phyre_engine.data.Sequence import Sequence

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
            A Sequence object representing the sequence from the FASTA file.

        Raises:
            IOError: Error reading the file.
            FastaInput.TooManySequencesError: The FASTA file contained multiple sequences.
        """

        input = self.get_vals(data)
        with open(input, "r") as fasta:
            #Count the number of identifiers so we can throw an errror if
            #there is more than one.
            identifiers = 0
            seq_lines = []
            for line in fasta:
                if line.startswith(">"):
                    identifiers = identifiers + 1
                    if identifiers > 1:
                        raise self.TooManySequencesError()
                else:
                    seq_lines.append(line.strip())
            data['sequence'] = Sequence("".join(seq_lines))
        return data


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass
