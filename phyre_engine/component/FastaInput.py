from phyre_engine.component import Component
from phyre_engine.data.Sequence import Sequence

class FastaInput(Component):
    """Read a FASTA file as input and output the sequence."""

    def run(self, input: str):
        """Read the sequence from a FASTA file.

        Reads a single sequence from the file with the path given by input.
        If the file can not be read or contains multiple sequences an exception
        will be thrown.

        Args:
            input: Path of the FASTA file to read.

        Returns:
            A Sequence object representing the sequence from the FASTA file.

        Raises:
            IOError: Error reading the file.
            FastaInput.TooManySequencesError: The FASTA file contained multiple sequences.
        """

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
            return Sequence("".join(seq_lines))


    class TooManySequencesError(Exception):
        """Indicates too many sequences were present in a FASTA file."""
        pass
