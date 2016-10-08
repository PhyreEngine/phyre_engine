import data

class Sequence:
    """Class representing an amino acid sequence.

    The string representation of objects of this class will give the textual
    representation of the sequence. Objects of this type should be treated as
    immutable, though this is not enforced.
    """

    def __init__(self, sequence, allowed_extra=set()):
        """Create a new Sequence from sequence.

        If the supplied sequence contains invalid characters, an exception will
        be raised. This method calls validate() internally.

        Args:
            sequence:A string containing the single-letter representation of
                an amino acid sequence.
            allowed_extra: Any allowed amino acids outside of the usual 21.
        """

        self.allowed_extra = allowed_extra
        self.validate(sequence)
        self._sequence = sequence

    def validate(self, sequence):
        """Validate sequence, throwing error if invalid.

        Args:
            sequence: Sequence to check.

        Raises:
            Sequence.UnknownAminoAcidError: The sequence contained disallowed
                amino acids.
            Sequence.IsNucleotideError: The sequence appears to be a DNA
                string.
        """
        for i, aa in enumerate(sequence):
            if aa not in data.AMINO_ACIDS and aa not in self.allowed_extra:
                raise Sequence.UnknownAminoAcidError(aa, i)

    def __str__(self):
        """Gives the string representation of the sequence."""
        return self._sequence

    class UnknownAminoAcidError(Exception):
        """Error thrown when an unknown amino acidn type is encountered."""
        def __init__(self, residue, pos):
            err_msg = "Unknown residue type '{}' at residue {}"
            super().__init__(err_msg.format(residue, pos))

