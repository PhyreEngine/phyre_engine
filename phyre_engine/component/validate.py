"""Module containing classes for validating parsed sequences."""
from phyre_engine.component import Component
import Bio.Alphabet

class SeqValidator(Component):
    """
    Ensure that a valid sequence was given.

    This component raises a :py:exc:`.InvalidSeqError` if an unknown amino acid
    is present in the ``sequence`` field of the pipeline state.

    :param bool extended: If true, allow ``X``, ``B`` and ``Z`` amino acids.
    """
    REQUIRED = ['sequence']
    ADDS     = []
    REMOVES  = []

    def __init__(self, extended=False):
        if extended:
            self.alphabet = Bio.Alphabet.IUPAC.ExtendedIUPACProtein()
        else:
            self.alphabet = Bio.Alphabet.IUPAC.IUPACProtein()

    def run(self, data, config=None, pipeline=None):
        """Validate the ``sequence`` key.

        The ``sequence`` key of ``data`` should be a string containing
        single-letter amino acid codes. This component simply raises an
        exception if the sequence contains characters that are not contained
        within its alphabet.

        :raises SeqValidator.InvalidSeqError: Error describing the alphabet
                mismatch.
        """
        seq = self.get_vals(data)
        for letter in seq:
            if letter not in self.alphabet.letters:
                raise SeqValidator.InvalidSeqError(seq, self.alphabet)
        return data


    class InvalidSeqError(Exception):
        """Indicates that a sequence contained invalid letters."""
        def __init__(self, seq, alphabet):
            #Find which are disallowed. First, build a set of all letters in
            #this seq. Then, check each element of that set against the
            #alphabet.
            letters = set(seq)
            invalid = letters - set(alphabet.letters)

            super().__init__("Sequence contained the invalid characters {}",
                    invalid)
            self.invalid = invalid
            self.seq = seq
            self.alphabet = alphabet


