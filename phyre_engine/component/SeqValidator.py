"""Module containing classes for validating parsed sequences."""
from phyre_engine.component import Component
import Bio.Alphabet

class SeqValidator(Component):
    """Ensure that a valid sequence was given."""
    #: :param sequence: Sequence to be checked.
    #: :type sequence: :class:`Bio.SeqRecord`
    REQUIRED = ['sequence']
    ADDS     = []
    REMOVES  = []

    def run(self, data, config=None, pipeline=None):
        """Validate the ``sequence`` key.

        The ``sequence`` key of ``data`` should be a ``Bio.SeqRecord``-like
        object associated with an alphabet. This component simply raises an
        exception if the sequence contains characters that are not contained
        within its alphabet.

        :raises SeqValidator.InvalidSeqError: Error describing the alphabet
                mismatch.
        """

        seq = self.get_vals(data)
        if not Bio.Alphabet._verify_alphabet(seq.seq):
            raise SeqValidator.InvalidSeqError(seq)
        return data


    class InvalidSeqError(Exception):
        """Indicates that a sequence contained invalid letters."""
        def __init__(self, seq):
            #Find which are disallowed. First, build a set of all letters in
            #this seq. Then, check each element of that set against the
            #alphabet.
            letters = set(str(seq.seq))
            invalid = set()
            for letter in letters:
                if not letter in seq.seq.alphabet.letters:
                    invalid.add(letter)
            super().__init__("Sequence contained the invalid characters {}",
                    invalid)
            self.invalid = invalid
            self.seq     = seq


