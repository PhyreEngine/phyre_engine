"""Module containing classes for validating parsed sequences."""
from phyre_engine.component import Component
import Bio.Alphabet.IUPAC

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

class SeqLenFilter(Component):
    """
    Raise an exception when a sequence has a length outside of the specified
    range. Sequences of length :math:`l` will be retained by this filter if
    :math:`L_\text{min} \le l L_\text{max}`, where :math:`L_\text{min}` and
    :math:`L_\text{max}` are the `min` and `max` parameters passed to this
    filter.

    :param int min: Minimum sequence length (inclusive).
    :param int max: Maximum sequence length (inclusive).
    """

    class SeqLenError(RuntimeError):
        """Raised when a sequence length is invalid."""

        _ERR_MSG = "Sequence of length {} outside of length range {}-{}"
        def __init__(self, seq_len, min_len, max_len):
            super().__init__(
                self._ERR_MSG.format(seq_len, min_len, max_len))

    REQUIRED = ["sequence"]
    ADDS = []
    REMOVES = []

    def __init__(self, min_len, max_len):
        self.seq_range = range(min_len, max_len + 1)

    def run(self, data, config=None, pipeline=None):
        """
        Filter sequences by length.

        :raises phyre_engine.component.validated.SeqLenFilter.SeqLenError: When
            the sequence falls outside of allowed range.
        """
        sequence = self.get_vals(data)
        if len(sequence) not in self.seq_range:
            raise self.SeqLenError(
                len(sequence),
                self.seq_range.start, self.seq_range.stop - 1)
        return data
