"""
Exceptions that may be raised from the :py:mod:`phyre_engine.component.pbs`
package.
"""

class UncompletedState(Exception):
    """
    Raised when a the output of a job does not contain a CompletedState object.
    """

    _ERR_MSG = "Pickle file `{}' did not contain a CompletedState object."

    def __init__(self, pickle):
        super().__init__(self._ERR_MSG.format(pickle))

class TooManyjobs(Exception):
    """
    Raised when trying to merge non-array jobs without a ``join_var`` key.
    """

    _ERR_MSG = "Could not merge jobs `{}' without a `join_var' key."

    def __init__(self, jobs):
        super().__init__(self._ERR_MSG.format([j.id for j in jobs]))
