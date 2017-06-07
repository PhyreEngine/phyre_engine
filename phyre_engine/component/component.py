"""Module containing the base class of components."""
from abc import ABC, abstractmethod

class Component(ABC):
    """Base class for all component classes."""
    pass

    def get_vals(self, data):
        """Get all required values from a key-value mapping.

        This method returns a list containing the values corresponding to the
        keys given in the class variable list `REQUIRED`. It is designed to make
        it easy to unpack all required variables from a dictionary:

        >>> required_1, required_2 = self.get_vals(data)

        >>> required = self.get_vals(data)

        Note that if only one variable is required, a single variable is
        returned rather than a list of length 1. This is to allow simple
        unpacking, as shown in the example above.

        :param dict data: Key-value mapping from which to retrieve values.

        :returns:
            If only one item is required, a scalar corresponding to that item in
            the data field. If multiple items are required, a list containing
            the corresponding values is returned.
        """
        if len(type(self).REQUIRED) == 1:
            return data[type(self).REQUIRED[0]]
        else:
            return [data[x] for x in type(self).REQUIRED]

    @abstractmethod
    def run(self, data):
        """Run this component.

        This method must be implemented by all subclasses of Component.

        :param dict data: Dictionary containing all parameters produced by the
            pipeline.
        """
        pass #pragma: no cover
