"""Module containing the base class of components."""
from abc import ABC, abstractmethod
import logging
import phyre_engine

log = lambda: logging.getLogger(__package__)

class Component(ABC):
    """Base class for all component classes."""

    @property
    @abstractmethod
    def REQUIRED(self):
        pass

    @property
    @abstractmethod
    def ADDS(self):
        pass

    @property
    @abstractmethod
    def REMOVES(self):
        pass

    #: A string specifying a configuration section, or ``None`` if no external
    #: configuration is used.
    #:
    #: All members of the corresponding pipeline configuration section will be
    #: passed into the component constructor. Explicitly-set values will
    #: override values from the pipeline configuration.
    CONFIG_SECTION = None

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
    def run(self, data, config=None, pipeline=None):
        """Run this component.

        This method must be implemented by all subclasses of Component.

        :param dict data: Dictionary containing all parameters produced by the
            pipeline.
        """
        pass #pragma: no cover

class PipelineComponent(Component):
    """
    Base class for components that can be configured to execute another
    pipeline.

    :param pipeline: Either a :py:class:`phyre_engine.pipeline.Pipeline` class
        or a list of component names and arguments to be passed to
        :py:meth:`phyre_engine.pipeline.Pipeline.load`.

    :param dict config: Pipeline configuration to be used for the new pipeline.
    """
    # pylint: disable=abstract-method

    def __init__(self, pipeline, config=None):
        if not isinstance(pipeline, phyre_engine.Pipeline):
            pipeline = phyre_engine.Pipeline.load(pipeline)

        self.pipeline = pipeline
        self.config = config if config is not None else {}

class Map(PipelineComponent):
    """
    Apply a child pipeline for each element in an array contained within the
    pipeline state.

    :param str field: Field over which to iterate.

    .. seealso::

        :py:class:`.PipelineComponent`
            For extra class parameters.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, field, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.field = field

    def run(self, data, config=None, pipeline=None):
        """Iterate over given field, applying a child pipeline."""

        pipeline = self.pipeline
        pipe_output = []
        for item in data[self.field]:
            pipeline.start = item
            pipe_output.append(pipeline.run())
        data[self.field] = pipe_output
        return data

class Conditional(PipelineComponent):
    """
    Apply a child pipeline if the pipeline state contains a particular field
    and that field evaluates as ``True``.

    :param str field: Field to check.

    .. seealso::

        :py:class:`.PipelineComponent`
            For extra class parameters.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, field, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.field = field

    def run(self, data, config=None, pipeline=None):
        """Run child pipeline if `self.field` is ``True``."""
        if self.field in data and data[self.field]:
            pipeline = self.pipeline
            pipeline.start = data
            pipe_output = pipeline.run()
            data.update(pipe_output)
        return data

class TryCatch(PipelineComponent):
    """
    Run a child pipeline, catching and logging any exceptions that are raised.

    .. seealso::

        :py:class:`.PipelineComponent`
            For class parameters.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, logger=None, log_level=logging.ERROR, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logger if logger is not None else log()

        # Allow string log levels that we look up via getattr
        if isinstance(log_level, str):
            log_level = getattr(logging, log_level)
        self.log_level = log_level

    def run(self, data, config=None, pipeline=None):
        """Run child pipeline, ignoring errors."""
        try:
            pipeline = self.pipeline
            pipeline.start = data
            pipe_output = pipeline.run()
            return pipe_output
        except Exception as error:
            self.logger.log(
                self.log_level,
                "TryCatch: Ignoring exception %s",
                error, exc_info=True)
        return data
