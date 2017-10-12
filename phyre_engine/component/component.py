"""Module containing the base class of components."""
from abc import ABCMeta, abstractmethod
import logging
import phyre_engine
import copy
import enum

class ComponentMeta(ABCMeta):
    """Metaclass for components.

    Currently, this class only exists to provide the :py:meth:`.qualname` class
    property.
    """

    @property
    def qualname(cls):
        """
        Qualified name of this component. See :py:meth:`.Component.qualname`.
        """
        qualname = cls.__qualname__
        modname = cls.__module__
        if "." in qualname:
            return (modname, qualname)
        else:
            return ".".join([modname, qualname])

class Component(metaclass=ComponentMeta):
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

    @property
    def logger(self):
        """Get a logger named for this component."""
        logger_name = ".".join((type(self).__module__, type(self).__qualname__))
        return logging.getLogger(logger_name)

    @property
    def qualname(self):
        """Fully qualified name of this component.

        For non-nested classes---that is, classes defined at the module
        level---this will return a single string. For nested classes, it will
        return a tuple containing the module and the qualified name of the
        class. The return value of this function is suitable for use in the
        ``components`` list of the dictionary passed to
        :py:meth:`phyre_engine.pipeline.Pipeline.load`.

        Returns:
            Either a single string or a tuple containing the fully qualified
            name of this component.
        """
        return type(self).qualname

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

    :param str config_mode: Value describing how to combine the parent and
        child pipeline configurations. This parameter may be passed as either a
        string (corresponding to an entry in the
        :py:class:`.ConfigurationPreference` enum) or as an enum value.
    """
    # pylint: disable=abstract-method

    class ConfigurationPreference(enum.Enum):
        """Order in which to load pipeline configurations."""

        #: Treat the runtime configuration (the configuration of the parent
        #: pipeline) as the "base" configuration, and override it with the
        #: values in the ``config`` section of the child pipeline.
        PREFER_CHILD = "child"

        #: Use the child pipeline (passed with the `pipeline` parameter of this
        #: class) as the base configuration, and override it with properties
        #: from the runtime (parent) configuration.
        PREFER_PARENT = "parent"

        #: Discard the runtime configuration completely, and only retain the
        #: configuration of the child.
        DISCARD_PARENT = "discard"

    def __init__(self, pipeline, config_mode="child"):
        self._pipeline = pipeline
        self.config_mode = type(self).ConfigurationPreference(config_mode)

    def pipeline(self, runtime_config):
        """
        Retrieve the child pipeline.

        If the child pipeline is an instance of
        :py:class:`phyre_engine.pipeline.Pipeline`, it will be returned as-is
        without the configuration being altered. If the child pipeline is a
        dictionary, it will be loaded by
        :py:meth:`phyre_engine.pipeline.Pipeline.load`. See
        :py:class`.ConfigurationPreference` for details on how the pipeline
        configurations will be merged.

        :return: Loaded pipeline.
        :rtype: :py:class:`phyre_engine.pipeline.Pipeline`
        """
        if not isinstance(self._pipeline, phyre_engine.Pipeline):
            if runtime_config is None:
                runtime_config = {}

            pipeline_definition = copy.deepcopy(self._pipeline)
            pipeline_definition["config"] = self.combine_configs(
                runtime_config,
                pipeline_definition.get("config", {}))

            return phyre_engine.Pipeline.load(pipeline_definition)
        return self._pipeline


    def combine_configs(self, parent_config, child_config):
        """
        Combine parent and child configurations.

        See :py:class:`.ConfigurationPreference` for details on how the
        configurations will be combined.

        :param dict parent_config: Configuration of the runtime (parent)
            pipeline.

        :param dict child_config: Configuration of the child pipeline.
        """
        # Alias to reduce typing
        conf_enum = PipelineComponent.ConfigurationPreference

        if self.config_mode == conf_enum.DISCARD_PARENT:
            return child_config
        elif self.config_mode == conf_enum.PREFER_CHILD:
            base_config = copy.deepcopy(parent_config)
            base_config.update(child_config)
            return base_config
        elif self.config_mode == conf_enum.PREFER_PARENT:
            base_config = copy.deepcopy(child_config)
            base_config.update(parent_config)
            return base_config
        else:
            raise ValueError(
                "Invalid value '{}' for config_mode.".format(self.config_mode))

class Map(PipelineComponent):
    """
    Apply a child pipeline for each element in an array contained within the
    pipeline state. If the child pipeline returns `None`, it will not be
    included in the results. If it returns a list, then each element of the list
    is added to the pipeline state.

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

        pipeline = self.pipeline(config)
        pipe_output = []
        for item in data[self.field]:
            pipeline.start = item
            pipeline_results = pipeline.run()
            if pipeline_results is not None:
                if isinstance(pipeline_results, list):
                    pipe_output.extend(pipeline_results)
                else:
                    pipe_output.append(pipeline_results)
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
            pipeline = self.pipeline(config)
            pipeline.start = data
            pipe_output = pipeline.run()
            data.update(pipe_output)
        return data

class TryCatch(PipelineComponent):
    """
    Run a child pipeline, catching and logging any exceptions that are raised.

    :param bool pass_through: If an exception is raised in the child pipeline
        and this is `True`, the original pipeline state is retained. Otherwise,
        if an exception is raised and this is `False` this component will return
        `None` as the pipeline state.

    :param logger: Optional logger object. The default value is the module-level
        logger for this module.

    :param int log_level: The level at which to log events. This may
        alternatively be specified as a string.

    .. seealso::

        :py:class:`.PipelineComponent`
            For class parameters.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(
            self, pass_through=False, logger=None, log_level=logging.ERROR,
            *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pass_through = pass_through
        self.err_logger = logger if logger is not None else self.logger

        # Allow string log levels that we look up via getattr
        if isinstance(log_level, str):
            log_level = getattr(logging, log_level)
        self.log_level = log_level

    def run(self, data, config=None, pipeline=None):
        """Run child pipeline, ignoring errors."""
        try:
            pipeline = self.pipeline(config)
            pipeline.start = data
            pipe_output = pipeline.run()
            return pipe_output
        except Exception as error:
            self.err_logger.log(
                self.log_level,
                "TryCatch: Ignoring exception %s",
                error, exc_info=True)
            if self.pass_through:
                return data
            return None

class Branch(PipelineComponent):
    """
    Run a pipeline on a copy of the pipeline state and then restore the original
    state. This can be useful when you want to run components with side-effects
    but don't want to actually alter the pipeline state.

    For example, you might wish to run a component that sanitises the PDB file
    pointed to by the ``structure`` key and then write some information derived
    from the sanitised file. The branch could alter the ``structure`` key to
    point at the new sanitised file without affecting the main pipeline.

    To keep certain fields, set the `keep` parameter to a list of those fields.

    :param list[str] keep: List of fields to keep from the branched pipeline.

    .. note::

        This component only operates on the pipeline state. Side effects are
        not tracked and will not be reverted when the branch finishes.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, pipeline, keep=(), *args, **kwargs):
        super().__init__(pipeline, *args, **kwargs)
        self.keep = keep

    def run(self, data, config=None, pipeline=None):
        """Run a branch of the pipeline state."""
        pipeline = self.pipeline(config)
        pipeline.start = copy.deepcopy(data)
        branch_results = pipeline.run()
        for field in self.keep:
            if field in branch_results:
                data[field] = branch_results[field]
        return data
