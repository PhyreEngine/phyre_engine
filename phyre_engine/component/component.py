"""Module containing the base class of components."""
from abc import ABCMeta, abstractmethod
import logging
import phyre_engine
import copy
import enum

import jmespath

from phyre_engine.tools.jmespath import JMESExtensions
from phyre_engine.tools.util import apply_dotted_key, deep_merge

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

    @classmethod
    def config(cls, params, pipeline_config):
        """
        Combine pipeline configuration options with constructor parameters.

        The pipeline configuration provides a means of configuring multiple
        similar components at once. The default strategy, implemented by this
        class method, is to look up a configuration section in the pipeline
        configuration, then update that dictionary with the component-specific
        parameters. These parameters are then passed to the constructor of the
        component as named arguments. The configuration section is chosen
        according to the `CONFIG_SECTION` class attribute.

        For example, consider this toy pipeline definition:

        .. code-block:: yaml

            pipeline:
              config:
                thing:
                  foo: 3.141
                  bar: {"x": 1, "y": 2}
              components:
              - .foo.Thing:
                  foo: 3
                  bar: {"x": 3}

        If we assume that the ``.foo.Thing`` component has ``CONFIG_SECTION =
        "thing"``, then the component will be instantiated with the arguments
        ``foo=3, bar={"x": 3, "y": 2}``. Notice that the parameters supplied
        specifically to the component override the parameters supplied in the
        pipeline configuration. This applies even in "child" dictionaries: the
        configurations are merged using
        :py:func:`phyre_engine.tools.util.deep_merge`.


        :param dict[str, any] params: Keyword parameters for the component
            constructor.

        :param phyre_engine.pipeline.PipelineConfig pipeline_config:
            Global pipeline configuration.

        .. note::

            Components that require more complex combinations of the
            configurations should override this method.
        """
        if params is None:
            params = {}
        if cls.CONFIG_SECTION is None:
            return params
        config_section = pipeline_config.section(cls.CONFIG_SECTION)
        return config_section.merge_params(params)

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
        keys given in the list property `REQUIRED`. It is designed to make
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
        if len(self.REQUIRED) == 1:
            return data[self.REQUIRED[0]]
        else:
            return [data[x] for x in self.REQUIRED]

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

    :param components: List of component names and arguments to be passed to
        :py:meth:`phyre_engine.pipeline.Pipeline.load`. *Either* this or
        `pipeline` must be supplied. This parameter exists purely for the sake
        of convenience: it is equivalent to supplying the `pipeline` parameter
        with the same components and no extra configuration.

    :param str config_mode: Value describing how to combine the parent and
        child pipeline configurations. This parameter may be passed as either a
        string (corresponding to an entry in the
        :py:class:`.ConfigurationPreference` enum) or as an enum value.

    :param bool lazy_load: Lazily load pipelines passed as dicts. If this is
        `True`, the child pipeline will be created when the `run` method is
        called, which will allow it to use the runtime configuration. If lazy
        loading is disabled, the pipeline will be loaded when this component is
        instantiated, which will cause an exception if configuration values are
        missing. If your pipeline is completely static, set this to `False` so
        that any errors are made obvious as quickly as possible.
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

    def __init__(self, pipeline=None, config_mode="child", lazy_load=True,
                 components=None):
        self.config_mode = type(self).ConfigurationPreference(config_mode)

        if (pipeline, components).count(None) != 1:
            raise ValueError(
                "Supply one and only one of 'pipeline' or 'components'")
        if components is not None:
            pipeline = {"components": components}

        if not lazy_load and not isinstance(pipeline, phyre_engine.Pipeline):
            self._pipeline = phyre_engine.Pipeline.load(pipeline)
        else:
            self._pipeline = pipeline

    def pipeline(self, runtime_config, start=None):
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
            pipeline_defn = self.pipeline_definition(runtime_config, start)
            return phyre_engine.Pipeline.load(pipeline_defn)
        return self._pipeline

    def pipeline_definition(self, runtime_config, start=None):
        """
        Retrieve the definition of the child pipeline, as would be specified in
        a YAML pipeline file.

        :raises TypeError: If the pipeline was passed as an instance of
        :py:class:`phyre_engine.pipeline.Pipeline`.
        """
        if isinstance(self._pipeline, phyre_engine.Pipeline):
            raise TypeError(
                ("Cannot get pipeline definition for instance of Pipeline "
                 "class."))

        if runtime_config is None:
            runtime_config = {}

        pipeline_definition = copy.deepcopy(self._pipeline)
        pipeline_definition["config"] = self.combine_configs(
            runtime_config,
            pipeline_definition.get("config", {}))

        if start is not None:
            pipeline_definition["start"] = start
        return pipeline_definition

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
    ADDS = []
    REMOVES = []

    @property
    def REQUIRED(self):
        return [self.field]

    def __init__(self, field, *args, discard=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.field = field
        self.discard = discard

    def run(self, data, config=None, pipeline=None):
        """Iterate over given field, applying a child pipeline."""

        pipeline = self.pipeline(config)
        pipe_output = []
        for i, item in enumerate(data[self.field]):
            self.logger.debug("Runing pipeline %d / %d",
                              i, len(data[self.field]))
            pipeline.start = item
            pipeline_results = pipeline.run()
            if pipeline_results is not None and not self.discard:
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

    def __init__(self, *args, keep=(), **kwargs):
        super().__init__(*args, **kwargs)
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


class ConfigLoader(PipelineComponent):
    """
    Convenience component for loading pipeline state into sub-pipeline
    configuration.

    It can sometimes be useful to use part of the pipeline state as
    configuration parameters to some components. This component provides this
    capability by converting parts of the pipeline state into configuration
    according to `config_map`, then loading and executing the child pipeline
    with the current pipeline state.

    For example, let's say we want to take a slice from the ``templates`` list.
    We can do this easily using the
    :py:class:`phyre_engine.component.jmespath.Replace` component by passing
    the slice as the `value_expr`. To take the first ten templates, we would
    write the following:

    .. code-block:: yaml

        # In pipeline configuration file
        pipeline:
          components:
          # ...
          - .jmespath.Replace:
              select_expr: templates
              value_expr: [0:10]

    This slice could be configured on the command line with the parameter
    ``--config 'jmespath.value_expr:[0:10]'``. However, this is quite obtuse:
    it's not immediately clear what the ``value_expr`` is selecting. If we
    wanted to configure two separate JMESPath components from the command line,
    we would be out of luck.

    This component will transfer fields from the pipeline state into the child
    configuration according to a predefined set of aliases. For the example
    given above, we could use the ``slice`` field of the pipeline state as
    follows:

    .. code-block:: yaml

        # In pipeline configuraton file
        pipeline:
          components:
          # ...
          - .component.ConfigLoader:
              mapping:
                slice: jmespath.value_expr
              components:
              - .jmespath.Replace:
                  select_expr: templates

    This pipeline can then be run with the parameters ``--start
    'slice:[0:10]'``. The `ConfigLoader` component will transfer the ``slice``
    field into the ``value_expr`` field of the ``jmespath`` section of the
    sub-pipeline config.

    The value of each field in the ``mapping`` parameter is a string, with
    field names separated by dots. Each dot indicates drilling down into a
    new dictionary. That is, ``jmespath.value_expr`` is the same as
    ``config["jmespath"]["value_expr"]``.

    :param dict mapping: Mapping of pipeline state fields to configuration
        fields. The keys of this mapping are JMESPath expressions applied to
        the pipeline state, which must return the data to be transferred into
        the pipeline configuration. The values of each field are similar to
        JMESPath expressions, but more limited: they may only be dot-separated
        strings, each drilling further into a dictionary.
    """
    @property
    def REQUIRED(self):
        return list(self.mapping.keys())

    ADDS = []
    REMOVES = []

    def __init__(self, mapping, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mapping = mapping

    def generate_config(self, data, config):
        """
        Generate child pipeline configuration from runtime configuration
        and pipeline state.
        """
        config = config if config is not None else {}
        jmes_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        for search_term, config_location in self.mapping.items():
            state_value = jmespath.search(search_term, data, jmes_opts)
            apply_dotted_key(config, config_location, state_value)
        return config

    def run(self, data, config=None, pipeline=None):
        """Alias pipeline state into child configuraton."""
        config = self.generate_config(data, config)
        pipeline = self.pipeline(config, data)
        return pipeline.run()
