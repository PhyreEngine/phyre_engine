"""Provides classes for building, checking and executing pipelines."""

import copy
import enum
import importlib
import json
from collections import namedtuple
import pickle
import pathlib
import time
import logging

#: Allows the state of the pipeline to be saved.
#: :param int current_component: Index of the currently-running component. This
#:     component will be started if the pipeline is resumed fromn this
#:     checkpoint.
#: :param dict state: Data in the pipeline.
Checkpoint = namedtuple("Checkpoint", ["current_component", "state"])

#: A point in time.
#: :param str component: Name of the most recently run component, or `None`.
#: :param float wall: Wall-clock time in fractional seconds since the epoch.
#: :param float process: Process time in fractional seconds.
TimeInfo = namedtuple("TimeInfo", "component wall process")

#: Default namespace for pipelines.
#:
#: Components may be loaded from this namespace by specifying them with a
#: leading dot: ``.foo.Bar`` is equivalent to
#: ``phyre_engine.component.foo.Bar``. This namespace may be overriden by
#: the ``namespace`` top-level key.
DEFAULT_NAMESPACE = "phyre_engine.component"

class Pipeline:
    """Pipeline containing a list of components to be executed.

    A "pipeline" is simply a set of components, each of which should implement
    the Component abstract base class from phyre_engine.component. Each
    component has a "run" method that accepts a blob (simply implemented as a
    dictionary) of key-value pairs, performs some computation, modifies and
    returns the blob.

    Each component requires a certain set of keys on which to operate. For
    example, a component that reads a sequence file might require a "file" key
    giving the file name. Each component may also optionally add and remove
    keys. In the case of our sequence file, we might add a "sequence" key and
    could remove the "file" key. The blob may contain many other keys;
    components should not modify keys that they do not care about.

    To allow for programmatic validation of pipelines, each component defines
    the set of keys that it requires, adds and removes. This class provides
    methods to validate a pipeline.


    Often, the start of a pipeline will require an input (for example, a
    sequence or the name of a sequence file). To supply this, the parameter
    "start" may be specified. Without this, validation will fail if the first
    component in the list requires any keys.

    :param components: Optional list of component objects.
    :type components: List of :class:`phyre_engine.component.Component` objects
    :param dict start: Starting elements to include in the key-value map.
    :param str checkpoint: Path of a checkpoint file in which to load and save
        pipeline state.
    :param str statusfile: If set, write a JSON-encoded status file to this
        location. See :py:meth:`.update_statusfile`.
    """

    class CurrentComponent(enum.Enum):
        """
        Passed to :py:meth:`.update_statusfile` to indicate that the pipeline
        is either initialising or finished.
        """
        INITIALISING = "INITIALISING"
        FINISHED = "FINISHED"

    def __init__(self, components=None, start=None, checkpoint=None,
                 config=None, statusfile=None):
        """Initialise a new pipeline with an optional list of components."""
        self.components = components if components else []
        self.start = start if start else {}
        self.checkpoint = checkpoint
        self.config = config
        self.statusfile = statusfile
        self.logger = logging.getLogger(
            self.__module__ + "." + type(self).__qualname__)

    def update_statusfile(self, current_component):
        """
        Write a JSON-encoded file that contains a description of all the
        components that have been run so far.

        If the pipeline was not configured to use a statusfile (by setting the
        `statusfile` parameter), nothing is written.

        The status file is JSON-encoded, and contains the following parameters:

        ``components``
            List of components in the pipeline. Includes two extra components,
            ``INITIALISING`` and ``FINISHED``. These are placed at the
            beginning and end of the list, respectively, and indicate that
            the pipeline is either being set up or is complete.

        ``index``
            The index of the component (beginning from zero) currently being
            processed by the pipeline.
        """
        if self.statusfile is None:
            return

        components = (
            [self.CurrentComponent.INITIALISING.value]
            + [type(cmpt).__name__ for cmpt in self.components]
            + [self.CurrentComponent.FINISHED.value])

        if current_component in self.CurrentComponent:
            if current_component == self.CurrentComponent.INITIALISING:
                current_component = 0
            else:
                current_component = len(components) - 1
        else:
            # Compensate for leading "INITIALISING"
            current_component = current_component + 1

        with open(self.statusfile, "w") as status_out:
            json.dump({
                "index": current_component,
                "components": components
                }, status_out, indent=2)


    def validate(self):
        """Validate that the inputs and outptuts of each component match.

        This method walks through the pipeline, keeping track of the keys that
        are required, added, and removed for each component. An exception is
        raised if the components do not match.

        :raises Pipeline.ValidationError: The pipeline could not be validated.
        """

        keys = set(self.start.keys())
        for component in self.components:
            #Build a list of missing keys rather than failing on the first.
            #Hopefully this will be useful for debugging.
            missing = []
            for reqd in type(component).REQUIRED:
                if reqd not in keys:
                    missing.append(reqd)
            if len(missing) > 0:
                raise Pipeline.ValidationError(component, missing)

            for added in type(component).ADDS:
                keys.add(added)

            for removed in type(component).REMOVES:
                if removed in keys:
                    keys.remove(removed)


    def load_checkpoint(self):
        """
        Load checkpoint file. Returns the deserialised checkpoint, or `None` if
        the checkpoint could not be found or should not be used.
        """
        if self.checkpoint is not None:
            checkpoint_file = pathlib.Path(self.checkpoint)
            if checkpoint_file.exists():
                self.logger.info("Checkpoint found at %s.", checkpoint_file)
                with checkpoint_file.open("rb") as check_in:
                    chk = pickle.load(check_in)
                if chk.current_component >= len(self.components):
                    self.logger.info("All components complete. Finishing now.")
                else:
                    component = self.components[chk.current_component]
                    self.logger.info(
                        "Checkpoint found at %s. Skipping to component %d (%s)",
                        checkpoint_file, chk.current_component,
                        self._component_name(component))
                return chk
        return None

    def save_checkpoint(self, current_component, state):
        """
        If the pipeline is configured to use checkpoints, save the current
        state.
        """
        if self.checkpoint is not None:
            self.logger.debug("Saving checkpoint %s with current component %d",
                              self.checkpoint, current_component)
            checkpoint = Checkpoint(current_component, state)
            checkpoint_file = pathlib.Path(self.checkpoint)
            with checkpoint_file.open("wb") as check_out:
                pickle.dump(checkpoint, check_out)

    def store_timings(self, state, component):
        """
        Save current time in the ``timer`` field of the pipeline state.
        This saves both the current wall-clock time and the "process" time.
        The ``timer`` field will contain a list of dictionaries, each containing
        the same fields as :py:class:`.TimeInfo`.

        If the ``timer`` field in the pipeline configuration is not a `True`
        value, then nothing is added to the pipeline state.

        If the state is not a dictionary, then nothing is added. Components are
        allowed to return `None` (indicating a non-fatal error) or a list
        (indicating an expansion of the state) when they are run as part of a
        sub-pipeline.
        """
        if state is None or not isinstance(state, dict):
            return

        if self.config is not None and self.config.get("timer", False):
            if "timer" not in state:
                state["timer"] = []

            if component is not None:
                component_name = self._component_name(component)
            else:
                component_name = None

            current_time = TimeInfo(component_name, time.time(),
                                    time.process_time())
            state["timer"].append(dict(current_time._asdict()))


    def _component_name(self, component):
        """A pretty name for a component based off the module and class name."""
        qualname = component.qualname
        if isinstance(qualname, tuple):
            qualname = ".".join(qualname)
        return qualname

    def run(self, start_index=0):
        """Run this pipeline, executing each component in turn.

        :returns: Modified data dictionary, with results added by components.

        :raises Pipeline.ValidationError: If a component does not fulfil its
            promises and add a key that the next component requires.
        """
        # Write "INITIALISING" state to status file.
        self.update_statusfile(self.CurrentComponent.INITIALISING)

        # Load checkpoint or start a new pipeline if no checkpoint exists.
        checkpoint = self.load_checkpoint()
        if checkpoint is None:
            checkpoint = Checkpoint(start_index, copy.copy(self.start))
        current_component, state = checkpoint

        # Store timing data. We pass "None" as a component to indicate that this
        # is the reference datum. This is necessary because the process time has
        # no absolute meaning and must be compared to a reference.
        self.store_timings(state, None)

        # Start from "current_component" to enable restarting from a checkpoint.
        for cmpt in self.components[current_component:]:
            self.logger.info("Running component %d: %s",
                             current_component, self._component_name(cmpt))

            # Raise an exception if the required keys are not present in the
            # pipeline state.
            self.validate_runtime(state, cmpt)

            # Write progress information to a file if enabled.
            self.update_statusfile(current_component)

            # Actually run the component
            state = cmpt.run(state, self.config, self)

            # At this point, "state" can be None. This can happen, for example,
            # when a component runs a child pipeline. If a component in the
            # child pipeline wants to indicate a non-fatal (i.e. non-exception)
            # error, it can return None. The parent component can then decide
            # how to handle that. If a component in the top-level pipeline
            # returns None, that is an error that can be handled by the run
            # module.
            if state is None:
                return None

            # Record the finish time of this component.
            self.store_timings(state, cmpt)

            # Finally, save the checkpoint file.
            current_component += 1
            self.save_checkpoint(current_component, state)

        # Save a statfile with None as the current component
        self.update_statusfile(self.CurrentComponent.FINISHED)

        return state


    def validate_runtime(self, blob, component):
        """Check that the required keys are present in the key-value blob.

        This function is called at runtime when running each component of the
        pipeline to check that the promised keys are actually present in the
        data blob.

        :param dict blob: Key-value blob to check.
        :param component: Component from which the required keys are taken.

        :raises Pipeline.ValidationError: If a component does not fulfil its
            promises and add a key that the next component requires.
        """
        missing = []
        for reqd in type(component).REQUIRED:
            if reqd not in blob:
                missing.append(reqd)

        if len(missing) > 0:
            raise Pipeline.ValidationError(component, missing, blob)


    @staticmethod
    def _load_component(dotted_name, arg_list, config, namespace):
        """
        Load a component specified by either a dotted string or a module name
        and dotted string.

        The parameter ``dotted_name`` may be a tuple or a string. If it is a
        tuple, then the first element is assumed to be a module and the second
        element the class name. Nested classes are allowed.

        If ``dotted_name`` is a string, it is assumed that the class name is the
        final component and the module name everything before that. Nested
        classes may not be used.

        As shorthand, components may be specified with a leading dot, in which
        case `namespace` is prepended.
        """

        if isinstance(dotted_name, str):
            mod_name, cls_name = dotted_name.rsplit(".", maxsplit=1)
        else:
            mod_name, cls_name = dotted_name

        if mod_name.startswith("."):
            mod_name = namespace + mod_name

        module = importlib.import_module(mod_name)
        nested_cls_names = cls_name.split(".")
        component_cls = getattr(module, nested_cls_names[0])
        for nested_cls_name in nested_cls_names[1:]:
            component_cls = getattr(component_cls, nested_cls_name)

        # Collect *args and **kwargs
        args = []
        kwargs = {}

        if component_cls.CONFIG_SECTION is not None:
            kwargs.update(config.get(component_cls.CONFIG_SECTION, {}))

        for arg in arg_list if arg_list else []:
            if isinstance(arg, dict):
                kwargs.update(arg)
            else:
                args.append(arg)
        return component_cls(*args, **kwargs)

    @classmethod
    def load(cls, pipeline_dict):
        """
        Create a pipeline from a dictionary describing a pipeline.

        The dictionary must contain the top-level field ``components``, which
        should contain a list of strings giving the absolute name of each
        component that should be loaded.

        Alternatively, a component may be specified by a dictionary, in which
        case each key of the dictionary is treated as a component name and the
        values as arguments to be passed to the component constructor.

        Any other arguments are passed to the constructor of the pipeline.

        Consider the following ``pipeline_dict``:

        .. highlight:: python

            {
                "checkpoint": "checkpoint_file.chk",
                "start": {"abc": 123, "xyz":789},
                "components": [
                    "phyre_engine.component.dummy.Foo",
                    "phyre_engine.component.dummy.Bar", {
                        "phyre_engine.component.dummy.Baz": [
                            "arg1", "arg2", {
                                "named_arg1": "value1",
                                "named_arg2": "value2",
                            }
                        ]
                    },
                    "phyre_engine.component.dummy.Qux",
                ]
            }

        This will load a pipeline containing the components ``Foo``, ``Bar``,
        ``Baz`` and ``Qux``, each from the ``phyre_engine.component.dummy``
        package. The pipeline will be initialised with the argument
        ``checkpoint="checkpoint_file.chk"`` and ``start={...}``.

        The ``Foo``, ``Bar`` and ``Qux`` components will be intitialised with
        a default empty constructor; the ``Baz`` component will be initialised
        like so:

        .. highlight:: python

            phyre_engine.component.dummy.Baz(
                "arg1", "arg2",
                named_arg1="value1", named_arg2="value2")

        If the :py:attr:`phyre_engine.component.Component.CONFIG_SECTION` class
        variable of a component is not ``None``, then the corresponding section
        of the pipeline configuration is passed into the constructor of the
        component. Explicitly-supplied arguments will override the arguments
        from the configuration.

        To reduce boilerplate, components may be loaded with a shorthand
        notation by supplying a top-level ``namespace`` key and specifying
        shortened component names with a leading dot (``.``). The pipeline
        given above could be specified as:

        .. code-block:: python

            {
                # ...
                "namespace": "phyre_engine.component",
                "components": [
                    ".dummy.Foo",
                    # ...
                ]
                # ...
            }

        The default namespace is ``phyre_engine.component``.

        .. warning::

            Components specified as dictionaries should be specified with *one*
            component per dictionary. Dictionaries do not preserve ordering, so
            passing multiple components in a dictionary can easily break your
            pipeline.

        """
        config = pipeline_dict.get("config", {})
        namespace = pipeline_dict.pop("namespace", DEFAULT_NAMESPACE)
        component_descriptions = pipeline_dict.pop("components", [])
        components = []
        for description in component_descriptions:
            if isinstance(description, dict):
                for cls_name, arg_list in description.items():
                    component = cls._load_component(cls_name, arg_list, config,
                                                    namespace)
            else:
                component = cls._load_component(description, None, config,
                                                namespace)
            components.append(component)
        return cls(components=components, **pipeline_dict)



    class ValidationError(Exception):
        """Raised when a pipeline is found to be invalid.

        The optional argument `data` can be set if this exception is raised at
        runtime rather than at validation. It should contain the current state
        of the system at the point of failure.

        :ivar component: The component that requires a missing key.
        :vartype component: `phyre_engine.component.Component`
        :ivar missing: List of the missing keys.
        :vartype missing: List of strings indicatig the missing keys.
        :ivar data: List of the missing keys.
        :vartype data: Key-value blob describing the state of the pipeline at
            the time this exception was thrown.
        """

        def __init__(self, component, missing, data = None):
            """Returns a new exception indicating an error in component."""
            err_msg = "Component {} was missing keys {}"
            super().__init__(err_msg.format(type(component).__name__, missing))

            self.component = component
            self.missing = missing
            self.data = data if data else {}
