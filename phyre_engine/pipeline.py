"""Provides classes for building, checking and executing pipelines."""

import copy

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
    """

    def __init__(self, components=None, start=None):
        """Initialise a new pipeline with an optional list of components."""
        self.components = components if components else []
        self.start = start if start else {}

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

    def run(self):
        """Run this pipeline, executing each component in turn.

        :returns: Modified data dictionary, with results added by components.

        :raises Pipeline.ValidationError: If a component does not fulfil its
            promises and add a key that the next component requires.
        """

        blob = copy.copy(self.start)
        for cmpt in self.components:
            self.validate_runtime(blob, cmpt)
            blob = cmpt.run(blob)
        return blob


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
