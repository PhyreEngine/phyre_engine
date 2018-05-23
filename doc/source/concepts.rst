.. _phyreengine-concepts:

====================
PhyreEngine concepts
====================

PhyreEngine is a tool for building bioinformatics applications. In most cases,
a bioinformatics application will run a series of tools one after another, with
steps in between each tool to convert between file formats, parse files, or fix
quirks. The application then usually terminates by writing output somewhere.

PhyreEngine treats pipelines as a set of linked *components* that each receive,
modify and output the *pipeline state*. For example, let's consider a trivial
fold recognition pipeline. We'll discuss the pipeline in general terms first
as a toy example, and then see how to turn it into a pipeline that PhyreEngine
can actually process.

A toy fold recognition pipeline
=================================

A fold recognition pipeline needs to read in a sequence, build an "evolutionary
fingerprint" of that sequence to allow remote homology detection, and then scan
that fingerprint against a database of known structures. We can write that down
as some simple steps:

1. **Read sequence**. Open a file containing the query sequence, and read the
   sequence.

2. **Build fingerprint**. Iteratively search the query sequence against a
   sequence database, building a hidden Markov model (HMM) describing sequence
   conservation at each position.

3. **Search known structures**. Align the query HMM with HMMs generated for a
   list of known protein structures. The alignment of HMMs provides very
   sensitve homolgoy detection, so this is the *fold recognition* stage.

4. **Read output of search tool**. The search tool generates results in a
   report file that must now be parsed into a list of hits, with associated
   hit data such as confidence scores, start and stop residues, and hit names.

5. **Write results**. Write the results of the search in some format that can
   easily be consumed for future analysis.

We've added steps 4 and 5 here because there is no point in going to all this
effort just to ignore the results. This is a common theme to PhyreEngine
pipelines: input and output are usually *explicit* rather than implicit, and
require their own steps in the pipeline. This is a deliberate decision that may
seem frustrating to begin with, but will (hopefully) pay dividends in
flexibility.  It can be very useful, especially when testing, to be able to
swap input and output steps in and out of the pipeline.

Each of the steps above represents a *component* of the pipeline, so-called
because we can slot them together like building blocks to form a useful tool.
The collection of data that is sent between components is the *pipeline state*.
Components can modify the pipeline state, but cannot *be* modified by the
pipeline state. That is to say, running a component twice with the same state
will produce the same results (unless a component specifically incoporates
randomness, such as for a molecular dynamics simulation).

.. note::

    A component will produce the same result each time it is called (again,
    except for components specifically designed to be nondeterministic), but it
    will not necessarily *do* the same operations. Many of the components in
    PhyreEngine will resuse results if possible. For example, step 2 above
    first checks to see whether the file containing the HMM is already present
    and reuse it unless specifically told to overwrite it. This is extremely
    useful, as it vastly speeds up turnaround when building pipelines.


.. _creating_pipeline:

Creating a PhyreEngine pipeline
-------------------------------

A PhyreEngine pipeline is specified as a `YAML <http://yaml.org>`_ file. YAML
stands for *YAML Aint Markup Language*. From the YAML website, "YAML is a human
friendly data serialization standard for all programming languages". For our
purposes, we can mostly think of YAML as a nice, human- and machine-readable
format for specifying mappings.

Let's jump straight in and put the pipeline we described above into a format
that can be understood by PhyreStorm. Then we'll go over the syntax of the
pipeline definition and how the individual components are configured.


.. code-block:: yaml
   :linenos:

    pipeline:
      components:
      - .input.ReadSingleSequence
      - .hhsuite.HHBlits:
          database: /data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04
          input_type: sequence
          bin_dir: /data/phyreenginedev/conda/envs/phyreengine
          options:
            output: hhblits_build.hhr
            iterations: 2
            oa3m: query.a3m
            cpu: 10
            verbose: 2
      - .hhsuite.HHSearch:
          database: /bmm/phyreengine/foldlib_mk2/foldlib
          bin_dir: /data/phyreenginedev/conda/envs/phyreengine
          input_type: a3m
          options:
            output: hhsearch_search.hhr
            alt: 1
      - .hhsuite.ReportParser
      - .dump.Csv:
          jmespath_expr: 'templates[].{ template: name, prob: prob }'
          file: folds.csv

This is a lot to take in at first, but it's fairly simple when we break it
down. Line 1 starts the pipeline definition: anything indented relative to the
``pipeline:`` declaration is in a *mapping* named ``pipeline``. It is perfectly
valid to define other top-level mappings: they will just be ignored (unless
speficially referred to by a YAML anchor, but forget that for the moment).

Line 2 starts the ``components`` list. You can tell that is is a list rather
than a mapping because each nested entry---indented relative to the
``components:`` declaration---begins with a hyphen (``-``).

.. note::

    List indentation is a quirk of YAML. The hyphens indicating list elements
    are ignored when determining the indentation of elements under the list.
    For example, the following YAML code blocks both define a scalar ``A``,
    list ``B`` with elements ``1``, ``2``, and ``3``, and another separate
    scalar ``C``:

    .. code-block:: yaml

        A: foo
        B:
          - 1
          - 2
          - 3
        C: bar

    .. code-block:: yaml

        A: foo
        B:
        - 1
        - 2
        - 3
        C: bar

    In this documentation, we prefer the second form, simply because it saves
    some horizontal space.

Each component in the pipeline---``.input.ReadSingleSequence``,
``.hhsuite.HHBlits``, ``.hhsuite.HHSearch``, ``.hhsuite.TabularParser`` and
``.dump.Csv``---is configured with zero or more *parameters*. The
first component, ``.input.ReadSingleSequence``, takes zero parameters.
The second component, ``.hhsuite.HHBlits``, takes four.

.. note::

    The eagle-eyed or those already familiar with YAML might have spotted that
    line 4 ends with a semicolon, but line 3 does not. Internally, this is
    because the first list item is simply a scalar string, which is all that is
    required to configure a component that takes no parameters. The second list
    item is actually a nested mapping, containing the key ``.hhsuite.HHBlits``,
    pointing to the another map containing the options.

    In JSON, the first two elements of the list look like this (cutting long
    lines for brevity):

    .. code-block:: json

        [
          ".input.ReadSingleSequence",
          {
            ".hhsuite.HHBlits":
            {
              "database": "...",
              "input_type": "sequence",
              "bin_dir": "...",
              "options": {
                "output": "hhblits_build.hhr",
                "iterations": 2,
                "oa3m": "query.a3m",
                "cpu": 10,
                "verbose": 2
              }
            }
          }
        ]

    You can see why YAML was preferred over JSON for this task. Interestingly,
    JSON is actually valid YAML, so the masochists among us could choose to use
    JSON-formatted pipelines.

    The end result of this is that if a component takes no parameters, specify
    its name without a trailing colon. If you need to pass parameters, you need
    to include the colon.

A description of the parameters available for each component should be given on
the corresponding component documentation page. See, for example,
:py:class:`~phyre_engine.component.hhsuite.HHBlits` for the parameters taken by
``.hhsuite.HHBlits``. If the parameters for a component are not adequately
explained in the documentation, *this is a bug*! Please let us know.

This pipeline corresonds exactly to the steps we outlined above. Let's go
through each component in detail, examining the items that are added to the
pipeline state at each step.

1. ``.input.ReadSingleSequence``: This component reads a single sequence from
   a file in any format that can be understood by BioPython. It requires the
   ``input`` field in the pipeline state to contain the filename, which we
   set on the command line when calling PhyreEngine. The following elements
   are added to the pipeline state:

   ``sequence``
      The sequence that is encoded in the file. This is a simple string of
      single-letter amino acid codes.

   ``id``
      An ID of this sequence, parsed from the sequence file. For a FASTA file,
      the ID of the sequence is the first word in the sequence description
      line.

   ``name``
      Similarly, the sequence name parsed from the sequence file. This is
      identical to the ``id`` for a FASTA file.

   ``description``
      The sequence description, again parsed from the sequence file. This is
      the entire description line for a FASTA file.

2. ``.hhsuite.HHBlits``: Uses a local installation of
   `hhblits <https://toolkit.tuebingen.mpg.de/#/tools/hhblits>`_ to iteratively
   align the query sequence to sequence database, producing a multiple
   sequence alignment (MSA) that can be be converted to an HMM and aligned with
   a fold library. This component adds the following fields to the pipeline
   state:

   ``report``
      Name of the report file produced by hhblits. This is ignored by the
      remainder of the pipeline, but can be useful for manual checks. The
      name of the report file is configured by the ``output`` option.

   ``a3m``
      The file name of the multiple sequence alignment generated by hhblits.
      This is controlled by the ``oa3m`` option, which corresponds exactly
      with the ``-oa3m`` command line switch of the ``hbblits`` executable.
      If this option is not given, the ``a3m`` key is not generated and added
      to the pipeline state, which will cause an error when the next component
      is called.

3. ``.hhsuite.HHSearch:``: Uses a local installation of hhsearch to align the
   MSA produced in step 2 with a fold library. The only key added to the
   pipeline state is ``report``, which will overwrite the field set in the
   previous step.

4. ``.hhsuite.ReportParser``: Parse the report file generated by hhsearch into
   a list of hits with associated scores and alignment ranges. Note that this
   does *not* read the alignments from the report file, as the report file is
   written mainly for use by humans and is not easily parsed. To read
   alignments, pass the ``atab`` option in the previous step, and use a
   :py:class:`~phyre_engine.component.hhsuite.TabularParser` after the
   ``ReportParser``. This component adds the ``templates`` list to the
   pipeline state. Each element of the ``templates`` list will contain the
   following fields:

   ``name``
      Name of the hit.

   ``prob``
      The confidence (i.e. likelihood of homology) of this hit.

   Several more fields are added, each containing a different score. See
   :py:class:`phyre_engine.tools.hhsuite.parser.Report.Hit` for a full list.

5. ``.dump.Csv``: Dump the pipeline state to a CSV file. The format of the
   CSV file is controlled by the parameter ``jmespath_expr``, which is a
   `JMESPath <http://jmespath.org/>`_ query. I won't cover JMESPath in
   detail here---`the tutorial <http://jmespath.org/tutorial.html>`_ is quite
   good---but this expression loops over each element of the ``templates``
   list and selects the ``name`` field (renamed as ``template``) and ``prob``
   field. The ``.dump.Csv`` component then prints each field as a column.

The output of the ``.dump.Csv`` component is stored in ``folds.csv``
(configured via the ``file`` parameter), and contains a list of hits:

.. code-block:: none

    prob,template
    100.0,5NWK_B
    100.0,5NWJ_A
    100.0,5NWK_A
    100.0,5NWK_F
    100.0,2br9_A
    ...

Running pipelines
-----------------

To run a pipeline, use the command :command:`phyre_engine`. Alternatively, you
can run the module explicitly with ``python -mphyre_engine.run``, which is
synonymous with ``phyre_engine``: this can be useful when you are using
``conda develop`` and scripts have not been installed.

Let's try running our fold recognition pipeline. Simply pass the name of the
file containing the pipeline to :command:`phyre_engine`. Do this in an empty
directory, because the files generated by the pipeline will be placed in the
current directory.

.. code-block:: none

    phyre_engine ../fold_recognition.yml

Oops: that will give a nasty-looking error message. Python error message are a
little strange, and should be read bottom-up. At the bottom, we can see the
error: ``phyre_engine.pipeline.ValidationError: Component ReadSingleSequence
was missing keys ['input']``. This is saying that validation of the pipeline
state failed, because the first component was missing the ``input`` field of
the pipeline state, which makes sense because at no point have we actually told
the pipeline where our query is.

This is where the ``--start`` option to :command:`phyre_engine` comes in. It
allows us to specify the starting pipeline state separating key names and
values with a colon. Assuming your query sequence is stored in ``query.fasta``,
you would run:

.. code-block:: none

    phyre_engine --start input:query.fasta ../fold_recognition.yml

If the locations of your databases and software are configured correctly in the
pipeline defition, the pipeline will run to completion and evenetually produce
``folds.csv``. More likely, it will fail because ``hblits``, ``hhsearch`` and
their supporting databases are not installed in the correct location.


.. _what_is_a_component:

.. index::
    single: component; definition

What is a component?
====================

In the previous section, "components" were treated in a fairly abstract way.
It's obvious that at some point components must actually execute some code. In
pthis section, we discuss how to define a component, and how components are
looked up when they are listed in a pipeline YAML file.

A component is simply a Python class that implements the interface defined by
:py:class:`~phyre_engine.component.component.Component`. Basically, this means
a Python class with a ``run`` method that takes the pipeline state as an
argument and returns the modified state. We can define a trivial component
easily:

.. code-block:: python

    from phyre_engine.component import Component
    class Multiply(Component):
        """
        Multiply a field in the pipeline state by a constant factor.

        :param str var: Name of the field to multiply.
        :param num by: Constant factor.
        """
        ADDS = []
        REMOVES = []

        @property
        def REQUIRED(self):
            return [self.var]

        def __init__(self, var, by):
            self.var = var
            self.by = by

        def run(self, data, config=None, pipeline=None):
            """Multiply field by constant factor."""
            data[self.var] = data[self.var] * self.by
            return data

This component accepts two parameters: ``var``, which is the name of a field in
the pipeline state, and ``by``, which is a constant factor by which ``var``
will be multiplied. In the ``run`` parameter, you can see the logic. The
pipeline state is manipulated like any other Python dictionary (because it is
like any other dictionary). Ignore the ``config`` and ``pipeline`` parameters
for now: those can sometimes be useful for defining components that call child
pipelines.

We have included a class docstring describing the function of the component.
Docstrings are parsed by Sphinx when generating documentation, and serve as the
main reference for each component. The idea is to keep documentation and code
close, so that there is at least a chance that the documentation actually
matches the reality. The docstring for the ``run`` method is trivial, because
the purpose of the class should already have been described in the class
docstring.

The only other interesting thing about this class are the ``ADDS``, ``REMOVES``
and ``REQUIRED`` properties. These are used by the pipeline when calling
components to check whether the pipeline state is missing any vital fields. As
you would expect, the ``ADDS`` property is the list of items added to the
pipeline state, ``REMOVES`` are removed, and ``REQUIRED`` required.

In most cases, these will be constant arrays like the empty ``ADDS`` and
``REMOVES`` list, which say that this component does not add or remove any
fields. In some cases, illustrated here, the lists of items might not be known
until the component is configured. It is possible to work around this by just
declaring ``REQUIRED = []``, which will allow any pipeline state to be passed
to the component. A better approach, shown here, is to define a Python
`property <https://docs.python.org/3/library/functions.html#property>`_ that
can be treated like a constant array (i.e. ``Multiply("a", 1).REQUIRED`` is
valid) but is determined at runtime.

.. warning::

    A common mistake is to forget to return the modified pipeline state at the
    end of the ``run`` method. In general, ``run`` should always end with
    ``return data``.


.. index::
    single: component; configuration
    single: pipeline; configuration
    single: pipeline; configuration; namespace

Components in YAML files
========================

So far, we know that components are just Python classes, and that pipelines are
defined as YAML files with a list of components, but it is not immediately
clear how a component is configured and executed based on the information in
the pipeline file. In this section, we will discuss how a component is looked
up and configured.

.. index::
    single: component; lookup
    single: pipeline; namespace

Component names
---------------

Component classes are simply looked up by name, with an optional "namespace"
prepended if the name starts with a ``.``. The default namespace is
``phyre_engine.component``. For example, in :ref:`creating_pipeline`, we used
the component ``.input.ReadSingleSequence``. The component name begins with a
dot, so it is actually
:py:class:`phyre_engine.component.input.ReadSingleSequence`, which is then
loaded using Python's usual :py:mod:`importlib` machinery. In this case, the
:py:class:`~phyre_engine.component.input.ReadSingleSequence` class is loaded
from the :py:mod:`phyre_engine.component.input` module.

.. container:: toggle

    .. container:: header

        **Show/Hide Advanced**

    Python allows nested classes. None of the components distributed with
    PhyreEngine are nested classes (except in some tests), but it *is* possible
    to specify nested classes in a YAML file by specifying a class name as a
    two-part list containing the module name and class name.

    For example, specifying a component as ``[mod, Parent.Child]`` will load
    the ``Child`` inner class from the ``Parent`` outer class in the ``mod``
    module. This is facilitated by some extra YAML-parsing machinery in
    :py:mod:`phyre_engine.tools.yaml` to parse lists of scalars as tuples,
    which is required so nested class specifications can be used as a map key.

If you are developing a large library of custom modules and only using
PhyreEngine as a generic way of running a pipeline, you may wish to change the
prefix that is prepended to dotted component names. You can do this by setting
the ``namespace`` field in the ``pipeline`` mapping:

.. code-block:: yaml

    pipeline:
      namespace: foo
      components:
      - .bar.A       # Loads foo.bar.A
      - full.path.B  # No leading dot, so loads full.path.B as usual.

.. _component_parameters:

Component parameters
--------------------

In :ref:`creating_pipeline`, we used the
:py:class:`phyre_engine.component.hhsuite.HHBlits` class to build a profile
from our sequence. The :py:class:`~phyre_engine.component.hhsuite.HHBlits`
component accepts several configuration parameters: the location of the profile
database, the number of CPUs to use, and so on. These were passed to the
component by supplying them in the YAML file as the elements of a mapping:

.. code-block:: yaml

    .hhsuite.HHBlits:
      database: /data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04
      input_type: sequence
      bin_dir: /data/phyreenginedev/conda/envs/phyreengine
      options:
        output: hhblits_build.hhr
        iterations: 2
        oa3m: query.a3m
        cpu: 10
        verbose: 2

We know that components are simply Python classes that are looked up by the
names given in the YAML file, so it is probably not surprising that the
parameters specified in the YAML file are passed directly to the class
constructor. The previous block of YAML code could be expressed in Python like
this:

.. code-block:: python

    import phyre_engine.component.hhsuite
    phyre_engine.component.hhsuite.HHBlits(
      database="/data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04",
      input_type="sequence",
      bin_dir="/data/phyreenginedev/conda/envs/phyreengine",
      options={
        "output": "hhblits_build.hhr",
        "iterations": 2,
        "oa3m": "query.a3m",
        "cpu": 10,
        "verbose": 2,
      })


.. warning::

    Knowing that arguments in the YAML file are passed directly to the class
    constructor, you might think that it should be possible to pass positional
    arguments by specifying an array:

    .. code-block:: yaml

        # This will NOT call .foo.Bar(1, 2, 3).
        # Positional arguments are not accepted.
        .foo.Bar: [1, 2, 3]

    This is not allowed, as it leads to (even more) dangerously obtuse
    pipelines. It is better to be explicit at the cost of a few keystrokes.


.. index::
    single: CONFIG_SECTION
    single: pipeline; configuration

.. _pipeline_config:

Pipeline configuration
======================

Passing parameters to individual components as described in
:ref:`component_parameters` can quickly become quite verbose. In many cases,
components will share related configuration parameters. For example, it is
likely that all of the programs used by the
:py:mod:`phyre_engine.component.hhsuite` module will share the same ``bin_dir``
parameter.

PhyreEngine allows for a pipeline-wide configuration to be defined. Components
may then use keys from this component as an initial configuration. For example,
the components in the :py:mod:`~phyre_engine.component.hhsuite` module take the
keys in the ``hhsuite`` section of the pipeline configuration. We can specify
parameters for the ``hhsuite`` classes like this:

.. code-block:: yaml

    pipeline:
      config:
        hhsuite:
          bin_dir: /data/phyreenginedev/conda/envs/phyreengine
          options:
            cpu: 10
            verbose: 2

      components:
      # ...
      - .hhsuite.HHBlits:
          database: /data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04
          input_type: sequence
          options:
            output: hhblits_build.hhr
            iterations: 2
            oa3m: query.a3m
      - .hhsuite.HHSearch:
          database: /bmm/phyreengine/foldlib_mk2/foldlib
          input_type: a3m
          options:
            output: hhsearch_search.hhr
            alt: 1
      # ...

The :py:class:`~phyre_engine.component.hhsuite.HHBlits` and
:py:class:`~phyre_engine.component.hhsuite.HHSearch` components now take the
``bin_dir`` parameter from the pipeline configuration, along with the
``options`` map.

There are two important things to note here:

* Parameters specified in the ``components`` section will always override
  parameters set in the ``config`` section. This allows for a general
  configuration to be set, and tweaked when necessary.

* Configuration parameters are "deep-merged": each map in the configuration is
  merged, rather than being overriden. In the example above, the ``options``
  parameter is set in both the ``config`` and for each component. Each
  ``options`` map is merged separately, so HHBlits sees the ``output``,
  ``iterations`` and ``oa3m`` options along with the ``cpu`` and ``verbose``
  options from the configuration. Likewise, HHSearch sees the ``output`` and
  ``alt`` options along with the ``cpu`` and ``verbose`` options.

Internally, components are configured by the
:py:class:`phyre_engine.pipeline.Pipeline` class. When a pipeline file is
loaded and the components are being created, the ``config`` method of that
component is called with the pipeline configuration and component arguments as
parameters. The ``config`` method must then merge the two and return the
arguments used to initialise the class.

The default :py:meth:`~phyre_engine.component.component.Component.config`
method first checks whether the component has defined a ``CONFIG_SECTION`` key.
If it has not, then the component parameters are returned without modification
and the pipeline configuration discarded for that component. If a
``CONFIG_SECTION`` key *is* defined, then the corresponding configuration
section is extracted from the pipeline configuration and merged with the
component parameters. In the case of the hhsuite components, each component
defines ``CONFIG_SECTION = "hhsuite"``. More complex ``config`` methods may be
defined (see :py:meth:`phyre_engine.component.foldlib.BuildProfiles.config`,
for example), but this is generally discouraged as it promotes spooky action at
a distance.

.. _global_config:

Global configuration
--------------------

Finally, there is one last source of configuration data. If a global
configuration file is found, then it is used as the first source of
configuration data, before being overwritten by a per-pipeline configuration.
The location of this global configuration file will depend on the platform: on
Linux, it is likely to be :file:`/home/{user}/.config/phyreengine/config`. The
exact path can be seen by executing the following in a python interpreter
(assuming that the PhyreEngine dependencies are installed correctly)::

    >>> import appdirs
    >>> print(appdirs.user_config_dir("phyreengine") + "/config")'

Summary
-------

To summarise, a component is configured like so:

1. PhyreEngine will read the global configuration file, if it exists.

2. PhyreEngine reads the pipeline file, and looks for a ``config`` section
   in the ``pipeline`` section. If a global configuration was found, it is
   merged with this per-pipeline configuration. The pipeline configuration
   takes precedence. Note that because the global configuration *only*
   supplies configuration data, there is no need to wrap it in a
   ``pipeline:`` and ``config:`` map.

3. For each component in the ``components`` section of the pipeline,
   PhyreEngine will:

   a. Parse the component parameters from the pipeline file.

   b. Call the component's ``config`` method, passing the entire pipeline
      configuration and the component parameters. The ``config`` method
      will return a dictionary of parameters by merging the pipeline
      configuration and the component parameters.

      The default behaviour of the ``config`` method is to look for a
      ``CONFIG_SECTION`` class attribute in the component.

      - If a ``CONFIG_SECTION`` attribute is found, the corresponding
        section is looked up in the pipeline configuration and merged with
        the component parameters if it is present.

      - If no ``CONFIG_SECTION`` attribute is present, the component parameters
        are returned without being touched.

    c. The dictionary returned from the ``config`` method is passed into the
       component's ``__init__`` method.
