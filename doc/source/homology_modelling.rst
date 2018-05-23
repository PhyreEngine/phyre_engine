======================================
A detailed homology modelling pipeline
======================================

In this guide, we will go through the process of running a full-featured
homology modelling pipeline. We will describe in detail how to run PhyreEngine,
how the pipeline and each component is configured, and then step through each
component examining how it affects the pipeline state.

Here is the pipeline that we will be examining. For long sections of code like
this, the full contents are shown collapsed by default. When we discuss smaller
sections in detail, we will include the expanded excerpts.

.. container:: toggle

    .. container:: header

        **Show/Hide pipeline**

    .. literalinclude:: /pipelines/homology_model.cluster-serial.yml
        :language: yaml
        :linenos:

As homology modelling pipelines go, this one is quite complex. It is not a toy
example, and is fully suited to running on a homology modelling web server.

Combined with this pipeline, we have the following :ref:`global configuration
<global_config>`, saved in :file:`~/.config/phyreengine/config.yml`.

.. container:: toggle

    .. container:: header

        **Show/Hide config**

    .. literalinclude:: /pipelines/config.yml.sample
        :language: yaml
        :linenos:

Files
-----

The pipeline listed above is saved in the current directory as
:file:`homology-model.yml`. In the same directory, we have our query sequence in
:file:`query.fasta`:

.. code-block:: console

    $ tree
    .
    ├── homology-model.yml
    └── query.fasta
    $ cat query.fasta
    >Q96NU7
    MASGHSLLLENAQQVVLVCARGERFLARDALRSLAVLEGASLVVGKDGFIKAIGPADVIQRQFSGETFEEIIDCSGKCIL
    PGLVDAHTHPVWAGERVHEFAMKLAGATYMEIHQAGGGIHFTVERTRQATEEELFRSLQQRLQCMMRAGTTLVECKSGYG
    LDLETELKMLRVIERARRELDIGISATYCGAHSVPKGKTATEAADDIINNHLPKLKELGRNGEIHVDNIDVFCEKGVFDL
    DSTRRILQRGKDIGLQINFHGDELHPMKAAELGAELGAQAISHLEEVSDEGIVAMATARCSAILLPTTAYMLRLKQPRAR
    KMLDEGVIVALGSDFNPNAYCFSMPMVMHLACVNMRMSMPEALAAATINAAYALGKSHTHGSLEVGKQGDLIIINSSRWE
    HLIYQFGGHHELIEYVIAKGKLIYKT

The query sequence arbitrarily chosen to be the UniProt entry `Q96NU7
<http://www.uniprot.org/uniprot/Q96NU7>`_.


.. index::
    single: logging; format

Running the pipeline
--------------------

I will assume that PhyreEngine has been installed as a :ref:`development package
<developing-phyreengine>` using :ref:`conda <installing-via-conda>`, and that
the environment in which PhyreEngine is installed has `been activated
<https://conda.io/docs/user-guide/tasks/manage-environments.html#activating-an-environment>`_.
When PhyreEngine is installed in this way, the Python modules that make up the
core of PhyreEngine are available to Python, but the :command:`phyre_engine`
script will not have been installed. We will work around this by calling the
:py:mod:`phyre_engine.run` module directly:

.. code-block:: console

    $ python -mphyre_engine.run
    usage: run.py [-h] [-v] [-e] [-s START] [-c CONFIG] pipeline
    run.py: error: the following arguments are required: pipeline

From the brief summary that is given when :py:mod:`phyre_engine.run` is called
without any arguments, we can see that PhyreEngine is expecting a pipeline.
Let's oblige it:

.. code-block:: console

    $ python -mphyre_engine.run homology-model.yml
    WARNING  : 2018-05-16 16:28:40,125 : phyreengine : phyre_engine.run : Some components might be missing input. This analysis is only an educated guess, so execution is continuing. Silence this warning with --no-static-check or the no_static_check config variable.
    WARNING  : 2018-05-16 16:28:40,125 : phyreengine : phyre_engine.run : Component phyre_engine.component.util.ChangeDir might be missing keys ['working_dir']
    WARNING  : 2018-05-16 16:28:40,125 : phyreengine : phyre_engine.run : Component phyre_engine.component.input.ReadSingleSequence might be missing keys ['input']
    INFO     : 2018-05-16 16:28:40,125 : phyreengine : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.util.ChangeDir
    ERROR    : 2018-05-16 16:28:40,125 : phyreengine : phyre_engine.run : Component phyre_engine.component.util.ChangeDir expected the missing keys ['working_dir'] in the pipeline state

That doesn't look good. These logging messages are telling that the pipeline is
missing some keys in the :ref:`pipeline state <phyreengine-concepts>`, but
before we examine that, it's worth looking taking a quick look at how
PhyreEngine logs errors or other information.

.. _logging:

Aside: Logging
~~~~~~~~~~~~~~

.. container:: aside

    PhyreEngine uses Python's :py:mod:`logging` module for logging. If you're
    familiar with Python's logging system, you may have recognised the
    ``logging`` section of the configuration file we included above:

    .. literalinclude:: /pipelines/config.yml.sample
        :start-after: logging:
        :prepend: logging:
        :language: yaml

    This mapping is passed directly to :py:func:`logging.config.dictConfig` in
    order to configure PhyreEngine's logging. If you've not used Python's
    logging functionality before, the logging configuration given above is a
    sensible default for development work but probably too verbose for general
    use. Change the line that says ``level: DEBUG`` to ``level: INFO`` or
    ``level: WARNING`` to quiet PhyreEngine down a little.

    In the console output shown above, we saw the following message:

    .. code-block:: none

        ERROR    : 2018-05-16 16:44:27,025 : phyreengine : phyre_engine.run : Component phyre_engine.component.util.ChangeDir expected the missing keys ['working_dir'] in the pipeline state

    This is split into five fields separated by colons and corresponding to the
    ``format`` string in the ``logging`` configuration.

    1. ``ERROR``: The severity of the logging event. ``ERROR`` indicates that
       the pipeline has failed. A ``CRITICAL`` error is shown in the event of an
       error so catastrophic the system has become unusable. The remaining
       categories, ``WARNING``, ``INFO`` and ``DEBUG`` provide progressively
       more fine-grained information about PhyreEngine's internal state.

    2. ``2018-05-16 16:44:27,025``: The date and time, in The One True Date
       Format (``YYYY-MM-DD``).

    3. ``phyreengine``: The hostname of the current system. In this case, the
       system is named ``phyreengine``, just to add some extra confusion.

    4. ``phyre_engine.run``: The name of the module or class emitting the
       message.  In this case, errors are being emitted directly from the
       :py:mod:`phyre_engine.run` module.

    5. ``Component phyre_engine....``: A free-text message describing the event
       that caused the logging message to be emitted.

    .. note::

        For the sake of horizontal space, I will omit the date and hostname from
        any future logging messages included in this document.

Returning from that aside about logging, let's return to examing the output of
PhyreEngine when we executed our homology modelling pipeline. It tells us that
``ChangeDir might be missing keys ['working_dir']`` and that
``ReadSingleSequence might be missing keys ['input']``. These warnings let us
know that the pipeline looks like it is missing some data in the pipeline state.
Sure enough, as soon PhyreEngine begins to execute the first component, we
see a fatal error telling us that the the ``working_dir`` key was, in fact,
required: ``ChangeDir expected the missing keys ['working_dir'] in the pipeline
state``.

I will discuss the :ref:`example-util-changedir` and
:ref:`example-util-changedir` components and their expected inputs later in this
document, but for the moment it is sufficient to know that these error messages
tell us that the pipeline state is missing some important data: a working
directory and the name of a sequence file. This is deliberate, because there is
no point in a homology modelling pipeline that can only operate on sequences in
the current directory with a fixed filename (I've seen lots of bioinformatics
tools that do this and I hate it).

We can kick-start the pipeline state with the :option:`--start` command-line
option, which takes a colon-separated key-value pair. The :option`--start`
option can be supplied as may times as necessary. Let's tell PhyreEngine to run
in the current directory, and process the file :file:`query.fasta`:

.. code-block:: console

    $ python -mphyre_engine.run \
    >    --start working_dir:. \
    >    --start input:query.fasta \
    >    homology-model.yml

By specifying the two :option:`--start` options here, the pipeline state
starts with the following contents:

.. literalinclude:: /homology_modelling_session/dump-00.pp
    :caption: Initial pipeline state
    :language: python

I will now go through each component in the pipeline. At the beginning of each
section, I will give a pretty-printed representation of any changes made to the
pipeline state by the component, and the console contents printed before
immediately before the component is run.

.. _example-util-changedir:

Component :py:class:`.util.ChangeDir`
-------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .util.ChangeDir

.. literalinclude:: /homology_modelling_session/dump-01.pp
    :diff: /homology_modelling_session/dump-00.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.util.ChangeDir
    INFO     : phyre_engine.component.util.ChangeDir : Moving to directory '.'.

The :py:class:`.util.ChangeDir <phyre_engine.component.util.ChangeDir>` is used
to change the working directory of the PhyreEngine process. Most operations will
write and read files from the working directory unless an absolute path is
explicitly specified. This component is often paired with
:py:class:`.util.MakeDir <phyre_engine.component.util.MakeDir>` to generate
directories based on, for example, the name of a sequence. In this case the
pipeline is intended to be called from a web server that initialises certain
metadata files in a working directory, so :py:class:`~.ChangeDir` should—and
will, by raising a :py:exc:`FileNotFoundError` exception—cause the pipeline to
fail.

This component reads the name of the directory into which it will change from
the ``working_dir`` key of the pipeline state, which is why we needed to set it
on the command line with the :option:`--start` option.

From the listing of the pipeline state given above, we can see that the
``directory_stack`` key was added to the pipeline state. It has been truncated
here for readability, but this key contains the previous working directory. If
the pipeline were to run :py:class:`~.ChangeDir` again, a new item would be
pushed onto the directory stack. The component :py:class:`:.util.PopDir
<phyre_engine.component.util.PopDir>` can be used to pop directories off of the
directory stack and walk back through directories.

.. note::

    When running components wrapped in a :py:class:`~.TryCatch` component, is
    possible for a section of a pipeline to fail but for execution of the
    pipeline as a whole to continue. If this is the case, wrapping groups of
    components in :py:class:`~.ChangeDir` and :py:class:`~.PopDir` components is
    not sufficient: if the pipeline fails before the :py:class:`~.PopDir` is
    executed, then :py:class:`~.ChangeDir` will be executed from the *new*
    working directory. Instead, use the :py:class:`~.SaveRootDir` component to
    set a "root" directory and a :py:class:`~.RestoreRootDir` component to
    restore that root directory.

So, after the :py:class:`~.ChangeDir` component has been executed, PhyreEngine
is running in the directory specified with ``--start working_dir``. For this
example, we are just using the current directory, ``.``.

.. _example-input-readsinglesequence:

Component :py:class:`.input.ReadSingleSequence`
-----------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .input.ReadSingleSequence

.. literalinclude:: /homology_modelling_session/dump-02.pp
    :diff: /homology_modelling_session/dump-01.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 1: phyre_engine.component.input.ReadSingleSequence

The :py:class:`.ReadSingleSequence` component reads a sequence and any
associated metadata file from the file given in the ``input`` section of the
pipeline state. Internally, PhyreEngine uses BioPython's `Bio.SeqIO
<http://biopython.org/wiki/SeqIO>`_ interface to read sequences, so this
component can read most common sequence formats. The format of the file must be
set with the ``format`` parameter of the component. The default format is
``fasta``, so our sequence file is read without any problems. If we had a
`Stockholm
<http://bioperl.org/formats/alignment_formats/Stockholm_multiple_alignment_format>`_
file, we would have configured the component like so:

.. code-block:: yaml

    - .input.ReadSingleSequence:
        format: stockholm

The :py:class:`~.ReadSingleSequence`, and its companion class
:py:class:`~.ReadMultipleSequences` also read the sequence ``id``, ``name`` and
``description`` from the sequence file, and add these keys to the pipeline
state. Our sequence is in FASTA format, which doesn't have well-defined metadata
fields, so BioPython has parsed some sensible values from the sequence
identifier. Metadata parsing can be disabled by setting the ``metadata``
parameter to ``False``.

.. _example-validate-seqvalidator:

Component :py:class:`.validate.SeqValidator`
--------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .validate.SeqValidator

.. code-block:: none
    :caption: Change in pipeline state

    No change

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 2: phyre_engine.component.validate.SeqValidator

Because this pipeline is designed to be run from a web server, it will be
subjected to abuse from people pasting DNA sequences, random text, attempts at
SQL injection, and so on. The
:py:class:`phyre_engine.component.validate.SeqValidator` component will raise an
exception and abort the pipeline if any of the letters in the protein sequence
are non-standard amino acid codes. The parameter ``extended`` could be set to
``True`` to also allow the letters ``X``, ``B`` and ``Z``. This component will
often be paired with :py:class:`phyre_engine.component.validate.SeqLenFilter`
to prevent processing of sequences that are considered too short or too long.

.. _example-hhsuite-hhblits:

Component :py:class:`.hhsuite.HHBlits`
--------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    # Generate HMM and search against the fold library.
    - .hhsuite.HHBlits:
        database: !template '{LOCAL_DATA_PREFIX}/hhsuite/uniclust30_2017_04/uniclust30_2017_04'
        input_type: sequence
        options:
          output: hhblits_build.hhr
          iterations: 2
          oa3m: query.a3m
          alt: 1

.. literalinclude:: /homology_modelling_session/dump-04.pp
    :diff: /homology_modelling_session/dump-03.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output (snipped and pretty-printed)

    INFO     : phyre_engine.pipeline.Pipeline : Running component 3: phyre_engine.component.hhsuite.HHBlits
    DEBUG    : phyre_engine.component.hhsuite.HHBlits : Running [
        '/data/phyreenginedev/conda/envs/phyreengine/bin/hhblits',
        '-v', '2',
        '-cpu', '16',
        '-o', 'hhblits_build.hhr',
        '-n', '2',
        '-oa3m', 'query.a3m',
        '-alt', '1',
        '-d', '/data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04',
        '-i', '/tmp/query-fnwgk8nu.fasta']

    - 19:20:44.126 INFO: Searching 12365503 column state sequences.

    - 19:20:44.184 INFO: /tmp/query-fnwgk8nu.fasta is in A2M, A3M or FASTA format

    - 19:20:44.184 INFO: Iteration 1

    - 19:20:44.268 INFO: Prefiltering database

    <✂✂✂>

    Query         Q96NU7
    Match_columns 426
    No_of_seqs    1688 out of 35249
    Neff          11.6019
    Searched_HMMs 5499
    Date          Wed May 16 19:21:48 2018
    Command       /data/phyreenginedev/conda/envs/phyreengine/bin/hhblits -v 2 -cpu 16 -o hhblits_build.hhr -n 2 -oa3m query.a3m -alt 1 -d /data/phyre3/data/hhsuite/uniclust30_2017_04/uniclust30_2017_04 -i /tmp/query-fnwgk8nu.fasta

     No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
      1 tr|G1SP40|G1SP40_RABIT Unchara 100.0 1.5E-69 5.4E-75  502.5   0.0  426    1-426    66-491 (491)
      2 tr|S7NJX7|S7NJX7_MYOBR Putativ 100.0 8.1E-66   3E-71  476.6   0.0  417    2-425    75-492 (493)
      3 tr|A0A0W1R2Z1|A0A0W1R2Z1_9EURY 100.0   2E-63 7.6E-69  462.7   0.0  404    4-425     1-404 (405)
      4 tr|A0A158RF71|A0A158RF71_HYDTA 100.0 3.2E-61 1.2E-66  428.5   0.0  420    5-426     1-474 (481)
      5 tr|B7RYF4|B7RYF4_9GAMM 5-methy 100.0   6E-60 2.3E-65  447.8   0.0  377    3-426    19-429 (455)

    <✂✂✂>

Now, we reach the meat of the pipeline. The
:py:class:`phyre_engine.component.hhsuite.HHBlits` component runs `HHBlits
<https://www.nature.com/articles/nmeth.1818>`_, part of the `hhsuite
<https://github.com/soedinglab/hh-suite>`_ set of tools, to generate a multiple
sequence alignment (MSA). The MSA will be used in the next step as an
"evolutionary fingerprint" to search the fold library.

This component requires more configuration than we have seen in the previous
components, so it's worth taking a close look at how it is configured. The
configuration of the component in the pipeline file is shown in the :ref:`usage
example <example-hhsuite-hhblits>` shown above.

This is, hopefully, fairly intuitive if you're familiar with hhblits. The
``database`` parameter points to a sequence database in hhsuite format. Here, we
are using the `uniclust30 <https://uniclust.mmseqs.com/>`_ database. The
``input_type`` parameter tells HHBlits that we are going to be passing a
sequence, rather than an MSA or hidden Markov Model, which is necessary as the
component needs to know which parameters are going to be read from the pipeline
state. The default ``input_type`` is ``sequence``, but it is worth specifying
explicitly to avoid confusion. The available options are given in the
:py:class:`~phyre_engine.component.hhsuite.QueryType` enum.

The ``options`` key passed to :py:class:`~.hhsuite.HHBlits` is passed to the
:command:`hhblits` command as command line arguments, with the appropriate
prefixing hyphen added. Some obtuse options have more descriptive aliases
associated with them—``output`` is an alias for ``-o`` and ``iterations`` for
``-n`` — but the single-character versions are perfectly acceptable, if less
readable. The command line used to actually run the :command:`hhblits`
executable is shown in the console output (and wrapped here for readability).

But how does :py:class:`~.hhsuite.HHBlits` know how to expand the
``{LOCAL_DATA_PREFIX}`` field in the ``database`` template? The YAML tag
``!template`` allows us to include parameters from the pipeline configuration
our custom :py:mod:`.tools.yaml <YAML loader>`, but the ``config`` section of
our pipeline is empty, so we think about how :ref:`pipelines are configured
<pipeline_config>` and look at the global configuration file,
:file:`~/.config/phyreengine/config.yml`.

In the configuration file, we note two likely-looking sections:

.. literalinclude:: /pipelines/config.yml.sample
    :end-before: qsub:
    :language: yaml
    :emphasize-lines: 4

.. literalinclude:: /pipelines/config.yml.sample
    :start-after: scwrl4:
    :end-before: By default
    :lines: 2-
    :language: yaml

The first section defines the ``LOCAL_DATA_PREFIX`` field, used in our pipeline
file to define the location of the HHBlits sequence database. The second section
defines the ``hhsuite`` section, which is merged with the parameters given in
our pipeline file. The ``bin_dir`` parameter tells the
:py:class:`.hhsuite.HHBlits` component where to find the :command:`hhblits`
executable, and the ``HHLIB`` parameter tells the component how to set the
:envvar:`HHLIB` environment variable, which is required by all the hh-suite
tools.

.. note::

    All of the components that wrap an external tool will accept a ``bin_dir``
    parameter telling that component how to find the executable. If the
    parameter is not supplied, components will search the system :envvar:`PATH`,
    but the system path relies on the current shell. Importantly, conda
    environments may not be activated if PhyreEngine is run from a non-login
    shell (e.g. from a cron job). It is recommended that the ``bin_dir``
    parameters be supplied. The boilerplate this requires was one of the big
    motivators for having PhyreEngine read a global configuration file!

If we take a look at the keys added to the pipeline state by this component, we
can see that the additions are quite lightweight. The only keys that have been
added are ``a3m`` and ``report``, which respectively contain the file names of
the a3m (MSA) and report file generated by :command:`hhblits`. The
:py:class:`~.hhsuite.HHBlits` component makes no effort to parse these files,
and we wouldn't want it to: we don't particularly care about the alignment of
sequences in the MSA, and the report file is solely for our records. If we
wished to parse the report file, we would include a
:py:class:`.hhsuite.ReportParser` component; we *do* include one of these
components later in this pipeline to parse the output of hhsearch when we scan
the fold library.

Aside: External tools
~~~~~~~~~~~~~~~~~~~~~

.. container:: aside

    Since so many components in a bioinformatics pipeline simply (or
    "complexly", as the case may be) wrap an existing third-party tool,
    PhyreEngine comes with tools that can be used to aid the wrapping of these
    tools. As of May 2018, PhyreEngine contains interfaces for 38 third-party
    tools.

    It is fairly common for the command-line interfaces for popular tools to be
    quite obtuse. Many bioinformatics tools were written when single-character
    command-line options were the norm, and then grown features over the years
    and been forced to retain a clumsy command-line interface for the sake of
    backwards compatibility.

    For example, the `Scwrl4 <http://dunbrack.fccc.edu/scwrl4/>`_ executable
    takes lots of options specified with single characters that are hard to
    remember. PhyreEngine provides the :py:class:`~.ExternalTool` class for
    alisaing these command line flags:

    .. code-block:: python

        SCWRL4 = ExternalTool(flag_map={
            "input": "i",
            "output": "o",
            "sequence": "s",
            "parameters": "p",
            "frame": "f",
            "graph": "g",
            "workspace": "w",
            "symmetry": "%",
            "crystal": "#",
            "disable_subrotamers": "v",
            "omit_hydrogens": "h",
            "disable_terminal_capping": "t",
        })

    To build a command line, call ``SCWRL4`` with the arguments to actually call
    the executable. The command line can the be run using
    :py:func:`subprocess.run` as usual:

    .. code-block:: python

        scwrl_cmd_line = SCWRL4(options={
            "input": "foo.pdb",
            "output": "foo.scwrl4",
            # ... etc
        })
        subprocess.run(scwrl_cmd_line)

    .. seealso::

        :py:class:`phyre_engine.tools.external.ExternalTool`
            Wrapper for third-party tools.

        :py:mod:`phyre_engine.tools.hhsuite.tool`
            Pre-defined wrappers for hhsuite tools.


.. _example-hhsuite-addpsipred:

Component :py:class:`.hhsuite.AddPsipred`
-----------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.AddPsipred

.. code-block:: none
    :caption: Change in pipeline state

    None

.. code-block:: none
    :caption: Console output (trimmed horizontally)

    INFO     : phyre_engine.pipeline.Pipeline : Running component 4: phyre_engine.component.hhsuite.AddPsipred
    DEBUG    : phyre_engine.component.hhsuite.AddPsipred : Running ['…/share/hhsuite/scripts/addss.pl', '-i', 'query.a3m']
    $ cp query.a3m /tmp/nmL6wlQYG9/FfN67_4d_6.in.a3m
    Filtering alignment to diversity 7 ...
    $ hhfilter -v 1 -neff 7 -i /tmp/nmL6wlQYG9/FfN67_4d_6.in.a3m -o /tmp/nmL6wlQYG9/FfN67_4d_6.in.a3m
    $ …/share/hhsuite//scripts/reformat.pl -v 1 -r -noss a3m psi /tmp/nmL6wlQYG9/FfN67_4d_6.in.a3m /tmp/nmL6wlQYG9/FfN67_4d_6.in.psi
    Predicting secondary structure with PSIPRED ... $ …/bin/blastpgp -b 1 -j 1 -h 0.001 …
    $ echo FfN67_4d_6.chk > /tmp/nmL6wlQYG9/FfN67_4d_6.pn

    $ echo FfN67_4d_6.sq  > /tmp/nmL6wlQYG9/FfN67_4d_6.sn

    $ /data/phyreenginedev/conda/envs/phyreengine/bin/makemat -P /tmp/nmL6wlQYG9/FfN67_4d_6
    $ /data/phyreenginedev/conda/envs/phyreengine/bin/psipred /tmp/nmL6wlQYG9/FfN67_4d_6.mtx …
    $ /data/phyreenginedev/conda/envs/phyreengine/bin/psipass2 …
    done

HHBlits is able to take advantage of predicted secondary structure information
to slightly improves its alignment accuracy, so this component uses `PSIPRED
<http://bioinf.cs.ucl.ac.uk/psipred/>`_ to add secondary structure information
to the top of the MSA generated in the previous step. This component is a
simple wrapper around the :command:`addss.pl` script distributed with hh-suite,
and relies on legacy NCBI blast (*groan*: it's 2018, people). No changes are
made to the pipeline state, because the existing A3M file is modified.

.. _example-hhsuite-hhsearch:

Component :py:class:`.hhsuite.HHSearch`
---------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.HHSearch:
        database: !template '{FOLDLIB_PREFIX}/foldlib'
        input_type: a3m
        options:
          output: hhsearch_search.hhr
          atab: hhsearch_search.atab
          alt: 1
          z: 100
          b: 100
          Z: 50000
          B: 50000
          E: 10
          Ofas: alignments.fasta

.. literalinclude:: /homology_modelling_session/dump-06.pp
    :diff: /homology_modelling_session/dump-05.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output (snipped and pretty-printed)

    INFO     : phyre_engine.pipeline.Pipeline : Running component 5: phyre_engine.component.hhsuite.HHSearch
    DEBUG    : phyre_engine.component.hhsuite.HHSearch : Running [
        '/data/phyreenginedev/conda/envs/phyreengine/bin/hhsearch',
        '-v', '2',
        '-cpu', '16',
        '-o', 'hhsearch_search.hhr',
        '-atab', 'hhsearch_search.atab',
        '-alt', '1',
        '-z', '100',
        '-b', '100',
        '-Z', '50000',
        '-B', '50000',
        '-E', '10',
        '-Ofas', 'alignments.fasta',
        '-d', '/bmm/phyreengine/foldlib_mk2/foldlib',
        '-i', 'query.a3m']
    - 19:21:51.979 INFO: query.a3m is in A2M, A3M or FASTA format

    - 19:21:52.073 WARNING: MSA Q96NU7 looks too diverse (Neff=11.23>11). Better check it with an alignment viewer for non-homologous segments. Also consider building the MSA with hhblits using the - option to limit MSA diversity.

    - 19:21:52.081 INFO: Searching 179354 database HHMs without prefiltering

    - 19:21:52.165 INFO: Iteration 1

    - 19:21:52.250 WARNING: database contains sequences that exceeds maximum allowed size (maxres = 20001). Maxres can be increased with parameter -maxres.

    - 19:21:52.388 INFO: HMMs passed 2nd prefilter (gapped profile-profile alignment)   : 179354

    - 19:21:52.388 INFO: HMMs passed 2nd prefilter and not found in previous iterations : 179354

    - 19:21:52.388 INFO: Scoring 179354 HMMs using HMM-HMM Viterbi alignment

    - 19:21:53.278 INFO: Alternative alignment: 0

    - 19:22:22.807 INFO: 179354 alignments done

    - 19:22:24.847 INFO: Realigning 579 HMM-HMM alignments using Maximum Accuracy algorithm

    - 19:22:25.495 INFO: 482 sequences belonging to 482 database HMMs found with an E-value < 0.001

    Query         Q96NU7
    Match_columns 426
    No_of_seqs    141 out of 1690
    Neff          11.23
    Searched_HMMs 179354
    Date          Wed May 16 19:22:25 2018
    Command       /data/phyreenginedev/conda/envs/phyreengine/bin/hhsearch -v 2 -cpu 16 -o hhsearch_search.hhr -atab hhsearch_search.atab -alt 1 -z 100 -b 100 -Z 50000 -B 50000 -E 10 -Ofas alignments.fasta -d /bmm/phyreengine/foldlib_mk2/foldlib -i query.a3m

     No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
      1 2bb0_A                         100.0 4.8E-42 2.7E-47  324.2  40.6  410    3-426     2-413 (413)
      2 2G3F_A                         100.0 5.8E-42 3.2E-47  323.7  41.0  410    3-426     3-414 (414)
      3 2PUZ_A                         100.0   7E-39 3.9E-44  302.0  36.8  399    6-424     2-404 (404)
      4 2OOF_A                         100.0 1.3E-37 7.4E-43  292.2  40.7  387    3-424     1-392 (393)
      5 2Q09_A                         100.0 1.9E-37 1.1E-42  291.1  41.6  385    5-425     2-392 (392)

    <✂✂✂>

After building an MSA for our query sequence, we align it to a library of HMMs
built from the sequences of "template" proteins with known structure. Later, we
will use the alignment between the query sequence and the template sequence to
map the query positions onto the known structure.

This component takes the same parameters as :ref:`example-hhsuite-hhblits`.
Internally, the two components components actually inherit from the same
:py:class:`base class <.hhsuite.HHSuiteTool>`. As before, we can see the
options that are passed to the component are reflected in command line
parameters passed to the :command:`hhsearch` tool, which are logged at the
``DEBUG`` level. This illustrates the way options are passed through to
external tools: the ``z``, ``b``, etc options are passed verbatim to hhsuite.


.. _example-hhsuite-a3mssparser:

Component :py:class:`.hhsuite.A3MSSParser`
------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.A3MSSParser

.. literalinclude:: /homology_modelling_session/dump-07.pp
    :diff: /homology_modelling_session/dump-06.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 6: phyre_engine.component.hhsuite.A3MSSParser

The :command:`addss.pl` script run by the :ref:`.hhsuite.AddPsipred component
<example-hhsuite-addpsipred>` writes secondary structure information directly
to the MSA, but that information could be of interest to a user. This component
reads secondary structure information from the A3M file and stores it in the
``secondary_structure_sequence`` field. Both the prediction confidence
(``ss_conf``) and predicted state (``ss_pred``) are stored.


.. _example-hhsuite-reportparser:

Component :py:class:`.hhsuite.ReportParser`
-------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.ReportParser

.. literalinclude:: /homology_modelling_session/dump-08.pp
    :diff: /homology_modelling_session/dump-07.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 7: phyre_engine.component.hhsuite.ReportParser

After running hhsearch, the report file :file:`hhsearch_search.hhr` contains the
results of the search. The :py:class:`.hhsuite.ReportParser` component will
examine that report file (pointed to by the ``report`` key in the pipeline
state) and create a ``templates`` list.

The names of the elements in the ``templates`` list are similar to the names of
the fields in the report file. :py:class:`~.hhsuite.ReportParser` will read data
from the summary section and alignment section of the report file, but it does
not read the query-template alignments from the report, because the alignments
are stored in a format that is non-trivial to parse. Thankfully, modern versions
of hhsearch include the ``-atab`` option to generate machine-readable alignments.

Component :py:class:`.hhsuite.TabularParser`
--------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.TabularParser

.. literalinclude:: /homology_modelling_session/dump-09.pp
    :diff: /homology_modelling_session/dump-08.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 8: phyre_engine.component.hhsuite.TabularParser

This component acts similarly to :ref:`the previous component
<example-hhsuite-reportparser>`, but instead parses machine-readable alignments
generated with the ``-atab`` option to hhsearch. The name of the file to be read
is taken from the ``atab`` field in the pipeline state.

The alignment between the query and the template is represented as a list of
:math:`ij` pairs, stored as (named) tuples along with any other information
parsed from the alignment file. In this case, we have a residue-residue
alignment score, a secondary structure similarity score (always reported as 0
because of a bug in hhsearch), a posterior alignment probability, and the
DSSP-assigned secondary structure state of the template.

For each aligned residue pair, the :math:`i` and :math:`j` indices can be
accessed either by name or by index:

.. code-block:: python

    i = residue_pair.i
    j = residue_pair.j
    i, j = residue_pair[0:2]

Component :py:class:`.hhsuite.FastaParser`
------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .hhsuite.FastaParser:
        ignore: ["Consensus"]


.. literalinclude:: /homology_modelling_session/dump-10.pp
    :diff: /homology_modelling_session/dump-09.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 9: phyre_engine.component.hhsuite.FastaParser

It is possible to generate textual representations of the alignments parsed by
the previous component, but it is easier to let hhsearch take the burden. When
we called :ref:`hhsearch earlier <example-hhsuite-hhsearch>`, we passed
the ``-Ofas`` option to tell hhsearch to write alignments in FASTA format. The
sequences were written to :file:`alignments.fasta`, which is pointed to by the
``pairwise_fasta`` field in the pipeline state. This component parses the FASTA
file and stores the textual sequence alginments in the ``sequence_alignments``
field.

The "consensus" sequence written by hhsearch is of no interest to us, so we
include it in the ``ignore`` list.

Component :py:class:`.hhsuite.PSSM`
-----------------------------------

.. code-block:: yaml
    :caption: Excerpt

    # Generate position-specific scoring matrix from query MSA.
    - .hhsuite.PSSM:
        hhsuite_dir: !template '{BIN_PREFIX}/bin'
        HHLIB: !template '{BIN_PREFIX}/share/hhsuite/'
        blast_dir: !template '{BIN_PREFIX}/bin'

.. literalinclude:: /homology_modelling_session/dump-11.pp
    :diff: /homology_modelling_session/dump-10.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output


    INFO     : phyre_engine.pipeline.Pipeline : Running component 10: phyre_engine.component.hhsuite.PSSM
    - 19:22:35.895 INFO: Input file = /tmp/phyreengine-dhmdzc4a-pssm/msa.a3m

    - 19:22:35.895 INFO: Output file = /tmp/phyreengine-dhmdzc4a-pssm/msa.a3m

    Done
    Reformatted /tmp/phyreengine-dhmdzc4a-pssm/msa.a3m with 111 sequences from a3m to psi and written to file /tmp/phyreengine-dhmdzc4a-pssm/msa.psi
    [blastpgp] WARNING: -t larger than 1 not supported when restarting from a checkpoint; setting -t to 1

    BLASTP 2.2.22 [Sep-27-2009]
    ...

Some components—notably :py:class:`.disorder.Disopred` and :py:class:`.modelling.LoopModel`
—require a PSSM as input, either as an :file:`.mtx` or in text format. We don't
currently have an elegant way of generating PSSMs, so this component resorts to
the old school and uses legacy :command:`blastpgp` and :command:`makemat` to
build PSSMs in multiple different formats. Yes, yes, I know. It's
``$CURRENT_YEAR``. Sorry.

.. note::

    This is the only tool in PhyreEngine that uses legacy BLAST (other than
    indirectly via :py:class:`~.AddPsipred`, which we can't be blamed for), and
    we should really look for a replacement. Because of this, there is no
    configuration section for legacy blast and this tool must be configured in
    the pipeline.

Component :py:class:`.disorder.Disopred`
----------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .disorder.Disopred

.. literalinclude:: /homology_modelling_session/dump-12.pp
    :diff: /homology_modelling_session/dump-11.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 11: phyre_engine.component.disorder.Disopred
    INFO     : phyre_engine.component.disorder.Disopred : Predicting disorder with DISOPRED2.
    DEBUG    : phyre_engine.component.disorder.Disopred : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/disopred2', '/tmp/tmp2ffup013/disopred2', 'profile.mtx', '/data/phyreenginedev/conda/envs/phyreengine/share/disopred/data/', '5']
    INFO     : phyre_engine.component.disorder.Disopred : Running neural network classifier.
    DEBUG    : phyre_engine.component.disorder.Disopred : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/diso_neu_net', '/data/phyreenginedev/conda/envs/phyreengine/share/disopred/data/weights.dat.nmr_nonpdb', 'profile.mtx'] > /tmp/tmp2ffup013/diso_neu_net
    INFO     : phyre_engine.component.disorder.Disopred : Running nearest neighbour classifier.
    DEBUG    : phyre_engine.component.disorder.Disopred : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/diso_neighb', 'profile.mtx', '/data/phyreenginedev/conda/envs/phyreengine/share/disopred/data/dso.lst'] > /tmp/tmp2ffup013/diso_neighb
    INFO     : phyre_engine.component.disorder.Disopred : Combining disordered residue predictions.
    DEBUG    : phyre_engine.component.disorder.Disopred : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/combine', '/data/phyreenginedev/conda/envs/phyreengine/share/disopred/data/weights_comb.dat', '/tmp/tmp2ffup013/disopred2.diso', '/tmp/tmp2ffup013/diso_neu_net', '/tmp/tmp2ffup013/diso_neighb'] > disorder.diso

Generates a disorder prediction using `DISOPRED
<http://bioinf.cs.ucl.ac.uk/web_servers/disopred/disopred_overview/>`_. Disorder
predictions are stored in the ``disorder`` field of the pipeline state. To allow
for the use of multiple predictors, the ``disorder`` field is further indexed by
the name of the tool that produced the prediction. Residue-level predictions are
stored as dictionaries with an ``assigned`` and ``confidence`` key, indicating
which state was assigned and the confidence of either state.

.. _example-sort-sort:

Component :py:class:`.sort.Sort`
--------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .sort.Sort:
        field: templates
        keys: [{key: 'sum_probs', reverse: True}]

.. literalinclude:: /homology_modelling_session/dump-13.pp
    :diff: /homology_modelling_session/dump-12.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 12: phyre_engine.component.sort.Sort

This component will sort a list in the pipeline state according to some key.
This component was configured in the same way as the example given above, which
means "sort the ``templates`` field in descending order, using the
``sum_probs`` key as the sort key." The component first retrieves the list to
be sorted (``templates``), then looks up the ``sum_probs`` key for each element
in the list. List elements are sorted according to the value of the
``sum_probs`` field, in descending order because ``reverse`` is ``True``.

If the sort key contains duplicates, ties can be broken by supplying more keys
in the ``keys`` list. For example, we could use the following to first sort by
PDB ID and then break ties on chain ID:

.. code-block:: yaml

    - .sort.Sort:
        field: templates
        keys:
        - {key: 'PDB'}
        - {key: 'chain'}

Component :py:class:`foldlib.Open`
----------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .foldlib.Open

.. literalinclude:: /homology_modelling_session/dump-14.pp
    :diff: /homology_modelling_session/dump-13.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 13: phyre_engine.component.foldlib.Open

At this point, we have a list of hits, but we lack any metadata, and the hits
are taken from our HMM library, which is grouped at 100% sequence identity. We
can extract the information we want from the :ref:`fold library <foldlib>`.
Opening and closing database connections can be an expensive task, so this
component is used to open a persistent database connection that is stored in the
``template_db`` field of the pipeline state.

.. warning::

    The ``template_db`` is a :py:class:`~.TemplateDatabase` object. Note that
    this object cannot be pickled, which means it cannot be transferred to a
    remote machine using the components in :py:mod:`phyre_engine.component.pbs`.
    This makes sense: how would we sensibly serialise an open database
    connection?

Component :py:class:`.foldlib.Map` (PDB and chain ID parsing)
-------------------------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .foldlib.Map:
        components:
        - .db.metadata.ParseField:
            field: name
            regex: '^(?P<PDB>\w{4})_(?P<chain>\S+)$'

.. literalinclude:: /homology_modelling_session/dump-15.pp
    :diff: /homology_modelling_session/dump-14.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 14: phyre_engine.component.foldlib.Map
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 0 / 579
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.db.metadata.ParseField
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 1 / 579
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.db.metadata.ParseField
    ...

Any component with the class name ``Map`` is a subclass of
:py:class:`~.component.Map`. A ``Map`` component loops over a list of elements
in the pipeline state, and applies a pipeline of components to each element in
the list. The :py:class:`~.foldlib.Map` component is a simple subclass that
always loops over the ``templates`` list and temporarily copies the
``template_db`` field into each child element:

.. code-block:: yaml

    # The follwowing two components do the same job:
    - .foldlib.Map:
        components: [...]

    - .component.Map:
        # Loop over the "templates" list
        field: templates
        # Temporarily copy the "template_db" field to each child
        copy: [template_db]
        components: [...]

Here, we are looping over each hit in the ``templates`` list and running the
:py:class:`.db.metadata.ParseField` component. The :py:class:`~.ParseField`
component is very useful: it takes a ``field`` to examine, and runs the regex
given in the ``regex`` parameter against that field. Any matching named captures
in the regular expression are copied into the pipeline state.

First, the :py:class:`.foldlib.Map` copies the ``template_db`` field into each
element of the ``templates`` list. Because :py:class:`.db.metadata.ParseField`
was run inside the :py:class:`.foldlib.Map`, it only sees a single element of
the ``templates`` list and treats it like the whole pipeline state. When the
:py:class:`.db.metadata.ParseField` component is run, it reads the ``name``
field of the pipeline state, and extracts the ``PDB`` and ``chain`` fields by
regex matching. Finally, when the loop is complete the :py:class:`.foldlib.Map`
components delete the temporary ``template_db`` field and joins the
``templates`` list back together.


Component :py:class:`foldlib.Map` (representative expansion)
------------------------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .foldlib.Map:
        components:
        - .foldlib.ExpandSequenceRepresentatives

.. literalinclude:: /homology_modelling_session/dump-16.pp
    :diff: /homology_modelling_session/dump-15.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 15: phyre_engine.component.foldlib.Map
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.foldlib.ExpandSequenceRepresentatives
    INFO     : phyre_engine.component.foldlib.ExpandSequenceRepresentatives : 2G3F_A expanded to 2 members
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 1 / 579
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.foldlib.ExpandSequenceRepresentatives
    INFO     : phyre_engine.component.foldlib.ExpandSequenceRepresentatives : 2bb0_A expanded to 2 members
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 2 / 579
    ...

Now that we have the PDB and chain ID for each template, we can expand the list
of sequence representatives. This is done by the
:py:class:`.foldlib.ExpandSequenceRepresentatives` component, which looks up
each template from the fold library according to its PDB and chain ID. The
:py:class:`.foldlib.ExpandSequenceRepresentatives` assumes that since the
sequences of the representative templates and the expanded members are
identical, the scores of the representative can be copied verbatim into each
member. In the console output above, we can see that both templates were
expanded into members.

Since the ``diff`` output given above is a bit mangled, we could take a look at
the expansions of these templates using the :py:class:`~.TemplateDatabase` class
directly:

.. code-block:: python

    >>> from phyre_engine.tools.template import TemplateDatabase
    >>> tdb = TemplateDatabase("/path/to/foldlib/db", "/path/to/foldlib/chains")
    >>> tdb.expand_seq_reps("2G3F", "A")
    [('2g3f', 'B')]
    >>> tdb.expand_seq_reps("2bb0", "A")
    [('2bb0', 'B')]

.. warning::

    The :py:class:`.foldlib.ExpandSequenceRepresentatives` returns a *list* of
    items, because a representative may expand into arbitrarily many chains
    (`3J3Q <https://www.rcsb.org/structure/3J3Q>`_ is a notorious example).
    **This means that it must be the *last* component in a Map**, so that the
    map can take care of merging the resulting lists.

.. note::

    HHSearch includes a secondary-structure score in its matching algorithm, so
    the assumption that two templates with identical sequence will have
    identical scores is not strictly true. A more rigorous approach would be to
    group chains by 100% identical sequence *and* secondary structure. In
    practice, the effect is likely to be minimal.

Component :py:class:`foldlib.Map` (metadata retrieval)
------------------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .foldlib.Map:
        components:
        - .foldlib.Metadata

.. literalinclude:: /homology_modelling_session/dump-17.pp
    :diff: /homology_modelling_session/dump-16.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 16: phyre_engine.component.foldlib.Map
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 0 / 1170
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.foldlib.Metadata
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 1 / 1170
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.foldlib.Metadata
    DEBUG    : phyre_engine.component.foldlib.Map : Runing pipeline 2 / 1170

This loop again uses the database connection stored in the ``template_db`` field
to look up metadata for each template. This metadata is not used in the
pipeline, but is valuable information for the user.

Component :py:class:`foldlib.Rollback`
--------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .foldlib.Rollback

.. literalinclude:: /homology_modelling_session/dump-18.pp
    :diff: /homology_modelling_session/dump-17.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 17: phyre_engine.component.foldlib.Rollback

At this point, we are finished with our connection to the fold library and
should get rid of it. There are two components that can do this:
:py:class:`~.foldlib.Commit` and :py:class:`~.foldlib.Rollback`. These are named
according to the SQL command they run: :py:class:`~.foldlib.Commit` will commit
the current transaction to the database, and :py:class:`~.foldlib.Rollback` will
roll back any changes we have made.

In our case, we haven't made any changes, but we still specify the
:py:class:`~.foldlib.Rollback` component to make it clear that this pipeline is
not intended to affect the fold library in any way. This is why there is no
``.foldlib.Close`` component: it is better to be explicit about the way in which
the connection is closed.

Component :py:class:`.cluster.cluster.EM4GMM`
---------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .cluster.EM4GMM:
        select_expr: templates
        dimensions_expr: '[start(query_range), stop(query_range)]'
        num_components: 5
        merge: 0.01
        bin_dir: !template '{BIN_PREFIX}/bin'

.. literalinclude:: /homology_modelling_session/dump-19.pp
    :diff: /homology_modelling_session/dump-18.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 18: phyre_engine.component.cluster.cluster.EM4GMM
    DEBUG    : phyre_engine.component.cluster.cluster.EM4GMM : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/gmmtrain', '-d', '/tmp/tmplem5bvci', '-m', 'model.gmm', '-r', 'train.log.json', '-n', '5', '-u', '0.01']
    Number of Components: 000005
    Iteration: 00001    Improvement: 100%    LogLikelihood: -13.410
    Iteration: 00002    Improvement:  19%    LogLikelihood: -10.832
    Iteration: 00003    Improvement:   6%    LogLikelihood: -10.159
    Iteration: 00004    Improvement:   9%    LogLikelihood: -9.276
    Iteration: 00005    Improvement:   8%    LogLikelihood: -8.535
    Iteration: 00006    Improvement:  10%    LogLikelihood: -7.705
    Iteration: 00007    Improvement:   4%    LogLikelihood: -7.421
    Iteration: 00008    Improvement:   1%    LogLikelihood: -7.356
    Iteration: 00009    Improvement:   0%    LogLikelihood: -7.345
    Iteration: 00010    Improvement:   1%    LogLikelihood: -7.285
    Iteration: 00011    Improvement:   7%    LogLikelihood: -6.754
    Iteration: 00012    Improvement:   5%    LogLikelihood: -6.440
    Iteration: 00013    Improvement:   2%    LogLikelihood: -6.329
    Iteration: 00014    Improvement:   1%    LogLikelihood: -6.272
    Iteration: 00015    Improvement:   0%    LogLikelihood: -6.252
    Iteration: 00016    Improvement:   0%    LogLikelihood: -6.225
    Iteration: 00017    Improvement:   0%    LogLikelihood: -6.225
    Number of Components: 000003
    Iteration: 00001    Improvement: 100%    LogLikelihood: -9.799
    Iteration: 00002    Improvement:  22%    LogLikelihood: -7.666
    Iteration: 00003    Improvement:   0%    LogLikelihood: -7.663
    Number of Components: 000002
    Iteration: 00001    Improvement: 100%    LogLikelihood: -10.558
    Iteration: 00002    Improvement:  11%    LogLikelihood: -9.400
    Iteration: 00003    Improvement:   0%    LogLikelihood: -9.400
    Number of Components: 000002
    DEBUG    : phyre_engine.component.cluster.cluster.EM4GMM : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/gmmclass', '-d', '/tmp/tmplem5bvci', '-m', 'model.gmm', '-r', 'classify.log.json']
    Score: -9.4007792108

Runs `em4gmm <https://github.com/juandavm/em4gmm>`_ via the
:py:class:`.cluster.EM4GMM` component. This component is configured to cluster
the elements in the ``templates`` list using the dimensions defined by the start
and end of the template coverage.

The ``select_expr`` field is a `JMESPath <http://jmespath.org/>`_ expression
evaluated against the pipeline state that must gives the list of items to
cluster. In this case, it is just saying "take the ``templates`` list". The
``select_expr`` parameter is another JMESPath expression that is evaluated on
each element of ``templates``. Here, it says to cluster based on the first and
last query residue covered by the hit.

.. note::

    PhyreEngine implements several extensions to JMESPath in the
    :py:class:`phyre_engine.tools.jmespath.JMESExtensions` class. In this
    example, the ``start`` and ``stop`` functions are extensions that operate on
    a python ``range``.

Component :py:class:`.set.Union`
--------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .set.Union:
        sets: [templates]
        key: cluster
        destination: cluster_reps

.. literalinclude:: /homology_modelling_session/dump-20.pp
    :diff: /homology_modelling_session/dump-19.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 19: phyre_engine.component.set.Union

At this point in the pipeline, we have a big list of ``templates``, each of
which has a ``cluster`` field assigned by em4gmm. We now want to pick the best
model assigned to each cluster so that we can operate only on those models. We
don't want to discard the other hits, as they will still be useful to the user.

This could probably benefit from a separate component (``TopNPerGroup``, maybe),
but in the simple case here where we only want to pick one model from each
cluster, we can use the :py:class:`.set.Union`. component.

Components in the :py:mod:`.set` module operate on "sets" in the pipeline state.
Here, a "set" is a *unique* subsection of a list indexed by some selection of
fields. For example, em4gmm defined two clusters for these models: 0 and 1. If
we define sets based on the ``cluster`` field of elements in the ``templates``
list, we will end up with only two items, one with a ``cluster`` of 0 and one
with a ``cluster`` of 1. Since set operations preferentially keep the first
unique item seen and :py:class:`.set.Union`'s job is to combine sets, we end up
with a single list containing the top-ranked models in each cluster!

To illusrate this, consider the following pipeline state:

.. code-block:: python

    {"templates": [{"cluster": 0, "name": "0A"},
                   {"cluster": 1, "name": "1A"},
                   {"cluster": 0, "name": "0B"},
                   {"cluster": 1, "name": "1B"}]}

In our usage (shown above), :py:class:`.set.Union` first selects the
``templates`` list by applying the JMESPath expresions given in the ``sets``
parameter. We could have specified multiple lists here, but we are only using
the set operation to pick an item.

Next, sets are generated based on the ``key`` parameter. This is again a
JMESPath expression, this time evaluated in the context of each item in the
``templates`` list. Here, we are just selecting the ``cluster`` field. At this
point, the ``templates`` list looks like this:

.. code-block:: python

    [{"cluster": 0, "name": "0A"},
     {"cluster": 1, "name": "1A"}]

Note that this is calculated internally and does not overwrite the ``templates``
list.

Finally, the resulting list is saved in the ``cluster_reps`` field, which was
chosen with the ``destination`` parameter. The pipeline state will now look like
this:

.. code-block:: python

    {"templates": [{"cluster": 0, "name": "0A"},
                   {"cluster": 1, "name": "1A"},
                   {"cluster": 0, "name": "0B"},
                   {"cluster": 1, "name": "1B"}],
     "cluster_reps": [{"cluster": 0, "name": "0A"},
                      {"cluster": 1, "name": "1A"}]}


Component :py:class:`.component.Map` (modelling)
------------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .component.Map:
        field: cluster_reps
        copy:
          pssm: pssm
          sequence: query_sequence
        components:
        - .modelling.HomologyModeller
        - .modelling.LoopModel
        - .sidechain.Scwrl4

.. literalinclude:: /homology_modelling_session/dump-21.pp
    :diff: /homology_modelling_session/dump-20.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 20: phyre_engine.component.component.Map
    DEBUG    : phyre_engine.component.component.Map : Runing pipeline 0 / 2
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.modelling.HomologyModeller
    DEBUG    : phyre_engine.component.modelling.HomologyModeller : Creating model file 93-2IMR_A.pdb
    INFO     : phyre_engine.pipeline.Pipeline : Running component 1: phyre_engine.component.modelling.LoopModel
    DEBUG    : phyre_engine.component.modelling.LoopModel : Loop modelling using tmpdir: /tmp/phyreengine-ln1i01ud-loop
    DEBUG    : phyre_engine.component.modelling.LoopModel : Running ['/data/phyreenginedev/conda/envs/phyreengine/bin/nova', '-c', '/bmm/phyreengine/share/loop/loop.config', '-pssm', '/tmp/phyreengine-ln1i01ud-loop/loop.pssm', '-f', '/tmp/phyreengine-ln1i01ud-loop/query.fasta', '-l', '/tmp/phyreengine-ln1i01ud-loop/model.list', '-d', '93-2IMR_A.loop']
    ________________________________________________________________________________

                                           Nova

                     Protein structure prediction by fragment assembly
                                       (Ver. 0.6.6)
    ________________________________________________________________________________

    INFO     : phyre_engine.pipeline.Pipeline : Running component 2: phyre_engine.component.sidechain.Scwrl4

Finally, we can actually do some modelling! It's amazing how much of a "homology
modelling" pipeline is bookkeeping. Here, we loop over each of the hits in the
``cluster_reps`` list. We first generate a "crude" model using
:py:class:`.modelling.HomologyModeller`, fill in as many loops as possible using
:py:class:`.modelling.LoopModel`, and then repack side-chains using
:py:class:`.sidechain.Scwrl4`.

.. note::

    These tools make no guarantees about *where* the modelled files will end up.
    Use the ``model`` field of the pipeline state to determine the file name. In
    the future, these components will likely be changed to produce more friendly
    file names.

Component :py:class:`.component.Conditional` (TMScore)
------------------------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .component.Conditional:
        field: native
        components:
        - .component.Map:
            field: cluster_reps
            copy: [native]
            components:
            - .strucaln.TMScore

.. code-block:: yaml
    :caption: Change in pipeline state

    No change

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 21: phyre_engine.component.component.Conditional

The :py:class:`~.component.Conditional` is similar to
:py:class:`~.component.Map` in the sense that they both run a child pipeline.
The :py:class:`~.component.Conditional` component takes a ``field`` parameter
that it examines for "`truthiness
<https://docs.python.org/3/library/stdtypes.html#truth-value-testing>`_": if the
field does not exist in the pipeline state or is a value that evaluates to
`False` in a boolean context, then the child pipeline is not executed.

In this case, if the ``native`` field is present in the pipeline state, the
embedded :py:class:`.component.Map` would loop over the ``templates`` list run
:py:class:`~.component.strucaln.TMScore`.  The
:py:class:`~.component.strucaln.TMScore` component would run the
:command:`TMScore` executable to find the structural similarity between each
model and the native structure. To supply a native structure, add an extra
:option:`--start` option.

Component :py:class:`set.Union`
-------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .set.Union:
        sets: [cluster_reps, templates]
        key: '[PDB, chain, rank]'
        destination: templates

.. literalinclude:: /homology_modelling_session/dump-23.pp
    :diff: /homology_modelling_session/dump-22.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 22: phyre_engine.component.set.Union

At this point in the pipeline, we have two lists containing models: the
``cluster_reps`` list contains the cluster representatives that have been
modelled, and the ``templates`` list contains the remainder of the hits. To
combine the two, we treat the two lists as sets and merge them with
:py:class:`.set.Union`.

This time, we use the combination of the fields ``PDB``, ``chain`` and ``rank``
to uniquely identify each hit. Then, we merge the ``cluster_reps`` and
``templates`` list together based on that unique key. Lists on the left take
precedence, so the hits in the ``cluster_reps`` list overwrite the correponding
hits in the ``templates`` list.

.. note::

    The ``key`` parameter is a single JMESPath expression, so it must be
    specified as a list. YAML treats square brackets as denoting a list, so be
    sure to enclose the ``key`` field in single quotes to force the YAML parser
    to treat it as a string. The ``sets`` parameter is a list of JMESPath
    expressions, so it must be specified as a list of strings.

Component :py:class:`sort.Sort`
-------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .sort.Sort:
        field: templates
        keys: [{key: 'sum_probs', reverse: True}]

.. literalinclude:: /homology_modelling_session/dump-24.pp
    :diff: /homology_modelling_session/dump-23.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 23: phyre_engine.component.sort.Sort

Set operations leave the resulting list unordered, so we again sort the
templates by the ``sum_probs`` field. This component is identical to the
:ref:`previous Sort component <example-sort-sort>`.

Component :py:class:`pymol.Init`
--------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .pymol.Init:
        # Load python functions we want to use
        pymol: !template '{BIN_PREFIX}/bin/pymol'
        command: !template 'run {CONFIG_PREFIX}/pymol/phyre_engine.py'
        quiet: False

.. literalinclude:: /homology_modelling_session/dump-25.pp
    :diff: /homology_modelling_session/dump-24.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 24: phyre_engine.component.pymol.Init
    INFO     : phyre_engine.component.pymol.Init : Starting pymol server on port 45800: (Command: ['/data/phyreenginedev/conda/envs/phyreengine/bin/pymol', '-c', '-K', '-d', 'import pymol.rpc; pymol.rpc.launch_XMLRPC(port=45800, nToTry=1); if not pymol.rpc.serv: cmd.quit(1);'])
    DEBUG    : phyre_engine.component.pymol.Init : Pymol exit status: None
    DEBUG    : phyre_engine.component.pymol.Init : Pymol connection error on port 45800
     PyMOL(TM) Molecular Graphics System, Version 2.1.0.
     Copyright (c) Schrodinger, LLC.
     All Rights Reserved.

        Created by Warren L. DeLano, Ph.D.

        PyMOL is user-supported open-source software.  Although some versions
        are freely available, PyMOL is not in the public domain.

        If PyMOL is helpful in your work or study, then please volunteer
        support for our ongoing efforts to create open and affordable scientific
        software by purchasing a PyMOL Maintenance and/or Support subscription.

        More information can be found at "http://www.pymol.org".

        Enter "help" for a list of commands.
        Enter "help <command-name>" for information on a specific command.

     Hit ESC anytime to toggle between text and graphics.

    PyMOL>import pymol.rpc; pymol.rpc.launch_XMLRPC(port=45800, nToTry=1); if not pymol.rpc.serv: cmd.quit(1);
    DEBUG    : phyre_engine.component.pymol.Init : Pymol exit status: None
    PyMOL>1 + 1
    DEBUG    : phyre_engine.component.pymol.Init : Connected to port 45800
    PyMOL>run /bmm/phyreengine/share/pymol/phyre_engine.py

Now that we have a list of hits, some of which have been modelled, we wish to
generate some pretty pictures. We can do this using the components in the
:py:mod:`.pymol` module.

This component starts a `PyMOL <https://pymol.org/2/>`_ process that then
listens for commands. We supply PyMol with a ``command`` parameter to tell it to
load a script containing customised rendering commands: you can see the commands
that PyMol runs in its console output above. The long command beginning with
``import`` tells PyMol to listen for commands on a port; the ``1 + 1`` command
is there to test for connectivity, and the ``run /bmm/...`` command is
set with the ``command`` parameter to the component.

Component :py:class:`.pymol.Run`
--------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .pymol.Run:
        # Load colour scheme
        commands:
        - !template '@{CONFIG_PREFIX}/pymol/setup_magma.pml'

.. code-block:: none
    :caption: Change in pipeline state

    No change

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 25: phyre_engine.component.pymol.Run

    PyMOL>load 02-2G3F_A.loop/model.1/model.1.scwrl4.pdb
     CmdLoad: PDB-string loaded into object "model.1.scwrl4", state 1.
    PyMOL>@/bmm/phyreengine/share/pymol/render_magma.pml
     Ray: render time: 1.24 sec. = 2900.4 frames/hour (1.24 sec. accum.).
    PyMOL>png 002-2G3F_A.png
     ScenePNG: wrote 400x400 pixel image to file "002-2G3F_A.png".
    PyMOL>delete all

We now have a persistent PyMol instance that we can use to run commands without
all the overhead of starting PyMol each time. The first command we want to run
just loads some extra colour definitions from a script.

Component :py:class:`.component.Map`
------------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .component.Map:
        field: templates
        copy: [pymol]
        components:
        - .component.Conditional:
            field: model
            components:
            - .pymol.Run:
                commands:
                - 'load {model}'
                - !template '@{CONFIG_PREFIX}/pymol/render_magma.pml'
                - 'png {rank:03d}-{PDB}_{chain}.png'
                - 'delete all'
            - .alter.Set:
                field: image
                value: '{rank:03d}-{PDB}_{chain}.png'
                reformat: True

.. literalinclude:: /homology_modelling_session/dump-27.pp
    :diff: /homology_modelling_session/dump-26.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 26: phyre_engine.component.component.Map
    DEBUG    : phyre_engine.component.component.Map : Runing pipeline 0 / 845
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.component.Conditional
    INFO     : phyre_engine.pipeline.Pipeline : Running component 1: phyre_engine.component.alter.Remove
    DEBUG    : phyre_engine.component.component.Map : Runing pipeline 1 / 845
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.component.Conditional
    INFO     : phyre_engine.pipeline.Pipeline : Running component 0: phyre_engine.component.pymol.Run
    PyMOL>load 02-2G3F_A.loop/model.1/model.1.scwrl4.pdb
     CmdLoad: PDB-string loaded into object "model.1.scwrl4", state 1. 
    PyMOL>@/bmm/phyreengine/share/pymol/render_magma.pml
     Ray: render time: 3.16 sec. = 1137.9 frames/hour (3.16 sec. accum.).
    PyMOL>png 002-2G3F_A.png
     ScenePNG: wrote 400x400 pixel image to file "002-2G3F_A.png".
    PyMOL>delete all
    INFO     : phyre_engine.pipeline.Pipeline : Running component 1: phyre_engine.component.alter.Set
    INFO     : phyre_engine.pipeline.Pipeline : Running component 1: phyre_engine.component.alter.Remove
    DEBUG    : phyre_engine.component.component.Map : Runing pipeline 2 / 845
    ...

Here, we iterate over each element in the ``templates`` list, temporarily
copying the ``pymol`` field into each. For each component, we check if it has
the field ``model`` (meaning that it has been homology-modelled), and generate
an image using :py:class:`.pymol.Run`.

The only new component here is :py:class:`.alter.Set`, which has the job of
adding the ``image`` field to each element. This is necessary because the
:py:class:`.pymol.Run` component is just dumbly running commands, and has no
idea that it has written an image. We use :py:class:`.alter.Set` to add an
``image`` field containing the name of the image written by PyMol.

.. note::

    The syntax for this portion of the pipeline might look a little obtuse at
    first glance. The thing to note is that both :py:class:`.pymol.Run` and
    :py:class:`.alter.Set` apply string formatting to their parameters based on
    the contents of the pipeline state. This should not be confused with the
    templates that are applied when the YAML file is loaded. If a string has a
    ``!teplate`` tag, the template is resolved at load time; otherwise, if you
    see braces in a string then it will be resolved at run time.

    In the example above, the commands passed to :py:class:`.pymol.Run` are
    formatted with Python's :py:meth:`str.format` method at run time. The line
    beginning with ``!teplate`` is resolved at load time. When strings are
    resolved at run time, they have all parameters in the pipeline state
    available to them. The line ``png {rank:03d}-{PDB}_{chain}.png`` is
    processed like this:

    .. code-block:: python

        "png {rank:03d}-{PDB}_{chain}.png".format(**pipeline_state)

    That is, the fields ``rank``, ``PDB`` and ``chain`` are taken from
    ``pipeline_state`` and formatted into the string.

We could alternatively (and perhaps more neatly) have run these components
without the :py:class:`~.component.Conditional` on the ``cluster_reps`` list
before merging it with the ``templates`` list.

Component :py:class:`.pymol.Quit`
---------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .pymol.Quit

.. literalinclude:: /homology_modelling_session/dump-28.pp
    :diff: /homology_modelling_session/dump-27.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 27: phyre_engine.component.pymol.Quit
    PyMOL>quit

We no longer need our PyMol server, so we should be good citizens and get rid of
it. If we forget, it will be closed when the Python process ends, but it is best
to close it more gracefully.

Component :py:class:`.dump.Pickle`
----------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .dump.Pickle:
        output: state.pickle

.. literalinclude:: /homology_modelling_session/dump-29.pp
    :diff: /homology_modelling_session/dump-28.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 28: phyre_engine.component.dump.Pickle

To allow for pipelines to continue on where this one leaves off, we will write a
:py:mod:`pickle` file containing the current pipeline state. We can use this
later with a refinment pipeline that can model unmodelled hits.

Component :py:class:`.dump.Csv`
-------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .dump.Csv:
        file: models.csv
        select_expr: >
          templates[].{
            length: length(root().sequence),
            model: model,
            template: name,
            sum_probs: sum_probs,
            prob: prob,
            evalue: evalue,
            qrange_start: start(query_range),
            qrange_stop: stop(query_range),
            score: score,
            similarity: similarity,
            TM: TM
          }

.. literalinclude:: /homology_modelling_session/dump-30.pp
    :diff: /homology_modelling_session/dump-29.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 29: phyre_engine.component.dump.Csv

We have finally finished our pipeline, so we should write our results. The
:py:class:`.dump.Csv` component is used to write a comma-separate file. We use
the JMESPath expression ``select_expr`` to choose the columns to write. The
expression should return a list of maps, each field of which becomes a column.

Component :py:class:`.dump.Json`
--------------------------------

.. code-block:: yaml
    :caption: Excerpt

    - .dump.Json:
        output: models.json
        select_expr: >
          merge(@, {templates: map(
              &except(@, ['alignment']),
              templates[])})

.. literalinclude:: /homology_modelling_session/dump-31.pp
    :diff: /homology_modelling_session/dump-30.pp
    :caption: Change in pipeline state
    :language: python

.. code-block:: none
    :caption: Console output

    INFO     : phyre_engine.pipeline.Pipeline : Running component 30: phyre_engine.component.dump.Json

The comma-separated file that we wrote is useful for analysis with tools such as
R or pandas, but a JSON file will be of more use for generating webpages of
results. We can use the :py:class:`.dump.Json` to write a JSON representation of
the pipeline state. The parameters are the same as for :py:class:`.dump.Csv`,
but the results of ``select_expr`` are less constrained because they do not have
to map to row-oriented data. If ``select_expr`` is not included, the entire
state is dumped.

Here, we are dumping everything except the ``alignment`` field of each template,
which is not necessary for the webpages and bulks out the files considerably.
The ``select_expr`` given here merges the root of the pipeline state (``@``)
with a mapping of ``templates`` with the ``alignment`` field excluded.
