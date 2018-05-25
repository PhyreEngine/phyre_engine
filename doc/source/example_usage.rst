.. _example_usage:

=============
Example usage
=============

The sampe pipeline files in the :file:`pipelines` directory can be used for
some common tasks. First, remember to copy :file:`pipelines/config.yml.sample`
to :file:`~/.config/phyreengine/config.yml` and edit it to match your system.
In this document I will use the path names given in the sample configuration
file.

Fold library
============

Building a fold library
-----------------------

To build a new fold library, first back up the existing template database
:file:`{FOLDLIB_PREFIX}/foldlib.db` and the profile databases
:file:`{FOLDLIB_PREFIX}/foldlib.ff*`. Then, delete the original files so
PhyreEngine has a clean slate. The :file:`profiles`, :file:`mmcif` and
:file:`chains` directory can stay: the pipeline will overwrite files in these
directories.

To generate the fold libraray run the pipeline
:file:`pipelines/fold_library-parallel.yml` from your PBS server:

.. code-block:: console

    $ python -mphyre_engine.run pipelines/fold_library-parallel.yml \
        --start keep_existing_cif:1

For debugging purposes, it may be more useful to the run the serial version
that does not farm work out to nodes via PBS:

.. code-block:: console

    $ python -mphyre_engine.run pipelines/fold_library-serial.yml \
        --start keep_existing_cif:true

.. note::

    The pipeline will attempt to download all CIF files thqat were updated
    since the last fold library update.  If :file:`{FOLDLIB_PREFIX}/foldlib.db`
    is deleted, PhyreEngine will attempt to download *all* CIF files unless
    ``overwrite`` is `True` for :py:class:`~.StructureRetriever`.

    For convenience, :py:class:`~.StructureRetriever` is surrounded by a
    :py:class:`~.ConfigLoader` mapping the pipeline state field
    ``keep_existing_cif`` to ``foldlib.overwrite``.


Updating a fold library
-----------------------

TO update the fold library, just run the pipeline with the same command as
above, but without erasing the template database. The pipeline will examine the
``meta`` table of the database to see when the last update happened, grab any
new CIF files, convert them to chains and rebuild HMMs for any new structures,
overwriting existing HMMs.

Troubleshooting
---------------

If there is an error while upating the database, the SQL database will not be
updated, but files on disk may have changed. Once you have fixed any errors, it
may be worth setting the various ``overwrite`` parameters in the pipeline to
`False` to prevent the pipeline from performing expensive operations like
regenerating profiles or chains that it has already built.

Homology modelling
==================

Running a homology modelling pipeline should just require the following command:

.. code-block:: console

    $ python -mphyre_engine.run pipelines/homology_model.cluster-serial.yml \
        --start working_dir:$DIR \
        --start input:query.fasta

The pipeline :file:`homology_model.cluster-serial.yml` may be replaced with
:file:`homology_model.cluster-parallel.yml` to detach the process onto a node
via PBS.

Troubleshooting
---------------

Remember that this pipeline is designed to be called from a web server. The
``input`` sequence file must exist, and already be present in the
``working_dir``.

Refining homology models
========================

The homology modelling pipeline will only build models for cluster
representatives. The pipelines :file:`pipeline/refine-serial.yaml` and
:file:`pipeline/refine-parallel.yml` can be used to refine a list of models by
their index in the ``models`` field.

Assuming that the homology modelling pipeline successfully wrote
:file:`state.pickle`, models 1–20 can be modelled by running:

.. code-block:: console

    $ python -mphyre_engine.run pipelines/refine-serial.yml \
        --start working_dir:$DIR \
        --start state:state.pickle \
        --start 'templates:@[1:20]'

Troubleshooting
---------------

The ``templates`` parameter supplied as a ``start`` value is a `JMESPath
<http://jmespath.org/>`_ expression evaluated relative to the ``templates``
list in the pipeline state.

Bulk modelling
==============

The bulk modelling pipelines, :file:`pipelines/pipelines/bulk_model-serial.yml`
and :file:`pipelines/bulk_model-parallel.yml`, can be used to model all protein
sequences in a file. The pipelines need to be supplied an ``input`` file and a
regex describing how to parse the sequence name (which is used to name
subdirectories and must be unique for each sequence).

.. code-block:: console

    $ python -mphyre_engine.run ~/phyre_engine/pipelines/bulk_model-parallel.yml \
        --start input:IWGSC_v1.1_HC_20170706_pep.fasta \
        --config 'metadata.regex:^(?P<name>.*)$'

Troubleshooting
---------------

.. warning::

    If you supply a regular expression that produces non-unique sequence names,
    the pipelines will become confused when the second sequence doesn't match
    cached MSAs and alignments generated from the first sequence. This can
    manifest itself in in very non-obvious ways when running the paralle
    pipeline, as processes fight over files with different sequences.

    **This should be the first thing to check when you see obscure errors in
    the logs!** In particular, :py:class:`~.modelling.SoedingSelect` does some
    alignment manipulation that will produce errors about index ranges if
    the size the query-template alignments do not match the size of the query.

Backphyre
=========

The sequence of a protein structure can be scanned against an hhsearch database
using :file:`pipelines/backphyre-serial.yml` and
:file:`pipelines/backphyre-parallel.yml`. 

.. code-block:: console

    $ python -mphyre_engine.run \
      ~/code/phyre_engine/configs/backphyre.yml \
      --start working_dir:$PWD  \
      --start structure:12as_A.pdb \
      --config hhsuite.database:/bmm/phyreengine/data/scop40_01Mar17/scop40 \
      --config hhsuite.database:/bmm/phyreengine/data/Pfam/hhsuite/pfam

Genome daabases
---------------

The bulk modelling pipelines do not currently generate hhseach databases when
the genomes are built. To do so, you would need to add
:py:class:`~.component.hhsuite.HHMake` and
:py:class:`~.component.hhsuite.CSTranslate` after
:py:class:`~.component.hhsuite.HHBlits`, and
:py:class:`~.component.hhsuite.BuildDatabase` later in the pipeline.
