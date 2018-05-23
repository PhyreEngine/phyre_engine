.. _foldlib:

================
The fold library
================

PhyreEngine's fold library contains every structure in the PDB, though we cheat
a little and only generate profiles for every unique sequence. Sequences are
parsed from the ATOM records of the mmCIF files that make up the PDB---*not*
the SEQRES records.

.. note::

    It's uncertain whether it is better to use the sequence from ATOM records
    or the SEQRES (or the mmCIF equivalent, ``_pdbx_poly_seq``) field. The
    advantage of using ATOM records is that there are no surprises when
    building models from an alignment: a long hit with high coverage and score
    won't be built on top of a template with half its residues missing, giving
    a misleadingly-scored model. On the other hand, users may be surprised to
    find that their sequence doesn't confidently hit anything in the fold
    library when they know part of it has a determined structure, if the PDB
    entry it ought to match is missing many residues. In most cases, the effect
    is probably subtle, but this should be benchmarked.

The fold library consists of individual chains, stored as PDB files; profile
data, stored as A3M alignments and HMMs; and an SQLite databases containing
metadata.

Individual chains
=================

The official file format of the PDB is now the mmCIF format, so all our
coordinates for the fold library are parsed from mmCIF files. Reading only PDB
files would dramatically reduce our coverage, as many structures are only
distributed as mmCIF files, but most tools still only read PDB files. We solve
this by parsing mmCIF files but storing each chain as a PDB file. Since PDB
files are only limited by the maximum number of atoms or chains, this won't
cause a problem unless the structures of some truly gigantic protein chains are
solved.

The bulk of the work to convert mmCIF files into single-chain PDB files is done
by :py:class:`phyre_engine.component.db.db.ChainPDBBuilder` and
:py:class:`phyre_engine.tools.template.Template`. The `ChainPDBBuilder` is
responsible for parsing mmCIF files, generating the correct directory
structure, and optionally re-using existing chain files instead of overwriting
the old ones. The `Template` class is PhyreEngine's representation of a
template, and contains the "canonical" sequence of the template, the atomic
coordinates, and a mapping to the original residue IDs before renumbering. See
:py:class:`phyre_engine.tools.template.Template` for more information.

Selecting conformations
-----------------------

PhyreEngine will always use the first model in an mmCIF file, but alternate
protein conformations may be specified within a single MODEL record using the
"alternate location" identifier. This 80-column kludge has persisted into the
world of mmCIF, so PhyreEngine needs to cope gracefully with alternate
locations.

The :py:class:`phyre_engine.component.db.db.ChainPDBBuilder` component accepts
a `conf_sel` argument that is used to filter residues and atoms with alternate
locations. Each element in the `conf_sel` list should implement
:py:class:`phyre_engine.tools.conformation.ConformationSelector`, and is used
to select a single conformation. The default selectors for
:py:class:`~phyre_engine.component.db.db.ChainPDBBuilder` are
:py:class:`~phyre_engine.tools.conformation.PopulationMutationSelector`
:py:class:`~phyre_engine.tools.conformation.PopulationMicroHetSelector`, which
both inherit from
:py:class:`~phyre_engine.tools.conformation.PopulationConformationSelector` and
are used to select between different point mutations ("Mutation") and different
locations for atoms within the same residue ("MicroHet" for microheterogenity).
See the documentation for the individual classes, beginning with
:py:class:`~phyre_engine.tools.conformation.PopulationConformationSelector`,
for more details.

Profiles
========

A hidden Markov model (HMM) is generated for each unique sequence in the PDB.
Each HMM is stored as a :file:`.hhm` (sic) file, and is generated from a
multiple sequence alignment in A3M format. The A3M file is kept along with the
HMM for future reference. A file containing a 219-letter-alphabet sequence
representation for each sequence is also generated and stored in a
:file:`.cs219` file, which should allow the creation of databases that can be
searched by hhblits rather than the slower hhsearch.

.. warning::

    In practice, I have not had much luck using hhblits to search a fold
    library. Feel free to give it a try, but using the slower hhsearch is a
    safer option.

The component responsible for generating these profiles is
:py:class:`phyre_engine.component.foldlib.BuildProfiles`, which internally
calls several other components to go through the entire profile-generation
process. Profiles are stored underneath the `profile_dir` subdirectory passed
as a parameter to the :py:class:`~phyre_engine.component.foldlib.BuildProfiles`
component (or extracted from the ``foldlib.profile_dir`` configuration value),
organised in the same way as PDB files. Fasta files containing the canonical
sequence are also stored for the sake of convenience:

.. code-block:: none

    profile_dir
    ├── a3m
    │   ├── 15
    │   │   ├── 115l_A.a3m
    │   │   └── 215l_A.a3m
    │   ├── 16
    │   │   ├── 216l_A.a3m
    │   │   ├── 216l_B.a3m
    │   │   └── 316d_C.a3m
    │   └ ...
    ├── cs219
    │   └ ...
    ├── fasta
    │   └ ...
    └── hhm
        └ ...

FFindex / FFdata databases
--------------------------

Modern versions of hhsuite expect their databases to be in ffindex/ffdata
format, which consist of six files: three pairs of ffindex and ffdata files,
one pair each for :file:`a3m`, :file:`hhm` and :file:`cs219` data. In each of
these files, the ffdata file contains all constituent files---all hhm, a3m or
cs219 files---concatenated and separated by a nul (``\0``, ASCII 0) character.
The ffindex file contains a sorted list of file names, and the byte offset of
that entry in the ffdata file. This allows for items to be looked up by a fast
binary search within the ffindex file followed by a single seek in the ffdata
file.

The :py:class:`phyre_engine.component.hhsuite.BuildDatabase` component can be
used to build an ffindex/ffdata database suitable for use with hhblits. When
updating a fold library, the
:py:class:`phyre_engine.component.hhsuite.FFDatabaseUnlink` component should be
used to first remove any existing entries from the database to prevent
duplicate entries.

.. note::

    "Unlinking" an entry from an ffindex/ffdata database simply removes it from
    the ffindex file without touching the ffdata file. If you frequently update
    an ffindex/ffdata database by repeatedly adding and unlinking entries, the
    ffdata file will continue to fill with useless, unreachable data.

SQL template database
=====================

Data describing the PDB entries and individual chains that make up the fold
library are stored in an SQLite database, which is just a single file on the
filesystem. This gives us the speed and flexibility of an SQL database,
allowing fast lookup of PDB- and template-related data, but avoids the overhead
of an SQL server. The downside is that SQLite does not support concurrent
writes, so updating the fold library when using multiple machines in parallel
can be a bit tricky.

The format of the SQL portion of the template database is described in detail
in the :py:mod:`phyre_engine.tools.template` module documentation.
