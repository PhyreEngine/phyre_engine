"""
This module contains tools for dealing with PhyreEngine's concept of template
structures.

Templates in PhyreEngine
------------------------

At their most basic, a template structure is a protein structure onto which an
aligned sequence may be threaded in order to generate a homology model.
PhyreEngine's idea of a template is a little more specific: a template is a PDB
file containing a single chain, starting at residue 1 and with residues numbered
consecutively, with no alternate conformations.

Templates contain metadata giving the "canonical" sequence of a PDB structure.
This sequence is generated from all residues in the ``ATOM`` records of the PDB
file that are standard amino acids and have a ``CA`` atom. Templates also
contain metadata giving the mapping between the canonical sequence and the
(renumbered) residue IDs of each residue contained within the canonical
sequence.

For example, we might start with the following original PDB file (with
everything right of the residue index excised for legibility):

.. code-block:: none

    ATOM      1 CA A  ALA A 10
    ATOM      2 CA A  GLY A 11
    HETATM    3 CA A  AMP A 11A
    ATOM      4 CB A  ALA A 12

This PDB file will be renumbered beginning from 1, and insertion codes will be
stripped. It will then look like this:

.. code-block:: none

    ATOM      1 CA A  ALA A 1
    ATOM      2 CA A  GLY A 2
    HETATM    3 CA A  AMP A 3
    ATOM      4 CB A  ALA A 4

The original residue mapping will be preserved as an array of residue IDs.
BioPython's idea of residue IDs is used, so each residue ID is a triplet of
values. The first value in the triplet indicates the type of ``HETATM``, if
any, the second gives the residue ID and the third gives the insertion code, if
any. Null values are represented as spaces. In this case, the mapping will look
like this:

.. code-block::

    [
        (' ', 10, ' '),
        (' ', 11, ' '),
        ('H_AMP', 11, 'A'),
        (' ', 12, ' ')
    ]

Next, the canonical sequence and sequence mapping are calculated. In this
example, residues 3 and 4 will not be included in the sequence: residue 3 is a
``HETATM`` and residue 4 does not contain a ``CA`` atom. The canonical sequence
is ``AG``, and the the *renumbered* residue IDs are ``[1, 2]]``. For the sake
of convenience, the canonical sequence and list of canonical indices will be
stoerd in ``REMARK`` fields:

.. code-block:: none

    REMARK 150 AG
    REMARK 156 [1, 2]
    ATOM      1 CA A  ALA A 1
    ATOM      2 CA A  GLY A 2
    HETATM    3 CA A  AMP A 3
    ATOM      4 CB A  ALA A 4

The canonical sequence is stored in ``REMARK 150``, named for the sum of the
ASCII values of ``C`` and ``S`` (Canonical Sequence). The residue IDs of the
residues in the canonical sequence are stored as a JSON array in ``REMARK 156``
(``S`` + ``I`` for Sequence Index).

.. warning::

    Templates may include JSON-encoded remarks in any format, not necessarily as
    it is shown in these examples. At the moment, JSON data is actually written
    on a single line. I hope that this doesn't break too many programs.

Template Objects
----------------

A template is represented by a :py:class:`~.Template` object. Each
:py:class:`phyre_engine.tools.template.Template` object has a with a
:py:class:`~Bio.PDB.Chain.Chain` object, canonical sequence and list of
canonical indices. Each object also contains an array of residue IDs giving the
ID of the original residue before renumbering and sanitisation.

.. seealso::

    :py:class:`phyre_engine.tools.template.Template`
        For more information about building and instantiating templates.

Template Database
-----------------

Template data are stored in a template database. The template database consists
of two components: an SQLite database to contain template metadata, and a
directory containing the PDB files for each template. PDB files are placed into
subdirectories named for the middle two letters of the PDB ID, then the full
PDB ID, and then into a file named for the PDB ID and chain ID. For example, if
the root is at :file:`{foldlib}`, then chain ``A`` of ``1XYZ`` will be stored
in :file:`{foldlib}/xy/1xyz/1xyz_A.pdb`.

The SQLite database portion of the template database contains metadata for an
entire PDB file such as deposition date, as well as metadata for the template
such as canonical sequence.

The metadata describing a PDB entry are the same as the metadata that can be
specified when adding a PDB entry via :py:meth:`~.TemplateDatabase.add_pdb`.
Metadata may be retrieved via the :py:meth:`~.TemplateDatabase.get_pdb`
method.

The metadata describing a template is the same as described in
:py:class:`~.Template`.

"""
import collections
import json
import operator
from pathlib import Path
import sqlite3

import phyre_engine.tools.pdb as pdb
import phyre_engine.tools.util as util
import Bio.PDB

class TemplateDatabase:
    """
    Database of templates. This is implemented internally as a relational
    (sqlite) database.

    :param str database: Path to the database file.

    :param str file_root: Root of the directory structure containing the
        template files.

    :param bool exclusive: If `True`, immediately begin an ``EXCLUSIVE``
        transaction, preventing any access to the database.

    :param callable trace: Callback used for printing SQL traces.

    .. warning::

        Changes will not be committed to the database until :py:meth:`~.commit`
        is called.
    """

    CREATE = """\
        CREATE TABLE meta (
            updated          DATE
        );
        INSERT INTO meta (updated) VALUES (NULL);

        CREATE TABLE pdbs (
            pdb_id           TEXT,
            deposition_date  DATE NOT NULL,
            last_update_date DATE NOT NULL,
            release_date     DATE NOT NULL,
            method           TEXT NOT NULL,
            resolution       REAL,
            organism_name    TEXT,
            organism_id      INTEGER,
            title            TEXT,
            descriptor       TEXT,
            PRIMARY KEY(pdb_id)
        );

        CREATE TABLE chains (
            pdb_id             TEXT,
            chain_id           TEXT,
            canonical_sequence TEXT NOT NULL,
            PRIMARY KEY (pdb_id, chain_id),
            FOREIGN KEY (pdb_id) REFERENCES pdbs(pdb_id)
                 ON DELETE CASCADE
        );

        CREATE TABLE sequence_reps (
            pdb_id             TEXT,
            chain_id           TEXT,
            canonical_sequence TEXT,
            PRIMARY KEY (pdb_id, chain_id, canonical_sequence),
            FOREIGN KEY (pdb_id) REFERENCES pdbs(pdb_id),
            FOREIGN KEY (pdb_id, chain_id) REFERENCES chains(pdb_id, chain_id)
                 ON DELETE CASCADE
        );

        CREATE TABLE canonical (
            pdb_id           TEXT,
            chain_id         TEXT,
            --
            -- Beginning from zero, gives the index in the canonical sequence.
            sequence_index   INTEGER,
            --
            -- Single-letter residue type.
            aa               TEXT,
            --
            -- Residue ID (beginning from 1) of corresponding residue in the
            -- PDB file. This is the ID of a residue after renumbering.
            residue_id       INTEGER,
            PRIMARY KEY (pdb_id, chain_id, sequence_index),
            FOREIGN KEY (pdb_id, chain_id) REFERENCES chains(pdb_id, chain_id)
                 ON DELETE CASCADE
        );

        CREATE TABLE original_residues (
            pdb_id           TEXT,
            chain_id         TEXT,
            -- Beginning from 0, the index of the residue as it appears in the
            -- original chain.
            sequence_index    INTEGER,
            --
            -- Hetero flag, original residue ID and insertion code. These fully
            -- specify a residue, and can be used with BioPython to look up a
            -- residue. These give the original residue ID before renumbering.
            hetero_flag      TEXT,
            orig_residue_id  INTEGER,
            insertion_code   TEXT,
            PRIMARY KEY (pdb_id, chain_id, sequence_index),
            FOREIGN KEY (pdb_id, chain_id) REFERENCES chains(pdb_id, chain_id)
                 ON DELETE CASCADE
        );

        CREATE INDEX seq_index chains(canonical_sequence);
        """

    CREATE_INDICES = """
    CREATE INDEX seq_index chains(canonical_sequence);
    """

    DROP_INDICES = """
    DROP INDEX seq_index;
    """

    SELECT_DB_META = """
    SELECT * FROM meta
    """

    UPDATE_DB_META = """
    UPDATE meta SET updated = :updated
    """

    SELECT_PDB = """
    SELECT * FROM pdbs WHERE pdb_id = :pdb_id
    """

    SELECT_TEMPLATE = """
    SELECT canonical_sequence
      FROM chains
     WHERE pdb_id = :pdb_id AND chain_id = :chain_id
    """

    SELECT_TEMPLATE_CANONICAL = """
    SELECT residue_id, aa
      FROM canonical
     WHERE pdb_id = :pdb_id AND chain_id = :chain_id
     ORDER BY sequence_index ASC
    """

    SELECT_TEMPLATE_ORIGINAL = """
    SELECT IFNULL(hetero_flag, ' ') as hetero_flag,
           orig_residue_id,
           IFNULL(insertion_code, ' ') as insertion_code
      FROM original_residues
     WHERE pdb_id = :pdb_id AND chain_id = :chain_id
     ORDER BY sequence_index ASC
    """

    SELECT_ALL_REP_SEQS = """
    SELECT canonical_sequence FROM chains GROUP BY canonical_sequence
    """

    SELECT_ALL_REPS = """
    SELECT pdb_id, chain_id FROM sequence_reps
    """

    SELECT_FIND_REP = """
    SELECT p.pdb_id, c.chain_id
      FROM chains c
     INNER JOIN pdbs p
        ON c.pdb_id = p.pdb_id
     WHERE c.canonical_sequence = :canonical_sequence
     ORDER BY p.resolution IS NULL, p.resolution, c.chain_id
     LIMIT 1
    """

    SELECT_SEQUENCE_CLUSTER = """
    SELECT c2.pdb_id, c2.chain_id
      FROM chains c1
     INNER JOIN chains c2
        ON c2.canonical_sequence = c1.canonical_sequence
     WHERE (c1.pdb_id = LOWER(:pdb_id) AND c1.chain_id = :chain_id)
       AND NOT (c2.pdb_id = LOWER(:pdb_id) AND c2.chain_id = :chain_id)
    """

    DELETE_ALL_SEQ_REPS = """
    DELETE FROM sequence_reps
    """

    DELETE_SEQ_REP_BY_SEQ = """
    DELETE FROM sequence_reps
     WHERE canonical_sequence = :canonical_sequence
    """

    INSERT_SEQ_REP = """
    INSERT INTO sequence_reps (pdb_id, chain_id, canonical_sequence)
    VALUES (:pdb_id, :chain_id, :canonical_sequence)
    """

    INSERT_PDB = """
    INSERT INTO pdbs (
           pdb_id, deposition_date, last_update_date, release_date,
           method, resolution,
           organism_name, organism_id,
           title, descriptor)
    VALUES (
           LOWER(:pdb_id), :deposition_date, :last_update_date, :release_date,
           :method, :resolution,
           :organism_name, :organism_id,
           :title, :descriptor)
    """

    UPDATE_PDB = """
    UPDATE pdbs
       SET deposition_date  = :deposition_date,
           last_update_date = :last_update_date,
           release_date     = :release_date,
           method           = :method,
           resolution       = :resolution,
           organism_name    = :organism_name,
           organism_id      = :organism_id,
           title            = :title,
           descriptor       = :descriptor
     WHERE pdb_id = LOWER(:pdb_id)
    """

    INSERT_TEMPLATE = """
    INSERT INTO chains (pdb_id, chain_id, canonical_sequence)
    VALUES (LOWER(:pdb_id), :chain_id, :canonical_sequence)
    """

    INSERT_CANONICAL = """
    INSERT INTO canonical (
           pdb_id, chain_id, sequence_index, aa, residue_id)
    VALUES (LOWER(:pdb_id), :chain_id, :sequence_index, :aa, :residue_id)
    """

    INSERT_MAPPING = """
    INSERT INTO original_residues (
           pdb_id, chain_id, sequence_index,
           hetero_flag, orig_residue_id, insertion_code)
    VALUES (
           LOWER(:pdb_id), :chain_id, :sequence_index,
           :hetero_flag, :orig_residue_id, :insertion_code)
    """

    DELETE_PDB = """
    DELETE FROM pdbs
     WHERE pdb_id = :pdb_id
    """

    DELETE_TEMPLATE = """
    DELETE FROM chains
     WHERE pdb_id = :pdb_id AND chain_id = :chain_id
    """

    class TemplateNotFoundException(Exception):
        ERR_MSG = "Template {} not found in database."

        def __init__(self, pdb_id, chain_id):
            super().__init__(self.ERR_MSG.format((pdb_id, chain_id)))

    class PdbNotFoundException(Exception):
        ERR_MSG = "PDB entry {} not found in database."

        def __init__(self, pdb_id):
            super().__init__(self.ERR_MSG.format(pdb_id))


    def __init__(self, database, file_root, exclusive=False, trace=None):
        self.database = database
        self.conn = sqlite3.connect(
            self.database,
            detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
        self.conn.row_factory = sqlite3.Row

        if trace is not None:
            self.conn.set_trace_callback(trace)

        self.conn.execute("PRAGMA foreign_keys = 1")
        if exclusive:
            self.conn.execute("BEGIN EXCLUSIVE TRANSACTION")
        self.file_root = Path(file_root)

    @property
    def updated(self):
        """Date of last database update."""
        return self.conn.execute(self.SELECT_DB_META).fetchone()["updated"]

    @updated.setter
    def updated(self, updated):
        """Set date of last update."""
        self.conn.execute(self.UPDATE_DB_META, {"updated": updated})
        self.conn.commit()

    @classmethod
    def create(cls, database):
        """Create an empty template database."""
        conn = sqlite3.connect(database)
        conn.executescript(cls.CREATE)
        conn.commit()

    def _canonical(self, pdb_id, chain_id):
        """Returns the canonical sequence and canonical indices."""
        where = {"pdb_id": pdb_id.lower(), "chain_id": chain_id}
        canonical = self.conn.execute(
            self.SELECT_TEMPLATE_CANONICAL, where
        ).fetchall()
        canonical_sequence = "".join([r["aa"] for r in canonical])
        canonical_indices = [r["residue_id"] for r in canonical]
        return canonical_sequence, canonical_indices

    def _original(self, pdb_id, chain_id):
        """Returns the IDs of the original residues."""
        where = {"pdb_id": pdb_id.lower(), "chain_id": chain_id}
        original_residues = self.conn.execute(
            self.SELECT_TEMPLATE_ORIGINAL, where
        ).fetchall()
        mapping = [
            (r["hetero_flag"], r["orig_residue_id"], r["insertion_code"])
            for r in original_residues]
        return mapping

    def sequence_reps(self):
        """
        Return a list of sequence representatives as `(pdb_id, chain_id)`
        tuples.
        """
        rows = self.conn.execute(self.SELECT_ALL_REPS).fetchall()
        return [(row["pdb_id"], row["chain_id"]) for row in rows]

    def get_pdb(self, pdb_id):
        """
        Retrieve metadata for the PDB with ID `pdb_id`.

        :returns: Dictionary containing metadata in the form described
            in :py:meth:`.add_pdb`.
        """
        result_row = self.conn.execute(
            self.SELECT_PDB,
            {"pdb_id": pdb_id.lower()}).fetchone()
        if result_row is None:
            raise self.PdbNotFoundException(pdb_id)

        result_dict = dict(result_row)
        del result_dict["pdb_id"]
        return result_dict

    def get_template(self, pdb_id, chain_id):
        """
        Returns the :py:class:`.Template` with the specified PDB and chain
        IDs.

        :param str pdb_id: PDB identifier for the structure.
        :param str chain_id: Identifier of the chain.
        """
        where = {"pdb_id": pdb_id.lower(), "chain_id": chain_id}
        template_results = self.conn.execute(self.SELECT_TEMPLATE,
                                             where).fetchone()
        if template_results is None:
            raise self.TemplateNotFoundException(pdb_id, chain_id)

        canon_seq, canon_idx = self._canonical(pdb_id, chain_id)
        mapping = self._original(pdb_id, chain_id)

        template_file = pdb.find_pdb(pdb_id, chain_id, self.file_root)
        if template_file is None:
            raise FileNotFoundError("Could not find template {} in {}".format(
                (pdb_id, chain_id), self.file_root))

        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        with pdb.open_pdb(template_file) as pdb_in:
            chain = pdb_parser.get_structure("", pdb_in)[0]["A"]

        return Template(pdb_id, chain_id, chain, mapping, canon_seq, canon_idx)

    def add_template(self, template):
        """
        Add template to database, including residue mappings.

        The corresponding PDB entry must already have been inserted. If it has
        not, a foreign key constraint error will be raised.

        :param str pdb_id: PDB identifier for the structure.
        :param str chain_id: Identifier of the chain.
        :param Template template: Template to add.
        """
        self.conn.execute(
            self.INSERT_TEMPLATE, {
                "pdb_id": template.pdb_id.lower(),
                "chain_id": template.chain_id,
                "canonical_sequence": template.canonical_seq})

        iterator = zip(template.canonical_seq, template.canonical_indices)
        for i, (aa, seq_idx) in enumerate(iterator):
            fields = {
                "pdb_id": template.pdb_id.lower(),
                "chain_id": template.chain_id,
                "sequence_index": i,
                "aa": aa,
                "residue_id": seq_idx,
            }
            self.conn.execute(self.INSERT_CANONICAL, fields)

        for i, orig_id in enumerate(template.mapping):
            hetero_flag = None if orig_id[0] == " " else orig_id[0]
            icode = None if orig_id[2] == " " else orig_id[2]
            fields = {
                "pdb_id": template.pdb_id.lower(),
                "chain_id": template.chain_id,
                "sequence_index": i,
                "hetero_flag": hetero_flag,
                "insertion_code": icode,
                "orig_residue_id": orig_id[1],
            }
            self.conn.execute(self.INSERT_MAPPING, fields)

    @staticmethod
    def _format_metadata(pdb_id, metadata):
        """
        Format PDB ID and metadata into dict suitable for use filling
        placeholders in an `execute` statement.
        """
        fields = collections.defaultdict(lambda: None)
        fields["pdb_id"] = pdb_id.lower()
        fields.update(metadata)
        for date in ("deposition_date", "last_update_date", "release_date"):
            fields[date] = fields[date].strftime("%Y-%m-%d")
        return fields

    def add_pdb(self, pdb_id, metadata):
        """
        Insert entry for PDB with ID `pdb_id` into the database.

        The following fields must be present in the `metadata` dictionary:

        ``deposition_date``, ``last_update_date``, ``release_date``:
            :py:class:`datetime.date` objects giving the deposition date,
            date of the last revision and the relase date, respectively.

        ``method``
            String describing the method used to resolve the structure.  This
            is taken from the ``_exptl.method`` field of the mmCIF file.

        The following fields are optional, but should be present in most
        cases:

        ``resolution``
            Real number giving the resolution of the structure, if applicable.

        ``organism_name``, ``organism_id``
            String and integer giving the name and NCBI taxonomy ID of the
            organism in which the protein was found.

        ``title``, ``descriptor``
            Title and description (``_struct.title`` and
            ``_struct.pdbx_descriptor`` fields) of the structure.

        :param str pdb_id: PDB identifier.

        :param dict metadata: Dictionary of metadata to be inserted into the
            database.
        """
        fields = self._format_metadata(pdb_id, metadata)
        self.conn.execute(self.INSERT_PDB, fields)

    def update_pdb(self, pdb_id, metadata):
        """
        Update a PDB entry in the fold library database.

        See :py:class:`.add_pdb` for information about the fields expected in
        `metadata`.
        """
        fields = self._format_metadata(pdb_id, metadata)
        self.conn.execute(self.UPDATE_PDB, fields)

    def del_pdb(self, pdb_id):
        """
        Delete the PDB with ID `pdb_id` and all templates based on that PDB
        entry from the database.

        :param str pdb_id: Identifier of the PDB entry to delete.
        """
        self.conn.execute(self.DELETE_PDB, {"pdb_id": pdb_id.lower()})

    def del_template(self, pdb_id, chain_id):
        """
        The the template with the given PDB and chain IDs, along with all
        metatadata associated with that template.

        :param str pdb_id: PDB ID of the template to delete.
        :param str chain_id: Chain ID of the template to delete.
        """
        self.conn.execute(
            self.DELETE_TEMPLATE,
            {"pdb_id": pdb_id.lower(), "chain_id": chain_id})

    def update_all_seq_reps(self):
        """
        Update all sequence representatives.

        This will first clear all sequence representatives from the
        ``sequence_reps`` table. It will then find all unique sequences
        from the ``chains`` table, and loop over them to find the PDB
        and chain IDs of template with that sequence with the lowest
        resolution. Ties are broken by sorting on PDB ID, then chain ID.
        """
        self.conn.execute(self.DELETE_ALL_SEQ_REPS)
        seq_rows = self.conn.execute(self.SELECT_ALL_REP_SEQS).fetchall()
        for seq_row in seq_rows:
            metadata_row = self.conn.execute(
                self.SELECT_FIND_REP, seq_row).fetchone()
            rep_data = {
                "pdb_id": metadata_row["pdb_id"],
                "chain_id": metadata_row["chain_id"],
                "canonical_sequence": seq_row["canonical_sequence"],
            }
            self.conn.execute(self.INSERT_SEQ_REP, rep_data)

    def update_seq_rep(self, pdb_id, chain_id):
        """
        Update the sequence representative for all templates with the same
        sequence as the template with PDB ID `pdb_id` and chain ID `chain_id`.

        :param str pdb_id: PDB ID of the template to update.
        :param str chain_id: Chain ID of the template to update.
        """
        where = {"pdb_id": pdb_id.lower(), "chain_id": chain_id}
        template_row = self.conn.execute(
            self.SELECT_TEMPLATE,
            {"pdb_id": pdb_id.lower(), "chain_id": chain_id}
        ).fetchone()
        metadata_row = self.conn.execute(
            self.SELECT_FIND_REP,
            template_row
        ).fetchone()
        rep_data = {
            "pdb_id": metadata_row["pdb_id"],
            "chain_id": metadata_row["chain_id"],
            "canonical_sequence": template_row["canonical_sequence"],
        }
        self.conn.execute(self.DELETE_SEQ_REP_BY_SEQ, template_row)
        self.conn.execute(self.INSERT_SEQ_REP, rep_data)

    def expand_seq_reps(self, pdb_id, chain_id):
        """
        Return a list of all chains with canonical sequence matching the given
        chain.

        :param str pdb_id: PDB ID of a cluster member.
        :param str chain_id: Chain ID of cluster member.
        :returns: List of 2-tuples containing the PDB IDs and chain IDs of the
            members of this cluster, *excluding* the original member.
        :rtype: list[tuple(str, str)]
        """
        placeholders = {"pdb_id": pdb_id, "chain_id": chain_id}
        rows = self.conn.execute(self.SELECT_SEQUENCE_CLUSTER,
                                 placeholders).fetchall()
        return [(row["pdb_id"], row["chain_id"]) for row in rows]

    def begin(self, exclusive=False):
        """
        Begin a new transaction. If `exclusive` is `True`, the database is
        fully locked.
        """
        if exclusive:
            self.conn.execute("BEGIN EXCLUSIVE TRANSACTION")
        else:
            self.conn.execute("BEGIN TRANSACTION")

    def commit(self):
        """Commit changes to the database."""
        self.conn.commit()

    def rollback(self):
        """Roll back changes made in this transaction."""
        self.conn.rollback()

    def drop_indices(self):
        """Drop all indices."""
        self.conn.executescript(self.DROP_INDICES)

    def create_indices(self):
        """Create all indices."""
        self.conn.executescript(self.CREATE_INDICES)


class Template:
    """
    Represents a template, including all metadata described in
    :py:mod:`phyre_engine.tools.template`.

    Templates may be generated from a :py:class:`Bio.PDB.Chain.Chain` object
    via the :py:meth:`.build` class method, or parsed from a template file using
    the :py:meth:`.load` class method.

    :ivar chain: Renumbered chain of residues.
    :vartype chain: :py:class:`Bio.PDB.Chain.Chain`.

    :ivar mapping: List of author-defined residue IDs.
    :vartype mapping: list[tuple(str, int, str)]

    :ivar canonical_seq: String of single-letter amino acid codes giving the
        canonical sequence of this template.
    :vartype canonical_seq: str

    :ivar canonical_indices: Indices of the residues (beginning from 1) included
        in the canonical sequence.
    :vartype canonical_indices: list[int]

    :ivar remarks: Dictionary containing extra ``REMARK`` to write, indexed by
        ``REMARK`` number.
    :vartype remarks: :py:class:`collections.defaultdict` with default value of
        ``[]``.
    """

    #: Number used in ``REMARK`` fields for the JSON-encoded mapping between
    #: the current residue index and the residue ID assigned by the author of
    #: the original structure.
    ORIG_MAPPING_REMARK_NUM = 161

    #: ``REMARK`` number of the canonical sequence.
    CANONICAL_SEQ_REMARK_NUM = 150

    #: ``REMARK`` number of the canonical sequence residue IDs.
    CANONICAL_INDICES_REMARK_NUM = 156

    def __init__(self, pdb_id, chain_id, chain, mapping, canonical_seq,
                 canonical_indices):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.chain = chain
        self.mapping = mapping
        self.canonical_seq = canonical_seq
        self.canonical_indices = canonical_indices
        self.remarks = collections.defaultdict(list)

    @classmethod
    def build(cls, pdb_id, chain_id, chain):
        """
        Built a template from a raw chain.

        :param Bio.PDB.Chain.Chain chain: Raw chain from which to generate a
            template.
        """
        mapping, chain = pdb.renumber(chain, "A")
        canonical_seq, canonical_res = pdb.atom_seq(chain)
        canonical_indices = [ r.get_id()[1] for r in canonical_res ]
        return cls(pdb_id, chain_id,
                   chain, mapping, canonical_seq, canonical_indices)

    def write(self, pdb_out):
        """
        Write this template in PDB format, including metadata.

        :param pdb_out: Output stream.
        """
        pdb.write_remark(
            pdb_out, [self.canonical_seq],
            self.CANONICAL_SEQ_REMARK_NUM)
        pdb.write_remark(
            pdb_out, [json.dumps(self.canonical_indices)],
            self.CANONICAL_INDICES_REMARK_NUM)

        for remark_num, remark_lines in self.remarks.items():
            pdb.write_remark(pdb_out, remark_lines, remark_num)
        pdbio = Bio.PDB.PDBIO()
        pdbio.set_structure(self.chain)
        pdbio.save(pdb_out)
