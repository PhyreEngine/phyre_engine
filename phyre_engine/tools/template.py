"""
This module contains tools for dealing with PhyreEngine's concept of template
structures.

At their most basic, a template structure is a protein structure onto which an
aligned sequence may be threaded in order to generate a homology model.
PhyreEngine's idea of a template is a little more specific: a template is a PDB
file containing a single chain, starting at residue 1 and with residues numbered
consecutively, with no alternate conformations and with some ``REMARK`` fields
indicating the residues containined within the template file that should be
included in the sequence of the template.

Each template must contain a ``REMARK 161`` (the sum of the decimal ascii codes
"T" and "M") field. The ``REMARK 161`` fields will contain a JSON-formatted list
of original author-assigned IDs in Biopython format. This will look like this:

.. code-block:: none

    REMARK 161 [
    REMARK 161     [' ', 3, ' '],
    REMARK 161     [' ', 3, 'A'],
    ...

And so on. Template residue IDs are given in three parts: a hetero flag, the
residue ID and an insertion code. In this case, residue 1 in the chain PDB is
mapped to residue 3 in the original template, and residue 2 maps to residue 3A.

Templates must also contain metadata giving the "canonical" sequence of a PDB
structure. This sequence is generated from all residues in the ``ATOM`` records
of the PDB file that are standard amino acids and have a ``CA`` atom. This
sequence is stored as a single line in ``REMARK 150`` (ASCII "C" + "S" for
Canonical Sequence).

Finally, a ``REMARK 156`` ("S" + "I" for Sequence Index) must be included
containing the *renumbered* residue ID corresponding to the canonical sequence.

For example, we might start with the following original PDB file (with
everything right of the residue index excised for legibility):

.. code-block:: none

    ATOM      1 CA A  ALA A 10
    ATOM      2 CA A  GLY A 11
    HETATM    3 CA A  AMP A 11A
    ATOM      4 CB A  ALA A 12

This PDB file will be renumbered beginning from 1, and insertion codes will be
stripped. With ``REMARK 161`` giving the original author-assigned residue
indices, the template will then look like this:

.. code-block:: none

    REMARK 161 [
    REMARK 161     [' ', 10, ' '],
    REMARK 161     [' ', 11, ' '],
    REMARK 161     ['H_AMP', 11, 'A'],
    REMARK 161     [' ', 12, ' ']
    REMARK 161 ]
    ATOM      1 CA A  ALA A 1
    ATOM      2 CA A  GLY A 2
    HETATM    3 CA A  AMP A 3
    ATOM      4 CB A  ALA A 4

Next, the canonical sequence and sequence mapping must be added. Residues 3 and
4 will not be included in the sequence: residue 3 is a ``HETATM`` and residue 4
does not contain a ``CA`` atom. The mapping is onto the
*renumbered* residue IDs:

.. code-block:: none

    REMARK 161 (unchanged)
    REMARK 150 AG
    REMARK 156 [1, 2]
    ATOM      1 CA A  ALA A 1
    ATOM      2 CA A  GLY A 2
    HETATM    3 CA A  AMP A 3
    ATOM      4 CB A  ALA A 4

.. warning::

    Templates may include JSON-encoded remarks in any format, not necessarily as
    it is shown in these examples. At the moment, JSON data is actually written
    on a single line. I hope that this doesn't break too many programs.
"""
import collections
import json
import phyre_engine.tools.pdb as pdb
import phyre_engine.tools.util as util
import Bio.PDB

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

    def __init__(self, chain, mapping, canonical_seq, canonical_indices):
        self.chain = chain
        self.mapping = mapping
        self.canonical_seq = canonical_seq
        self.canonical_indices = canonical_indices
        self.remarks = collections.defaultdict(lambda: [])

    @classmethod
    def build(cls, chain):
        """
        Built a template from a raw chain.

        :param Bio.PDB.Chain.Chain chain: Raw chain from which to generate a
            template.
        """
        mapping, chain = pdb.renumber(chain, "A")
        canonical_seq, canonical_res = pdb.atom_seq(chain)
        canonical_indices = [ r.get_id()[1] for r in canonical_res ]
        return cls(chain, mapping, canonical_seq, canonical_indices)

    @classmethod
    def load(cls, file):
        """
        Load a template from a stream, file or :py:class:`pathlib.Path` object.
        """
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        with util.Stream(file, "r") as stream:
            canonical_seq = "".join(
                pdb.read_remark(stream, cls.CANONICAL_SEQ_REMARK_NUM))
            stream.seek(0)

            index_remark = "".join(
                pdb.read_remark(stream, cls.CANONICAL_INDICES_REMARK_NUM))
            stream.seek(0)
            mapping_remark = "".join(
                pdb.read_remark(stream, cls.ORIG_MAPPING_REMARK_NUM))
            stream.seek(0)

            canonical_indices = json.loads(index_remark)
            # Load residue mappings, converting lists to tuples (not a JSON
            # concept).
            mapping = [tuple(res_id) for res_id in json.loads(mapping_remark)]

            structure = pdb_parser.get_structure("template", stream)
            chain = structure[0]["A"]
            return cls(chain, mapping, canonical_seq, canonical_indices)

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
        pdb.write_remark(
            pdb_out, [json.dumps(self.mapping)],
            self.ORIG_MAPPING_REMARK_NUM)

        for remark_num, remark_lines in self.remarks.items():
            pdb.write_remark(pdb_out, remark_lines, remark_num)
        pdbio = Bio.PDB.PDBIO()
        pdbio.set_structure(self.chain)
        pdbio.save(pdb_out)
