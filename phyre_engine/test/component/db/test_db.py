import json
import tempfile
import unittest
import phyre_engine.component.db.db as db
import phyre_engine.tools.pdb as pdb
import phyre_engine.test
import phyre_engine.test.data
import phyre_engine.test.tools.test_pdb
from pathlib import Path
import shutil
import textwrap
import Bio.SeqIO
import Bio.PDB.PDBParser

@phyre_engine.test.requireFields("net_tests")
class TestStructureRetriever(unittest.TestCase):
    """Test STructureRetriever component."""
    _PIPE_STATE = {"PDB": "12as"}

    def test_retrieve(self):
        """Try and download some files."""
        types = ("pdb", "cif", db.StructureType.PDB, db.StructureType.MMCIF)

        for struc_type in types:
            with self.subTest("Getting structure", type=struc_type):
                with tempfile.TemporaryDirectory() as tmpdir:
                    retriever = db.StructureRetriever(struc_type, tmpdir)
                    retriever.run(self._PIPE_STATE)

                    if isinstance(struc_type, str):
                        suffix = struc_type
                    else:
                        suffix = struc_type.value

                    pdb_12as = Path(tmpdir, "2a/12as.{}.gz".format(suffix))
                    self.assertTrue(
                        pdb_12as.exists(),
                        "{!s} should exist".format(pdb_12as))


class TestChainPDBBuilder(unittest.TestCase):
    """Test ChainPDBBuilder class"""

    _MINIMAL_MMCIF = """\
loop_
_atom_type.symbol
C
N
O
P
S
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   2    C CA    . ALA A 1 4   ? 12.501 39.048 28.539  1.00 30.68 ? 4   ALA A CA    1
ATOM   7    C CA    . TYR A 1 5   A 15.552 39.410 26.282  1.00 8.51  ? 4   TYR A CA    1
ATOM   19   C CA    . ILE A 1 6   ? 17.791 39.281 29.375  1.00 23.27 ? 6   ILE A CA    1
ATOM   27   C CA    . ALA A 1 7   ? 16.004 36.186 30.742  1.00 5.20  ? 7   ALA A CA    1
ATOM   32   C CA    . LYS A 1 8   ? 16.137 34.313 27.425  1.00 4.96  ? 8   LYS A CA    1
ATOM   41   C CA    . GLN A 1 9   ? 19.794 35.327 26.885  1.00 8.15  ? 9   GLN A CA    1
ATOM   50   C CB    . ARG A 1 10  ? 20.706 33.656 30.141  1.00 19.12 ? 10  ARG A CA    1
HETATM 5128 P P     . AMP D 3 .   ? 24.511 18.911 8.472   1.00 27.61 ? 332 AMP A P     1
"""

    _MINIMAL_ORIG_MAPPING = [
        [' ', 4, ' '],
        [' ', 4, 'A'],
        [' ', 6, ' '],
        [' ', 7, ' '],
        [' ', 8, ' '],
        [' ', 9, ' '],
        [' ', 10, ' '],
        ['H_AMP', 332, ' ']
    ]
    _MINIMAL_SEQ = "AYIAKQ"
    _MINIMAL_SEQ_INDICES = [1, 2, 3, 4, 5, 6]

    def setUp(self):
        """Create a temporary directory for tests to use."""
        self.tmpdir = Path(tempfile.mkdtemp("-test", "chain-"))
        self.mmcif_dir = self.tmpdir / "mmcif"
        self.pdb_dir = self.tmpdir / "pdb"

        (self.mmcif_dir / "2a").mkdir(parents=True)
        self.pdb_dir.mkdir()

        with (self.mmcif_dir / "2a" / "12as.cif").open("w") as mmcif_out:
            mmcif_out.write(self._MINIMAL_MMCIF)

    def _build(self):
        builder = db.ChainPDBBuilder(str(self.mmcif_dir), str(self.pdb_dir))
        results = builder.run({"PDB": "12as", "chain": "A"})
        return builder, results

    def _get_remark(self, results, num):
        """Get a remark from some results."""
        pdb_path = results["structure"]
        with open(pdb_path, "r") as pdb_in:
            mapping_str = "".join(pdb.read_remark(pdb_in, num))
            return mapping_str

    def tearDown(self):
        """Remove temporary directory."""
        #shutil.rmtree(str(self.tmpdir))

    def test_create_chain(self):
        """Create a PDB file from a minimal MMCIF."""
        _, results = self._build()
        pdb_12as_A_path = Path(results["structure"])
        self.assertTrue(
            pdb_12as_A_path.exists(),
            "PDB file containing chain was created.")

    def test_renumbered_mapping(self):
        """Ensure that a mapping of renumbered -> original AAs was written."""
        _, results = self._build()
        mapping_str = self._get_remark(
            results,
            db.ChainPDBBuilder.ORIG_MAPPING_REMARK_NUM)
        mapping = json.loads(mapping_str)
        self.assertListEqual(mapping, self._MINIMAL_ORIG_MAPPING)

    def test_canonical_seq_indices(self):
        """Canonical index mapping should point to standard AAs with CAs."""
        _, results = self._build()
        mapping_str = self._get_remark(
            results,
            db.ChainPDBBuilder.CANONICAL_INDICES_REMARK_NUM)
        mapping = json.loads(mapping_str)
        self.assertListEqual(mapping, self._MINIMAL_SEQ_INDICES)

    def test_canonical_seq(self):
        """Canonical seq should include only standard AAs with CA atoms."""
        _, results = self._build()
        mapping = self._get_remark(
            results,
            db.ChainPDBBuilder.CANONICAL_SEQ_REMARK_NUM)
        self.assertEqual(mapping, self._MINIMAL_SEQ)

class TestPDBSequence(unittest.TestCase):
    """Tests for the PDBSequence component."""

    _MINIMAL_TEMPLATE = "REMARK 150 AAG\n"

    def setUp(self):
        _, self.template_name = tempfile.mkstemp(".pdb", "template-", text=True)
        with open(self.template_name, "w") as template_fh:
            template_fh.write(self._MINIMAL_TEMPLATE)
            template_fh.flush()

    def test_sequence(self):
        """Sequence should match atom_seq function."""
        seq_parser = db.PDBSequence()
        results = seq_parser.run({ "structure": self.template_name })
        self.assertEqual(
            results["sequence"], "AAG",
            "Sequence read correctly")

class TestReduceExpand(unittest.TestCase):
    """Test the Reduce and Expand components."""

    _FULL_LIST = [
        {"a": 1, "b": 1, "x": 123},
        {"a": 1, "b": 1, "x": 456},

        {"a": 1, "b": 2, "x": 789},

        {"a": 2, "b": 1, "x": 321},
        {"a": 2, "b": 1, "x": 654},

        {"a": 2, "b": 2, "x": 987},
    ]

    _REDUCTION = [
        [
            {"a": 1, "b": 1, "x": 123},
            {"a": 1, "b": 1, "x": 456},
        ],
        [
            {"a": 1, "b": 2, "x": 789},
        ],
        [
            {"a": 2, "b": 1, "x": 321},
            {"a": 2, "b": 1, "x": 654},
        ],
        [
            {"a": 2, "b": 2, "x": 987},
        ]
    ]
    _REDUCED_LIST = [
        {"a": 1, "b": 1, "x": 123},
        {"a": 1, "b": 2, "x": 789},
        {"a": 2, "b": 1, "x": 321},
        {"a": 2, "b": 2, "x": 987},
    ]

    def test_reduce(self):
        """Reduce according to a tuple of fields."""
        reducer = db.Reduce("list", ("a", "b"))
        result = reducer.run({"foo": "bar", "list": self._FULL_LIST.copy()})
        self.assertEqual(result["foo"], "bar", "Untouched extra data")
        self.assertListEqual(result["list"], self._REDUCED_LIST)
        self.assertListEqual(result["reduction"], self._REDUCTION)

    def test_expand(self):
        """Expand according to a tuple of fields."""
        expander = db.Expand("list", ("a", "b"))
        result = expander.run({
            "list": self._REDUCED_LIST.copy(),
            "reduction": self._REDUCTION.copy()
        })
        self.assertListEqual(
            sorted(result["list"], key=lambda e: e["x"]),
            sorted(self._FULL_LIST, key=lambda e: e["x"]))

class TestAnnotateCATrace(unittest.TestCase):
    """Test AnnotateCATrace component."""

    _CA_TRACE = textwrap.dedent("""\
        REMARK 150 AA
        REMARK 156 [1, 2]
        REMARK 161 [[" ", 4, " "], [" ", 5, " "], [" ", 6, " "]]
        ATOM      2  CA  ALA A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      7  CA  ALA A   2      15.552  39.410  26.282  1.00  8.51
        ATOM     32  CB  ALA A   3      16.137  34.313  27.425  1.00  4.96
    """)

    _NON_CA_TRACE = textwrap.dedent("""\
        REMARK 150 AA
        REMARK 156 [1, 2]
        REMARK 161 [[" ", 4, " "], [" ", 5, " "], [" ", 6, " "]]
        ATOM      2  CA  ALA A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      5  CB  ALA A   1      12.902  39.919  29.730  1.00 16.77
        ATOM      7  CA  ALA A   2      15.552  39.410  26.282  1.00  8.51
        ATOM     32  CB  ALA A   3      16.137  34.313  27.425  1.00  4.96
    """)

    def setUp(self):
        """Write mock structures to PDB files."""
        self.template_dir = Path(tempfile.mkdtemp())
        self.ca_pdb = (self.template_dir / "ca_trace.pdb")
        self.bb_pdb = (self.template_dir / "backbone.pdb")

        with self.ca_pdb.open("w") as pdb_out:
            pdb_out.write(self._CA_TRACE)

        with self.bb_pdb.open("w") as pdb_out:
            pdb_out.write(self._NON_CA_TRACE)

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(str(self.template_dir))

    def test_annotate_ca_trace(self):
        """Check that a CA trace is marked with ``ca_trace = True``."""
        annotator = db.AnnotateCATrace()
        result = annotator.run({"structure": str(self.ca_pdb)})
        self.assertTrue(result["ca_trace"], "Annotated as CA trace")

    def test_annotate_non_ca_trace(self):
        """Check that a PDB with backbone is marked ``ca_trace = False.``"""
        annotator = db.AnnotateCATrace()
        result = annotator.run({"structure": str(self.bb_pdb)})
        self.assertFalse(result["ca_trace"], "Annotated as not a CA trace")


if __name__ == "__main__":
    unittest.main()
