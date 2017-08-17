import json
import tempfile
import unittest
import phyre_engine.component.db.db as db
import phyre_engine.test
import phyre_engine.test.data
import phyre_engine.test.tools.test_pdb
from pathlib import Path
import shutil
import Bio.SeqIO

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

    _12ASA_SEQ = (
        "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
        "AVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDE"
        "DRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSE"
        "EFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGG"
        "KLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMG"
        "IRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTML"
        "LLQLPHIGQVQAGVWPAAVRESVPSLLN")

    def setUp(self):
        """Create a temporary directory for tests to use."""
        self.tmpdir = Path(tempfile.mkdtemp("-test", "chain-"))
        self.mmcif_dir = Path(phyre_engine.test.data.__file__).parent / "mmcif"

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.tmpdir))

    def test_build(self):
        """
        Try and extract a chain. Test that the extracted chain exists and that
        it contains REMARK 149 and ATOM lines matching the required sequence.
        """
        pdb_dir = self.tmpdir

        builder = db.ChainPDBBuilder(str(self.mmcif_dir), str(pdb_dir))
        results = builder.run({"PDB": "12as", "chain": "A"})

        pdb_12as_A_path = Path(results["structure"])
        self.assertTrue(
            pdb_12as_A_path.exists(),
            "PDB file containing chain was created.")

        with pdb_12as_A_path.open("r") as pdb_in:
            self.assertEqual(
                str(Bio.SeqIO.read(pdb_in, "pdb-atom").seq), self._12ASA_SEQ,
                "Atom sequence matches read sequence")

        with pdb_12as_A_path.open("r") as pdb_in:
            json_str = ""
            for line in pdb_in:
                if not line.startswith("REMARK 149"):
                    continue
                json_str += line.replace("REMARK 149 ", "").rstrip()
            mapping = json.loads(json_str)
            map_dict = {m[0]: tuple(m[1]) for m in mapping}
            self.assertEqual(map_dict[1], (' ', 4, ' '))


    def test_build_with_selectors(self):
        """Write a structure with an active conformation selector."""
        class FilterHetatms:
            def select(self, chain):
                for res in list(chain.get_residues()):
                    if res.get_id()[0] != ' ':
                        del chain[res.get_id()]
                return chain

        pdb_dir = self.tmpdir

        builder = db.ChainPDBBuilder(
            str(self.mmcif_dir), str(pdb_dir),
            conf_sel=[FilterHetatms()])
        results = builder.run({"PDB": "12as", "chain": "A"})

        with Path(results["structure"]).open("r") as pdb_in:
            # Results should match the same sequence without the terminating
            # residue, which is marked as a HETATM.
            self.assertEqual(
                str(Bio.SeqIO.read(pdb_in, "pdb-atom").seq),
                self._12ASA_SEQ[0:-1],
                "Atom sequence matches read sequence")

class TestPDBSequence(unittest.TestCase):
    """Tests for the PDBSequence component."""

    def test_sequence(self):
        """Sequence should match atom_seq function."""

        with tempfile.NamedTemporaryFile("w") as pdb_fh:
            pdb_fh.write(phyre_engine.test.tools.test_pdb.ATOM_SEQ_TEST_PDB)
            pdb_fh.flush()
            seq_parser = db.PDBSequence()
            results = seq_parser.run({
                "structure": pdb_fh.name
            })
            self.assertEqual(
                results["sequence"], "AG",
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


if __name__ == "__main__":
    unittest.main()
