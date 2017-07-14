import json
import os
import tempfile
import unittest
import phyre_engine.component.db.db as db
import phyre_engine.test.data
from pathlib import Path
import shutil
import Bio.SeqIO

class TestStructureRetriever(unittest.TestCase):
    """Test STructureRetriever component."""
    _PIPE_STATE = {"templates": [
        {"PDB": "12as"},
        {"PDB": "4HHB"}
    ]}

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
                    pdb_4hhb = Path(tmpdir, "hh/4hhb.{}.gz".format(suffix))
                    self.assertTrue(
                        pdb_12as.exists(),
                        "{!s} should exist".format(pdb_12as))
                    self.assertTrue(
                        pdb_4hhb.exists(),
                        "{!s} should exist".format(pdb_4hhb))

class TestFunctions(unittest.TestCase):
    """Test utility functions."""

    def setUp(self):
        self.mmcif_dir = Path(phyre_engine.test.data.__file__).parent / "mmcif"

    def test_pdb_path(self):
        """Test that our PDB path function works."""
        self.assertEqual(
            db.pdb_path("1ABC", ".pdb.gz"),
            Path("ab/1abc.pdb.gz"))

        self.assertEqual(
            db.pdb_path("1AbC", ".cif"),
            Path("ab/1abc.cif"))

        self.assertEqual(
            db.pdb_path("1abC", ".cif.gz", "A"),
            Path("ab/1abc/1abc_A.cif.gz"))

    def test_find_pdb(self):
        """Test "find_pdb" function to find a structure."""
        self.assertEqual(
            db.find_pdb("12as", base_dir=self.mmcif_dir),
            self.mmcif_dir / "2a/12as.cif")

        self.assertEqual(
            db.find_pdb("4n6v", base_dir=self.mmcif_dir),
            self.mmcif_dir / "n6/4n6v.cif.gz")

        self.assertIsNone(db.find_pdb("1abc"))

    def test_open_pdb_uncompressed(self):
        """Open an uncompressed mmcif file."""
        pdb_path = db.pdb_path("12as", ".cif", base_dir=self.mmcif_dir)
        with db.open_pdb(pdb_path) as pdb_in:
            self.assertGreater(len(pdb_in.readlines()), 0)

    def test_open_pdb_compressed(self):
        """Open a compressed mmcif file."""
        pdb_path = db.pdb_path("4n6v", ".cif.gz", base_dir=self.mmcif_dir)
        with db.open_pdb(pdb_path) as pdb_in:
            self.assertGreater(len(pdb_in.readlines()), 0)


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
        results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

        pdb_12as_A_path = Path(results["templates"][0]["structure"])
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
        results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

        with Path(results["templates"][0]["structure"]).open("r") as pdb_in:
            # Results should match the same sequence without the terminating
            # residue, which is marked as a HETATM.
            self.assertEqual(
                str(Bio.SeqIO.read(pdb_in, "pdb-atom").seq),
                self._12ASA_SEQ[0:-1],
                "Atom sequence matches read sequence")

if __name__ == "__main__":
    unittest.main()
