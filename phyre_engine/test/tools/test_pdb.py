"""Tests for the phyre_engine.tools.pdb module."""
import unittest
import phyre_engine.test.data
from pathlib import Path
import phyre_engine.tools.pdb as pdb
import io
import Bio.PDB
import random
import tempfile

ATOM_SEQ_TEST_PDB = (
    "HETATM    1  CA  FOO A   1      11.751  37.846  29.016  1.00 44.65\n"
    "ATOM      2  CA  ALA A   2      12.501  39.048  28.539  1.00 30.68\n"
    "ATOM      3  C   ALA A   2      13.740  38.628  27.754  1.00 24.74\n"
    "ATOM      4  CB  GLU A   3      14.207  37.495  27.890  1.00 25.59\n"
    "ATOM      5  CA  GLY A   4      14.207  37.495  27.890  1.00 25.59\n")

class TestFunctions(unittest.TestCase):
    """Test utility functions."""

    def setUp(self):
        self.mmcif_dir = Path(phyre_engine.test.data.__file__).parent / "mmcif"

    def test_pdb_path(self):
        """Test that our PDB path function works."""
        self.assertEqual(
            pdb.pdb_path("1ABC", ".pdb.gz"),
            Path("ab/1abc.pdb.gz"))

        self.assertEqual(
            pdb.pdb_path("1AbC", ".cif"),
            Path("ab/1abc.cif"))

        self.assertEqual(
            pdb.pdb_path("1abC", ".cif.gz", "A"),
            Path("ab/1abc/1abc_A.cif.gz"))

    def test_find_pdb(self):
        """Test "find_pdb" function to find a structure."""
        self.assertEqual(
            pdb.find_pdb("12as", base_dir=self.mmcif_dir),
            self.mmcif_dir / "2a/12as.cif")

        self.assertEqual(
            pdb.find_pdb("4n6v", base_dir=self.mmcif_dir),
            self.mmcif_dir / "n6/4n6v.cif.gz")

        self.assertIsNone(pdb.find_pdb("1abc"))

    def test_open_pdb_uncompressed(self):
        """Open an uncompressed mmcif file."""
        pdb_path = pdb.pdb_path("12as", ".cif", base_dir=self.mmcif_dir)
        with pdb.open_pdb(pdb_path) as pdb_in:
            self.assertGreater(len(pdb_in.readlines()), 0)

    def test_open_pdb_compressed(self):
        """Open a compressed mmcif file."""
        pdb_path = pdb.pdb_path("4n6v", ".cif.gz", base_dir=self.mmcif_dir)
        with pdb.open_pdb(pdb_path) as pdb_in:
            self.assertGreater(len(pdb_in.readlines()), 0)

    def test_open_nonexistent_pdb(self):
        """Try and open a non-existent PDB file."""
        with tempfile.TemporaryDirectory() as tempdir:
            with self.assertRaises(FileNotFoundError):
                with pdb.open_pdb(Path(tempdir) / "bad") as _:
                    pass

    def test_open_nonexistent_pdb_gz(self):
        """Try and open a non-existent gzipped PDB file."""
        with tempfile.TemporaryDirectory() as tempdir:
            with self.assertRaises(FileNotFoundError):
                with pdb.open_pdb(Path(tempdir) / "bad.gz") as _:
                    pass

    def test_open_str(self):
        """Should be able to pass a string, not just a Path object."""
        pdb_path = str(pdb.pdb_path("4n6v", ".cif.gz", base_dir=self.mmcif_dir))
        with pdb.open_pdb(pdb_path) as pdb_in:
            self.assertGreater(len(pdb_in.readlines()), 0)

    def test_atom_seq(self):
        """Get atom sequence from a PDB file."""
        with io.StringIO(ATOM_SEQ_TEST_PDB) as pdb_in:
            structure = Bio.PDB.PDBParser().get_structure("test", pdb_in)
            seq, res = pdb.atom_seq(structure[0]["A"])
            self.assertEqual(
                str(seq), "AG",
                "Filter HETATMs and residues without a CA atom")
            self.assertEqual(res[0].get_id(), (' ', 2, ' '))
            self.assertEqual(res[1].get_id(), (' ', 4, ' '))

    def test_remark_read_write(self):
        """Test that we can read and write custom REMARK fields."""

        # Dump random array to REMARK line and see if we can read it back.
        data = random.sample(range(0, 1000), 1000)
        with io.StringIO() as buffer:
            pdb.write_json_remark(buffer, data, 240)
            buffer.seek(0)
            self.assertListEqual(
                pdb.read_json_remark(buffer, 240), data,
                "Read and write data to REMARK field.")
