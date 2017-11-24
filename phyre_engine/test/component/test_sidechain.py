"""
Test sidechain reconstruction tools in
:py:mod:`phyre_engine.component.sidechain`.
"""
import pathlib
import shutil
import tempfile
import unittest

import Bio.PDB

import phyre_engine.component.sidechain as sidechain
import phyre_engine.test

BACKBONE_ONLY = """\
ATOM      1  N   ALA A   4      11.751  37.846  29.016
ATOM      2  CA  ALA A   4      12.501  39.048  28.539
ATOM      3  C   ALA A   4      13.740  38.628  27.754
ATOM      4  O   ALA A   4      14.207  37.495  27.890
ATOM      6  N   TYR A   5      14.235  39.531  26.906
ATOM      7  CA  TYR A   5      15.552  39.410  26.282
ATOM      8  C   TYR A   5      16.616  38.913  27.263
ATOM      9  O   TYR A   5      17.187  37.844  27.068
ATOM     18  N   ILE A   6      16.789  39.630  28.369
ATOM     19  CA  ILE A   6      17.791  39.281  29.375
ATOM     20  C   ILE A   6      17.598  37.844  29.863
ATOM     21  O   ILE A   6      18.538  37.050  29.854
ATOM     26  N   ALA A   7      16.368  37.519  30.261
ATOM     27  CA  ALA A   7      16.004  36.186  30.742
ATOM     28  C   ALA A   7      16.371  35.097  29.741
ATOM     29  O   ALA A   7      17.143  34.202  30.049
"""

# "CB" is not a backbone atom, but it should appear in all of the
# reconstructed residues.
BACKBONE = {"N", "CA", "C", "O", "CB"}
ATOMS = {
    4: BACKBONE,
    5: BACKBONE | set("H CG HH OH CZ CD1 CE1 CD2 CE2".split()),
    6: BACKBONE | set("H CD1 CG1 CG2".split()),
    # Residue 7 has an "H" because it is at the C-terminus
    7: BACKBONE | {"H"},
}

class BaseSidechainTest(unittest.TestCase):
    """Base test class for writing temporary files and checking results."""

    def setUp(self):
        """Create a temporary directory for input/output files."""
        self.tmpdir = pathlib.Path(tempfile.mkdtemp())
        self.bb_only = self.tmpdir / "bb.pdb"
        with self.bb_only.open("w") as pdb_out:
            pdb_out.write(BACKBONE_ONLY)
            pdb_out.flush()

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(str(self.tmpdir))

    def has_sidechains(self, pdb_file):
        """Check that our simple PDB file has side-chains."""
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("bb", pdb_file)
        residues = list(structure.get_residues())

        # Sanity check: make sure we read in all residues
        self.assertEqual(
            [r.get_resname() for r in residues],
            ["ALA", "TYR", "ILE", "ALA"])

        for res in residues:
            with self.subTest("Contains sidechain", residue=res.get_id()):
                self.assertEqual(
                    set([a.get_name() for a in res]),
                    ATOMS[res.get_id()[1]])

@phyre_engine.test.requireFields("scwrl4", ["tools"])
class TestScwrl4(BaseSidechainTest):
    """Rebuild side-chains with SCWRL4."""

    def test_run(self):
        """SCWRL4 runs and builds sidechains."""
        bin_dir = phyre_engine.test.config["tools"]["scwrl4"]["bin_dir"]
        cmpt = sidechain.Scwrl4(bin_dir=bin_dir)
        results = cmpt.run({"model": str(self.bb_only)})

        # A new file was created:
        self.assertNotEqual(results["model"], str(self.bb_only))
        self.has_sidechains(results["model"])
