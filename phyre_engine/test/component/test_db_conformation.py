import io
import unittest
import phyre_engine.conformation as conformation
from Bio.PDB import PDBParser

from Bio.PDB import PDBIO
import sys

class TestMutationSelectors(unittest.TestCase):
    """Test point mutation selectors"""

    # Adapted from 3JQH
    MUTATED_PDB = """
ATOM      2  CA APRO A   1       3.746  20.507  21.289  0.83 65.19           C  
ATOM      9  CA BSER A   1       3.772  20.496  21.302  0.17 64.89           C  
ATOM     15  CA  GLU A   2       6.861  18.244  21.377  1.00 60.65           C  
ATOM     24  CA ALYS A   3       7.680  14.952  23.094  0.50 60.74           C  
ATOM     25  CA BLYS A   3       7.674  14.952  23.095  0.50 60.74           C  
"""
    def setUp(self):
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            self.pdb = PDBParser().get_structure("3jqh", string_fh)

    def test_arbitrary(self):
        """Test ArbitraryMutationSelector."""
        selector = conformation.ArbitraryMutationSelector()
        conf = selector.select(self.pdb[0]["A"])[0]

        self.assertEqual(len(list(conf.get_residues())), 3, "Got 3 residues.")
        # Check all residues are ordered
        for residue in conf:
            self.assertTrue(
                residue.is_disordered() != 2,
                "All residues ordered")
            self.assertEqual(
                len(list(residue.get_atom())),
                1,
                "Residue contains 1 atom")

    def test_all(self):
        """Test AllMutationSelector."""
        selector = conformation.AllMutationSelector()
        all_confs = selector.select(self.pdb[0]["A"])

        self.assertEqual(len(all_confs), 2, "Got two conformations.")

        for conf in all_confs:
            self.assertEqual(
                len(list(conf.get_residues())),
                3,
                "Got 3 residues for conformation 1.")

            # Check all residues are ordered
            for residue in conf:
                self.assertTrue(
                    residue.is_disordered() != 2,
                    "All residues ordered")
                self.assertEqual(
                    len(list(residue.get_atom())),
                    1,
                    "Residue contains 1 atom")

class TestMicroConformationSelectors(unittest.TestCase):
    """Test microheterogeneity selectors."""

    # Adapted from 3JQH
    MUTATED_PDB = """
ATOM     23  N   LYS A   3       7.510  16.054  22.150  1.00 60.70           N  
ATOM     24  CA ALYS A   3       7.680  14.952  23.094  0.50 60.74           C  
ATOM     25  CA BLYS A   3       7.674  14.952  23.095  0.50 60.74           C  
ATOM     26  C   LYS A   3       8.383  15.368  24.391  1.00 66.77           C  
ATOM     27  O   LYS A   3       7.975  14.975  25.493  1.00 59.94           O  
ATOM     28  CB ALYS A   3       8.460  13.807  22.435  0.50 58.21           C  
ATOM     29  CB BLYS A   3       8.440  13.800  22.434  0.50 58.22           C  
ATOM     30  CG ALYS A   3       7.594  12.729  21.806  0.50 60.30           C  
ATOM     31  CG BLYS A   3       7.602  12.930  21.510  0.50 60.56           C  
ATOM     32  CD ALYS A   3       8.403  11.460  21.582  0.50 56.61           C  
ATOM     33  CD BLYS A   3       8.400  11.726  21.036  0.50 57.10           C  
ATOM     34  CE ALYS A   3       7.514  10.241  21.575  0.50 50.44           C  
ATOM     35  CE BLYS A   3       7.503  10.670  20.433  0.50 54.84           C  
ATOM     36  NZ ALYS A   3       8.308   8.989  21.730  0.50 63.46           N  
ATOM     37  NZ BLYS A   3       8.293   9.541  19.864  0.50 65.48           N
"""

    def _verify_num_elements(self, conf):
        self.assertEqual(len(list(conf.get_atoms())), 9, "Got 9 atoms.")
        for residue in conf:
            self.assertEqual(
                residue.is_disordered(),
                0,
                "Residue is not disordered.")

            for atom in residue:
                self.assertEqual(
                    atom.is_disordered(),
                    0,
                    "Atom is not disordered.")

    def test_arbitrary(self):
        """Test ArbitraryConformationSelector."""
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            selector = conformation.ArbitraryConformationSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)

    def test_population(self):
        """Test PopulationConformationSelector."""

        # We should get conformation A, because it has a lower cumulative B
        # factor.
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            selector = conformation.PopulationConformationSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)
            # Verify that this is conformation A
            self.assertAlmostEqual(
                conf[3]["CB"].get_coord()[0],
                8.460, 3, "Chose conformation A")

        # Remove a CA atom from conformation A. We should pick conformation B.
        no_ca_a = "\n".join([
            ln
            for ln in self.MUTATED_PDB.split("\n")
            if " CA ALYS" not in ln])

        with io.StringIO(no_ca_a) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            selector = conformation.PopulationConformationSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)
            # Verify that this is conformation B
            self.assertAlmostEqual(
                conf[3]["CB"].get_coord()[0],
                8.440, 3, "Chose conformation B")

        # Alter the occupancy of conformation A slightly. We should pick
        # conformation B.
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            pdb[0]["A"][3]["CA"].disordered_get("A").set_occupancy(0.49)
            selector = conformation.PopulationConformationSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)
            # Verify that this is conformation B
            self.assertAlmostEqual(
                conf[3]["CB"].get_coord()[0],
                8.440, 3, "Chose conformation B")

        # Alter the bfactors of conformation A slightly. We should pick
        # conformation B.
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            pdb[0]["A"][3]["CA"].disordered_get("A").set_bfactor(70.74)
            selector = conformation.PopulationConformationSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)
            # Verify that this is conformation B
            self.assertAlmostEqual(
                conf[3]["CB"].get_coord()[0],
                8.440, 3, "Chose conformation B")

if __name__ == "__main__":
    unittest.main()
