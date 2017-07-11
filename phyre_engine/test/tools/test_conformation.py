import io
import unittest
import phyre_engine.tools.conformation as conformation
from Bio.PDB import PDBParser
from phyre_engine.tools.conformation import (PopulationConformationSelector,
                                             PopulationMutationSelector)
import Bio.PDB.Atom

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

    def test_population(self):
        sel = PopulationMutationSelector()
        sanitised = sel.select(self.pdb[0]["A"])

        # We should have chosen the proline, because it has the highest
        # occupancy
        self.assertEqual(sanitised[1].get_resname(), "PRO", "Highest occupancy")

        # This should not have affected the disordered status of the lysine
        self.assertTrue(sanitised[3].is_disordered(), "LYS still disordered")

        # If we add a CB atom to the SER, it should increase the population and
        # we should select it.
        add_cb = Bio.PDB.Atom.DisorderedAtom("CB")
        add_cb.disordered_add(
            Bio.PDB.Atom.Atom("CB", [0] * 3, 65, 0.5, "B", "CB", 3, "C"))
        self.pdb[0]["A"][1].disordered_get("SER").add(add_cb)
        sanitised = sel.select(self.pdb[0]["A"])
        self.assertEqual(sanitised[1].get_resname(), "SER", "Highest weight")

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

    def _compare_scores(self, scores1, scores2, message):
        """
        Check various operators for ConformationScores.

        The first parameter *must* be less than the second. Each parameter
        should be a tuple of values that will be passed to the
        :py:class:`phyre_engine.tools.conformation.PopulationConfirmationSelector.ConformationScore`
        constructor.
        """

        # Alias for shorter class name
        Score = PopulationConformationSelector.ConformationScore
        conf_scores1 = Score(*scores1)
        conf_scores2 = Score(*scores2)

        # For constructing diagnostic messages
        def _msg(operator):
            return "{}: {} operator".format(message, operator)

        self.assertLess(conf_scores1, conf_scores2, _msg("<"))
        self.assertLessEqual(conf_scores1, conf_scores2, _msg("<="))
        self.assertGreater(conf_scores2, conf_scores1, _msg(">"))
        self.assertGreaterEqual(conf_scores2, conf_scores1, _msg(">="))
        self.assertNotEqual(conf_scores1, conf_scores2, _msg("!="))
        self.assertEqual(conf_scores1, Score(*scores1), _msg("=="))
        self.assertEqual(conf_scores2, Score(*scores2), _msg("=="))

        # Try and sort each score, and check that the first score is less than
        # the second
        sorted_scores = sorted([conf_scores1, conf_scores2])
        self.assertLess(
            sorted_scores[0],
            sorted_scores[1],
            "{}: sorting".format(message))

    def test_conformation_sort(self):
        """Test comparison and sorting of ConformationScore objects."""
        self._compare_scores((1, 2, 3), (2, 2, 3), "Different population")
        self._compare_scores((1, 1, 3), (1, 2, 3), "Different occupancy")
        self._compare_scores((1, 2, 4), (1, 2, 3), "Different b-factor")

    def test_conformation_add(self):
        """Ensure that addition of ComparisonScore objects works properly."""
        Score = PopulationConformationSelector.ConformationScore
        self.assertEqual(
            Score(1, 2, 3) + Score(1, 2, 3),
            Score(2, 4, 6),
            "Add two scores together")

    def test_population(self):
        """Test PopulationMicroHetSelector."""

        # We should get conformation A, because it has a lower cumulative B
        # factor.
        with io.StringIO(self.MUTATED_PDB) as string_fh:
            pdb = PDBParser().get_structure("3jqh", string_fh)
            selector = conformation.PopulationMicroHetSelector()
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
            selector = conformation.PopulationMicroHetSelector()
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
            selector = conformation.PopulationMicroHetSelector()
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
            selector = conformation.PopulationMicroHetSelector()
            conf = selector.select(pdb[0]["A"])
            self._verify_num_elements(conf)
            # Verify that this is conformation B
            self.assertAlmostEqual(
                conf[3]["CB"].get_coord()[0],
                8.440, 3, "Chose conformation B")

if __name__ == "__main__":
    unittest.main()
