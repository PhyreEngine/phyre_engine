import os.path
import unittest
import tempfile
import Bio.Seq
import Bio.Alphabet.IUPAC
import Bio.SeqUtils
import phyre_engine.test
from phyre_engine.component.modelling import HomologyModeller;
from phyre_engine.tools.hhsuite.parser import Hit;

DATA_DIR = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')

class TestHomologyModeller(unittest.TestCase):

    def _test_12as(self, mmcif_dir):
        query_seq = Bio.Seq.Seq("AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL",
                Bio.Alphabet.IUPAC.protein)

        hit = Hit(name="12AS_A")
        hit.aln = [
            {"i": 1, "j":4},
            {"i": 2, "j":5},
            {"i": 3, "j":6},
            {"i": 7, "j":7},
            {"i": 8, "j":8},
        ]

        data = {"hits": [hit], "sequence": query_seq}
        modeller = HomologyModeller(mmcif_dir)
        data = modeller.run(data)

        self.assertTrue("models" in data, "models key added")
        self.assertEqual(len(data["models"]), 1, "Got 1 model for 1 hit")

        model = data["models"][0]
        residues = list(model.get_residues())
        self.assertEqual(len(residues), 5, "Got 5 residues")

        query_three = [Bio.SeqUtils.seq3(aa).upper() for aa in query_seq]

        # Make sure template residues were named correctly
        self.assertEqual(residues[0].get_resname(), query_three[0], "Residue 1")
        self.assertEqual(residues[1].get_resname(), query_three[1], "Residue 2")
        self.assertEqual(residues[2].get_resname(), query_three[2], "Residue 3")
        self.assertEqual(residues[3].get_resname(), query_three[6], "Residue 4")
        self.assertEqual(residues[4].get_resname(), query_three[7], "Residue 5")

        # Check that atoms have the correct ID
        self.assertIs(residues[0], model[1]["A"][1], "ID of residue 1")
        self.assertIs(residues[1], model[1]["A"][2], "ID of residue 2")
        self.assertIs(residues[2], model[1]["A"][3], "ID of residue 3")
        self.assertIs(residues[3], model[1]["A"][7], "ID of residue 7")
        self.assertIs(residues[4], model[1]["A"][8], "ID of residue 8")

        # Check coordinates of a few atoms to verify. Check to 3 decimal places,
        # which is all that is given in the MMCIF file.
        atom_1_n_coords = residues[0]["N"].get_coord()
        self.assertAlmostEqual(atom_1_n_coords[0], 11.751, 3, msg="X of 1N")
        self.assertAlmostEqual(atom_1_n_coords[1], 37.846, 3, msg="Y of 1N")
        self.assertAlmostEqual(atom_1_n_coords[2], 29.016, 3, msg="Z of 1N")

        atom_8_ca_coords = residues[4]["CA"].get_coord()
        self.assertAlmostEqual(atom_8_ca_coords[0], 16.137, 3, msg="X of 8CA")
        self.assertAlmostEqual(atom_8_ca_coords[1], 34.313, 3, msg="Y of 8CA")
        self.assertAlmostEqual(atom_8_ca_coords[2], 27.425, 3, msg="Z of 8CA")

    def test_auth_seq_id(self):
        """Build a model given an alignment.

        The MMCIF template file does contain ``auth_seq_id`` fields, so this
        test tests the conversion from sequence index to ``auth_seq_id``.
        """
        self._test_12as(os.path.join(DATA_DIR, "mmcif"))


    # We currently expect a failure when parsing MMCIF files without proper
    # auth fields, because MMCIFParser does not handle them correctly.
    @unittest.expectedFailure
    def test_label_seq_id(self):
        """Build a model given an alignment.

        The MMCIF template file has no auth IDs, so this tests the ability to
        read label ids only.
        """
        self._test_12as(os.path.join(DATA_DIR, "mmcif_label_only"))


if __name__ == "__main__":
    unittest.main()
