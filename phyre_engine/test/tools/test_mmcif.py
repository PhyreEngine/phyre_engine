import os
import unittest
import tempfile
import Bio.PDB

import phyre_engine.test
from phyre_engine.tools.mmcif import MMCIFToPDBChains

DATA_DIR = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')

class TestMMCIFToPDBChains(unittest.TestCase):
    """Test conversion of MMCIF file to multiple chains."""

    def test_convert(self):
        """Try and convert 12as.cif."""
        with tempfile.TemporaryDirectory() as tmpdir:
            conv = MMCIFToPDBChains(tmpdir)
            cif_file = os.path.join(DATA_DIR, "mmcif", "2a", "12as.cif")
            conv.convert("12AS", cif_file)

            chain_a_file = os.path.join(tmpdir, "2a", "12as_A.pdb")
            chain_b_file = os.path.join(tmpdir, "2a", "12as_B.pdb")

            self.assertTrue(
                    os.path.isdir(os.path.join(tmpdir, "2a")),
                    "Created /2a/ directory")
            self.assertTrue(
                    os.path.isfile(chain_a_file),
                    "Created /2a/12as_A.pdb file")
            self.assertTrue(
                    os.path.isfile(chain_b_file),
                    "Created /2a/12as_B.pdb file")

            # Read PDB files back
            pdb_parser = Bio.PDB.PDBParser()
            chain_a = pdb_parser.get_structure("12AS_A", chain_a_file)
            chain_b = pdb_parser.get_structure("12AS_B", chain_b_file)

            model_a = next(chain_a.get_models())
            model_b = next(chain_b.get_models())

            # Each chain should be labelled "A"
            self.assertIn(" ", model_a, "Chain A in B")
            self.assertIn(" ", model_b, "Chain A in B")

            self.assertAlmostEqual(
                    model_a[" "][4]["N"].get_coord()[0],
                    11.751,
                    places=3, msg="X coord of 4N is 11.751")
            self.assertAlmostEqual(
                    model_b[" "][4]["N"].get_coord()[0],
                    51.203,
                    places=3, msg="X coord of 4N is 11.751")
