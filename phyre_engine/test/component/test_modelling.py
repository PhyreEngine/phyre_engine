"""Test components in the :py:mod:`phyre_engine.component.modelling` module."""

import copy
import os.path
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import phyre_engine.tools.pdb as pdb
import phyre_engine.test
from phyre_engine.component.modelling import HomologyModeller
import Bio.PDB.PDBParser

DATA_DIR = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')

class TestHomologyModeller(unittest.TestCase):
    """
    Test HomologyModeller component.

    This simple test takes a minimal query sequence, template and alignment and
    generates a model. The query sequence ``PAG`` is aligned to a template with
    sequence ``IY`` as follows:

    .. code-block:: none

        >query/template
        PAG
        -IY

    The template has an additional ``AMP`` ``HETATM`` that should not affect the
    modelling.

    The resulting model should contain the residues ``AG`` numbered 2 and 3,
    inheriting coordinates from residue 1 and 2 of the template structure.
    """

    _QUERY_SEQ = "PAG"
    _TEMPLATE_PDB = textwrap.dedent("""\
        REMARK 150 PY
        REMARK 156 [1, 2]
        REMARK 161 [[' ', 4, ' '], [' ', 5, ' '], ['H_AMP', 331, ' ']]
        ATOM      2  CA  ILE A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      7  CA  TYR A   2      15.552  39.410  26.282  1.00  8.51
        HETATM 5121  N   AMP A 331      29.722  14.604   9.802  1.00 12.15
    """)
    _ALIGNMENT = [(2, 1), (3, 2)]

    def setUp(self):
        """Create a temporary pipeline state in `self.state`."""

        # Save dummy template to library.
        self.pdb_lib_dir = pathlib.Path(tempfile.mkdtemp())
        template_file = (self.pdb_lib_dir / "xy" / "1xyz" / "1xyz_A.pdb")
        template_file.parent.mkdir(parents=True)
        with template_file.open("w") as template_out:
            template_out.write(self._TEMPLATE_PDB)

        # Create empty directory for storing models and chdir there.
        self.model_dir = pathlib.Path(tempfile.mkdtemp())
        self.orig_dir = os.getcwd()
        os.chdir(str(self.model_dir))

        # Create pipeline state that can be modified without consequence.
        self.state = {
            "sequence": copy.deepcopy(self._QUERY_SEQ),
            "templates": [
                {
                    "PDB": "1xyz",
                    "chain": "A",
                    "alignment": copy.deepcopy(self._ALIGNMENT),
                    "rank": 1
                }
            ]
        }

        self.modeller = HomologyModeller(str(self.pdb_lib_dir))
        self.results = self.modeller.run(self.state)

    def tearDown(self):
        """Remove temporary directories."""
        os.chdir(self.orig_dir)
        shutil.rmtree(str(self.pdb_lib_dir))
        shutil.rmtree(str(self.model_dir))

    def test_model_exists(self):
        """Generated model must exist in pipeline state and on disk."""
        template = self.results["templates"][0]
        self.assertIn(
            "model", template,
            "Added 'model' field.")
        self.assertTrue(
            pathlib.Path(template["model"]).exists(),
            "Created a file for model.")

    def test_model_seq(self):
        """Model must contain residues 2 & 3, mapped to 1 & 2 from template."""
        template = self.results["templates"][0]
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        model_structure = pdb_parser.get_structure("model", template["model"])
        model_chain = model_structure[0]["A"]
        model_seq, model_residues = pdb.atom_seq(model_chain)
        self.assertEqual(
            model_seq, "AG",
            "Alignment obeyed: model contains residues 2 & 3.")
        self.assertEqual(
            model_residues[0].get_id(),
            (' ', 2, ' '),
            "First residue is residue 2")
        self.assertEqual(
            model_residues[1].get_id(),
            (' ', 3, ' '),
            "Second residue is residue 3")
        self.assertAlmostEqual(
            model_residues[0]["CA"].get_coord()[0], 12.501, 3,
            msg="X coord of residue 2 inherited from template")
        self.assertAlmostEqual(
            model_residues[1]["CA"].get_coord()[0], 15.552, 3,
            msg="X coord of residue 3 inherited from template")

if __name__ == "__main__":
    unittest.main()
