"""Test components in the :py:mod:`phyre_engine.component.modelling` module."""

import copy
import datetime
import io
import os
import os.path
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import unittest.mock
import phyre_engine.tools.pdb as pdb
from phyre_engine.tools.template import Template, TemplateDatabase
import phyre_engine.test
from phyre_engine.component.modelling import (HomologyModeller, SoedingSelect,
                                              LoopModel, ModRefiner)
import Bio.PDB.PDBParser
import collections

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
        REMARK 161 [[" ", 4, " "], [" ", 5, " "], ["H_AMP", 331, " "]]
        ATOM      2  CA  ILE A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      7  CA  TYR A   2      15.552  39.410  26.282  1.00  8.51
        HETATM 5121  N   AMP A 331      29.722  14.604   9.802  1.00 12.15
    """)
    _ALIGNMENT = [(2, 1), (3, 2)]


    def template(self):
        parser = Bio.PDB.PDBParser(QUIET=True)
        chain_buf = io.StringIO(self._TEMPLATE_PDB)
        chain = parser.get_structure("", chain_buf)[0]["A"]

        return Template(
            "1xyz", "A",
            chain,
            [[" ", 4, " "], [" ", 5, " "], ["H_AMP", 331, " "]],
            "PY",
            [1, 2])

    def setUp(self):
        """Create a temporary pipeline state in `self.state`."""

        # Save dummy template to library.
        self.pdb_lib_dir = pathlib.Path(tempfile.mkdtemp())
        template_file = (self.pdb_lib_dir / "xy" / "1xyz" / "1xyz_A.pdb")
        template_file.parent.mkdir(parents=True)
        template = self.template()
        with template_file.open("w") as template_out:
            template.write(template_out)

        # Create a temporary template database and add dummy template
        _, db_file = tempfile.mkstemp("-foldlib", "phyreengine-")
        self.template_db = pathlib.Path(db_file)
        TemplateDatabase.create(str(self.template_db))
        template_db = TemplateDatabase(
            str(self.template_db),
            str(self.pdb_lib_dir))
        template_db.add_pdb("1xyz", {
            "deposition_date": datetime.date(2000, 1, 1),
            "release_date": datetime.date(2000, 1, 1),
            "last_update_date": datetime.date(2000, 1, 1),
            "method": "X-RAY CRYSTALLOGRAPHY",
        })
        template_db.add_template(template)
        template_db.commit()

        # Create empty directory for storing models and chdir there.
        self.model_dir = pathlib.Path(tempfile.mkdtemp())
        self.orig_dir = os.getcwd()
        os.chdir(str(self.model_dir))

        # Create pipeline state that can be modified without consequence.
        self.state = {
            "query_sequence": copy.deepcopy(self._QUERY_SEQ),
            "PDB": "1xyz",
            "chain": "A",
            "alignment": copy.deepcopy(self._ALIGNMENT),
            "rank": 1
        }

        self.modeller = HomologyModeller(
            str(self.template_db),
            str(self.pdb_lib_dir))
        self.results = self.modeller.run(self.state)

    def tearDown(self):
        """Remove temporary directories."""
        os.chdir(self.orig_dir)
        shutil.rmtree(str(self.pdb_lib_dir))
        shutil.rmtree(str(self.model_dir))
        os.unlink(str(self.template_db))

    def test_model_exists(self):
        """Generated model must exist in pipeline state and on disk."""
        template = self.results
        self.assertIn(
            "model", template,
            "Added 'model' field.")
        self.assertTrue(
            pathlib.Path(template["model"]).exists(),
            "Created a file for model.")

    def test_model_seq(self):
        """Model must contain residues 2 & 3, mapped to 1 & 2 from template."""
        template = self.results
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

class TestSoedingSelect(unittest.TestCase):
    """Test the SoedingSelect component."""

    # Stripped down for our purposes
    Hit = collections.namedtuple("Hit", "i probab")

    def test_trivial_nooverlap(self):
        """
        Choose structures in trivial alignment with no overlap.

        Neither structure overlaps, so both should be chosen.
        """
        pipeline_state = {
            "sequence": "AAAAA",
            "templates": [
                {"prob": 100.0, "alignment": [
                    self.Hit(1, 1.0),
                    self.Hit(2, 1.0)]},
                {"prob": 100.0, "alignment": [
                    self.Hit(4, 1.0),
                    self.Hit(5, 1.0)]},
            ]
        }
        selector = SoedingSelect()
        results = selector.run(copy.deepcopy(pipeline_state))
        self.assertEqual(results["templates"], pipeline_state["templates"])
        self.assertEqual(
            results["template_at_residue"], [
                (1.0, pipeline_state["templates"][0]),
                (1.0, pipeline_state["templates"][0]),
                None,
                (1.0, pipeline_state["templates"][1]),
                (1.0, pipeline_state["templates"][1])])

    def test_trivial_overlap(self):
        """
        Choose structures with high- and low-scoring overlapping templates.

        The first template has a score of 1.0 for each residue, and the second
        has a score of 0.1 The first template should be chosen and the second
        discarded.
        """

        pipeline_state = {
            "sequence": "AAAAA",
            "templates": [
                {"prob": 100.0,
                 "alignment": [self.Hit(i, 1.0) for i in range(1, 6)]},
                {"prob": 100.0,
                 "alignment": [self.Hit(i, 0.1) for i in range(1, 6)]},
            ]
        }
        selector = SoedingSelect()
        results = selector.run(copy.deepcopy(pipeline_state))
        self.assertEqual(results["templates"],
                         [pipeline_state["templates"][0]])
        self.assertEqual(
            results["template_at_residue"],
            [(1.0, pipeline_state["templates"][0])] * 5)

@phyre_engine.test.requireFields(["bin_dir"], ["tools", "sbg"])
@phyre_engine.test.requireFields(["executable", "config"],
                                 ["tools", "sbg", "loop"])
class TestLoopModel(unittest.TestCase):
    """
    Test :py:class:`phyre_engine.component.modelling.LoopModel` component.

    This test will write a small PDB file missing a single residue to a
    temporary file, and attempt to fill the loop.

    In order to run, this test case needs the following keys to be set in
    the ``tools.sbg`` section of the test configuration (with dots denoting
    subsections):

    ``bin_dir``
        Directory containing the SBG binaries.

    ``loop.executable``
        Name of the loop-modeller executable (e.g ``nova`` or
        ``loop.assembler``).

    ``loop.config``
        Configuration file required by the loop modeller
    """

    SAMPLE_PDB = textwrap.dedent("""\
    ATOM      1  N   ALA A   1      11.751  37.846  29.016
    ATOM      2  CA  ALA A   1      12.501  39.048  28.539
    ATOM      3  C   ALA A   1      13.740  38.628  27.754
    ATOM      4  O   ALA A   1      14.207  37.495  27.890
    ATOM      6  N   TYR A   2      14.235  39.531  26.906
    ATOM      7  CA  TYR A   2      15.552  39.410  26.282
    ATOM      8  C   TYR A   2      16.616  38.913  27.263
    ATOM      9  O   TYR A   2      17.187  37.844  27.068
    ATOM     18  N   ILE A   3      16.789  39.630  28.369
    ATOM     19  CA  ILE A   3      17.791  39.281  29.375
    ATOM     20  C   ILE A   3      17.598  37.844  29.863
    ATOM     21  O   ILE A   3      18.538  37.050  29.854
    ATOM     31  N   LYS A   5      15.825  35.211  28.535
    ATOM     32  CA  LYS A   5      16.137  34.313  27.425
    ATOM     33  C   LYS A   5      17.626  34.197  27.104
    ATOM     34  O   LYS A   5      18.101  33.107  26.798
    ATOM     40  N   GLN A   6      18.356  35.312  27.173
    ATOM     41  CA  GLN A   6      19.794  35.327  26.885
    ATOM     42  C   GLN A   6      20.558  34.436  27.845
    ATOM     43  O   GLN A   6      21.576  33.835  27.497
    ATOM     49  N   ARG A   7      20.058  34.382  29.069
    ATOM     50  CA  ARG A   7      20.706  33.656  30.141
    ATOM     51  C   ARG A   7      20.295  32.189  30.141
    ATOM     52  O   ARG A   7      21.084  31.308  30.497
    """)

    PSSM = textwrap.dedent("""\

    Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
         1 A    1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0.00    0.00
         2 Y    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0.00    0.00
         3 I    0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0    0   0   0   0   0   0   0   0   0  50   0   0   0   0  50   0   0   0   0   0  0.86 0.07
         4 A    1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0   0   0   0   0   0   0   0   0   0   0   0   0   0 100   0   0   0   0   0  2.13 0.07
         5 K    0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0    0   0   0   0   0   0   0   0   0   0   0  29   0   0   0  44   0   0   0  27  0.27 0.14
         6 Q    0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0  44   0   0   0   0   0   0   0   0   0   0   0   0   0   0  56   0   0   0  0.54 0.14
         7 R    0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   26   0   0   0   0  21   0   0   0   0   0   0   0   0   0  53   0   0   0   0  0.49 0.20
    """)

    SAMPLE_SEQ = "AYIAKQR"

    def setUp(self):
        """Create temporary directory containing files for pipeline state."""
        self.tmpdir = pathlib.Path(tempfile.mkdtemp("-test", "phyreengine-"))
        pssm = (self.tmpdir / "pssm.ascii")
        model = (self.tmpdir / "model.pdb")

        with pssm.open("w") as pssm_out:
            pssm_out.write(self.PSSM)
            pssm_out.flush()
        with model.open("w") as model_out:
            model_out.write(self.SAMPLE_PDB)
            model_out.flush()
        self.pipeline = {
            "pssm": {"ascii": str(pssm)},
            "model": str(model),
            "query_sequence": self.SAMPLE_SEQ
        }

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.tmpdir))

    def test_model(self):
        """Model a 1-residue gap."""
        # pylint: disable=unsubscriptable-object
        sbg_config = phyre_engine.test.config["tools"]["sbg"]
        loop_modeller = LoopModel(
            bin_dir=sbg_config["bin_dir"],
            config=sbg_config["loop"]["config"],
            executable=sbg_config["loop"]["executable"])
        results = loop_modeller.run(self.pipeline)

        model = pathlib.Path(results["model"])
        self.assertTrue(model.exists(), "model exists")

        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure("model", str(model))
        residues = list(structure.get_residues())
        self.assertEqual(len(residues), 7)


class TestModRefiner(unittest.TestCase):
    """Test ModRefiner. This is VERY SLOW."""

    BACKBONE = textwrap.dedent("""\
    ATOM      2  CA  ALA A   4      -5.395  15.042   8.562  1.00 29.18
    ATOM      6  CA  ARG A   5      -1.926  13.569   8.780  1.00 22.77
    ATOM     10  CA  TRP A   6      -2.947  11.171  11.533  1.00 22.32
    ATOM     14  CA  CYS A   7      -1.910   7.640  10.705  1.00 22.05
    ATOM     18  CA  SER A   8      -0.780   4.478  12.505  1.00 20.18
    ATOM     22  CA  LEU A   9       2.916   3.961  13.161  1.00 17.82
    ATOM     26  CA  GLU A  10       2.865   0.972  10.783  1.00 26.95
    ATOM     30  CA  GLU A  11       1.248   2.744   7.876  1.00 22.83
    ATOM     34  CA  ALA A  12       3.336   5.907   8.317  1.00 20.91
    ATOM     38  CA  VAL A  13       6.676   4.138   8.535  1.00 19.66
    ATOM     42  CA  ALA A  14       4.826   4.984   3.368  1.00 23.16
    ATOM     46  CA  SER A  15       7.548   7.517   4.253  1.00 12.28
    ATOM     50  CA  ILE A  16      10.522   5.281   5.118  1.00  7.49
    ATOM     54  CA  PRO A  17      12.358   3.030   2.632  1.00  6.40
    ATOM     58  CA  ASP A  18      15.619   1.121   2.504  1.00 10.15
    ATOM     62  CA  GLY A  19      18.333   3.193   0.823  1.00 12.73
    ATOM     66  CA  ALA A  20      17.106   6.351   2.603  1.00  7.12
    """)

    def test_normalise_chains(self):
        """Test chain name normalisation."""
        # Replace chain ID of BACKBONE with a space
        unnormalised = []
        for bb_line in self.BACKBONE.split("\n"):
            unnormalised.append(bb_line[:21] + " " + bb_line[22:])

        # Renormalise chain
        with unittest.mock.patch("Bio.PDB.PDBIO") as pdbio:
            initial = io.StringIO(self.BACKBONE)
            refined = io.StringIO("\n".join(unnormalised))
            ModRefiner.normalise_chains(initial, refined)

            set_structure = pdbio.return_value.set_structure
            set_structure.assert_called_once()
            refined_struc = set_structure.call_args[0][0]
            self.assertEqual(
                list(refined_struc[0].get_chains())[0].get_id(),
                "A")

    @phyre_engine.test.requireFields(["bin_dir", "data_dir"],
                                     ["tools", "modrefiner"])
    def test_modrefiner(self):
        """Run ModRefiner over our sample backbone."""
        orig_dir = os.getcwd()
        with tempfile.TemporaryDirectory("-modrefiner",
                                         "phyreengine-test-") as tmpdir:

            try:
                os.chdir(tmpdir)
                init = pathlib.Path(tmpdir, "init.pdb")

                with init.open("w") as init_out:
                    init_out.write(self.BACKBONE)

                config = phyre_engine.test.config["tools"]["modrefiner"]
                modrefiner = ModRefiner(
                    bin_dir=config["bin_dir"],
                    data_dir=config["data_dir"])
                result = modrefiner.run({"model": init})

                # Model was changed
                self.assertNotEqual(result["model"], init)
                # Model includes non-CA atoms
                struc = Bio.PDB.PDBParser().get_structure(
                        "", result["model"])[0]["A"]
                self.assertIn("CB", struc[4])
            finally:
                os.chdir(orig_dir)

if __name__ == "__main__":
    unittest.main()
