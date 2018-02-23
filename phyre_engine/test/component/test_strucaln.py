"""Test components in the "phyre_engine.component.strucaln" module."""
from pathlib import Path
import shutil
import tempfile
import unittest

import phyre_engine.component.strucaln
import phyre_engine.test
from phyre_engine.test.data.minimal_template import MINIMAL_PDB, CANONICAL_SEQ

class AlnTester(unittest.TestCase):
    """Base class for testing alignment components."""

    def setUp(self):
        """
        Create a temporary directory in which we store a query.pdb file.
        Superpositions should be stored in the temporary directory.
        """
        self.tmpdir = Path(tempfile.mkdtemp("-tmalign", "phyreengine-"))
        self.query = Path(self.tmpdir, "query.pdb")
        with self.query.open("w") as pdb_out:
            pdb_out.write(MINIMAL_PDB)

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(str(self.tmpdir))


@phyre_engine.test.requireFields(["bin_dir"], ["tools", "tmalign"])
class TestTMAlign(AlnTester):
    """Tests for :py:class:`phyre_engine.component.strucaln.TMAlign`."""
    # pylint: disable=unsubscriptable-object

    @classmethod
    def setUpClass(cls):
        """Extract config values."""
        config = phyre_engine.test.config["tools"]["tmalign"]
        cls.bin_dir = config["bin_dir"]
        cls.executable = config.get("executable", "TMalign")

    def align_minimal(self, superposition=None):
        """Align a minimal PDB file with itself."""
        pipeline = {"model": str(self.query), "native": str(self.query)}
        cmpt = phyre_engine.component.strucaln.TMAlign(
            superposition=superposition,
            bin_dir=self.bin_dir,
            executable=self.executable)
        results = cmpt.run(pipeline)
        return results

    def test_alignment(self):
        """Align a minimal PDB with itself."""
        results = self.align_minimal()
        self.assertEqual(results["TM"], (1.0, 1.0))
        self.assertEqual(results["structural_alignment"],
                         (CANONICAL_SEQ, CANONICAL_SEQ))
        self.assertEqual(results["seqid"], 1.0)
        self.assertEqual(results["rmsd"], 0.0)

    def test_superpositions(self):
        """Test that superposition files exist with the correct prefix."""
        sup_prefix = str(self.tmpdir / "sup")
        results = self.align_minimal(sup_prefix)
        expected_supers = {
            "trace": sup_prefix,
            "all": sup_prefix + "_all",
            "atm": sup_prefix + "_atm",
            "all_atm": sup_prefix + "_all_atm",
            "all_atm_lig": sup_prefix + "_all_atm_lig",
        }
        self.assertEqual(results["superposition"], expected_supers)

@phyre_engine.test.requireFields(["bin_dir"], ["tools", "tmascore"])
class TestTMScore(AlnTester):
    """Tests for :py:class:`phyre_engine.component.strucaln.TMAlign`."""
    # pylint: disable=unsubscriptable-object

    @classmethod
    def setUpClass(cls):
        """Extract config values."""
        config = phyre_engine.test.config["tools"]["tmscore"]
        cls.bin_dir = config["bin_dir"]
        cls.executable = config.get("executable", "TMscore")

    def align_minimal(self, superposition=None):
        """Align a minimal PDB file with itself."""
        pipeline = {"model": str(self.query), "native": str(self.query)}
        cmpt = phyre_engine.component.strucaln.TMScore(
            superposition=superposition,
            bin_dir=self.bin_dir,
            executable=self.executable)
        results = cmpt.run(pipeline)
        return results

    def test_alignment(self):
        """Align a minimal PDB with itself."""
        results = self.align_minimal()
        self.assertAlmostEqual(results["TM"], 1.0, 4)
        self.assertAlmostEqual(results["maxsub"], 1.0, 4)
        self.assertAlmostEqual(results["GDT_TS"], 1.0, 4)
        self.assertAlmostEqual(results["GDT_HA"], 1.0, 4)
        self.assertEqual(results["structural_alignment"],
                         (CANONICAL_SEQ, CANONICAL_SEQ))

    def test_superpositions(self):
        """Test that superposition files exist with the correct prefix."""
        sup_prefix = str(self.tmpdir / "sup")
        results = self.align_minimal(sup_prefix)
        expected_supers = {
            "trace": sup_prefix,
            "atm": sup_prefix + "_atm",
        }
        self.assertEqual(results["superposition"], expected_supers)
