import json
import os
import tempfile
import unittest
import phyre_engine.component.db.db as db
import phyre_engine.conformation as conformation
from Bio.PDB import PDBParser
from pathlib import Path

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

class TestChainPDBBuilder(unittest.TestCase):
    """Test ChainPDBBuilder class"""

    @unittest.skipUnless("MMCIF" in os.environ, "MMCIF env var not set")
    def test_build(self):
        """Try and extract a chain"""
        with tempfile.TemporaryDirectory() as base_dir:
            out_dir = Path(base_dir, "pdb")
            map_dir = Path(base_dir, "map")

            out_dir.mkdir(exist_ok=True)
            map_dir.mkdir(exist_ok=True)

            builder = db.ChainPDBBuilder(
                os.environ["MMCIF"],
                out_dir, map_dir,
                conformation.ArbitraryMutationSelector(),
                conformation.ArbitraryConformationSelector())
            results = builder.run({"templates": [{"PDB": "12as", "chain": "A"}]})

            self.assertTrue(Path(out_dir, "2a").exists(), "Created 2a subdir")
            self.assertTrue(
                Path(out_dir, "2a/12as_A.pdb").exists(),
                "Created 12as_A"
            )
            self.assertTrue(Path(map_dir, "2a").exists(), "Created 2a map dir")
            self.assertTrue(
                Path(map_dir, "2a/12as_A.json").exists(),
                "Created 12as_A map"
            )
            self.assertEqual(
                results["templates"][0]["sequence"].seq, (
                "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
                "AVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDE"
                "DRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSE"
                "EFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGG"
                "KLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMG"
                "IRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTML"
                "LLQLPHIGQVQAGVWPAAVRESVPSLL"),
                "Got correct sequence")

            # Try and read the map file
            with results["templates"][0]["map"].open("r") as map_fh:
                res_map = json.load(map_fh)
            self.assertEqual(res_map[0][0], " ", "Residue 1 het flag")
            self.assertEqual(res_map[0][1], 4,   "Residue 1 ID")
            self.assertEqual(res_map[0][2], " ", "Residue 1 icode")

    @unittest.skipUnless("MMCIF" in os.environ, "MMCIF env var not set")
    def test_disordered_build(self):
        """Try and extract a chain for a structure with disordered atoms."""
        with tempfile.TemporaryDirectory() as base_dir:
            out_dir = Path(base_dir, "pdb")
            map_dir = Path(base_dir, "map")

            out_dir.mkdir(exist_ok=True)
            map_dir.mkdir(exist_ok=True)

            builder = db.ChainPDBBuilder(
                os.environ["MMCIF"],
                out_dir, map_dir,
                conformation.ArbitraryMutationSelector(),
                conformation.ArbitraryConformationSelector())
            results = builder.run({"templates": [{"PDB": "4n6v", "chain": "0"}]})

            with results["templates"][0]["structure"].open("r") as pdb_fh:
                structure = PDBParser().get_structure("result", pdb_fh)
                chain = next(structure.get_models())[' ']

            self.assertAlmostEqual(
                chain[1]['N'].get_coord()[0], -5.252,
                3, "X coord of residue 1 correct")
            self.assertAlmostEqual(
                chain[1]['N'].get_coord()[1], -0.725,
                3, "Y coord of residue 1 correct")
            self.assertAlmostEqual(
                chain[1]['N'].get_coord()[2], 4.200,
                3, "Z coord of residue 1 correct")

            self.assertAlmostEqual(
                chain[30]['CA'].get_coord()[0], -31.521,
                3, "X coord of residue 30 correct")
            self.assertAlmostEqual(
                chain[30]['CA'].get_coord()[1], -2.596,
                3, "Y coord of residue 30 correct")
            self.assertAlmostEqual(
                chain[30]['CA'].get_coord()[2], 2.172,
                3, "Z coord of residue 30 correct")

def skipUnlessEnv(*env):
    for e in env:
        if not e in os.environ:
            return unittest.skip("Environment variable {} required.".format(e))
    return lambda func: func

class TestSimpleDBPipeline(unittest.TestCase):
    """Try and run a pipeline to build a simple hhblits db."""

    @skipUnlessEnv("HHLIB", "MMCIF", "HHBLITS_DB")
    def test_pipeline(self):
        """Build and run a pipeline."""

        with tempfile.TemporaryDirectory() as pdb_dir, \
                tempfile.TemporaryDirectory() as database_dir, \
                tempfile.TemporaryDirectory() as map_dir:

            pipeline = [
                db.ChainPDBBuilder(
                    os.environ["MMCIF"],
                    pdb_dir, map_dir,
                    conformation.ArbitraryMutationSelector(),
                    conformation.ArbitraryConformationSelector("A")),
                db.MSABuilder(
                    os.environ["HHBLITS_DB"],
                    basedir=database_dir, cpu="23"),
                db.AddSecondaryStructure(),
                db.HMMBuilder(basedir=database_dir),
                db.CS219Builder(basedir=database_dir),
            ]
            data = {
                "templates":[
                    {"PDB": "12AS", "chain": "A"},
                    {"PDB": "3NGG", "chain": "A"},
                ]
            }
            for component in pipeline:
                data = component.run(data)

            db_builder = db.DatabaseBuilder(
                "test_prefix", overwrite=True,
                basedir=database_dir)

            results = db_builder.run(data)
            self.assertIn("database", results, "database key added")
            for file_type in ("a3m", "hhm", "cs219"):
                db_prefix = "{}_{}".format(results["database"], file_type)
                self.assertTrue(
                    Path("{}.ffindex".format(db_prefix)).exists(),
                    "ffindex file exists")
                self.assertTrue(
                    Path("{}.ffdata".format(db_prefix)).exists(),
                    "ffdata file exists")

if __name__ == "__main__":
    unittest.main()
