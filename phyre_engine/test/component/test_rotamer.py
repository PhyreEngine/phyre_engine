"""Unit tests for phyre_engine.component.rotamer package"""
import tempfile
import textwrap
import unittest

from Bio.PDB.Residue import Residue

from phyre_engine.component.rotamer.extract import AngleExtractor,\
    AssignRotamers
from phyre_engine.component.rotamer.parse import CsvParser
from phyre_engine.tools.rotamer.rotamer import Sidechain, Rotamer
from phyre_engine.tools.rotamer.data import molprobity, dunbrack
from phyre_engine.tools.rotamer.data.molprobity import FINAL_CHI_RANGE
from phyre_engine.component.rotamer.db import GroupRotamers


class TestAngleExtractor(unittest.TestCase):
    """Test AngleExtractor component."""

    PDB_STRING = textwrap.dedent("""
        ATOM      1  N   MET A   1      21.005  50.503  57.298  1.00 22.87           N
        ATOM      2  CA  MET A   1      19.791  50.519  58.057  1.00 21.80           C
        ATOM      3  C   MET A   1      18.895  49.476  57.441  1.00 18.90           C
        ATOM      4  O   MET A   1      19.052  49.212  56.248  1.00 18.19           O
        ATOM      5  CB  MET A   1      19.089  51.837  57.937  1.00 24.73           C
        ATOM      6  CG  MET A   1      19.802  52.903  58.721  1.00 30.68           C
        ATOM      7  SD  MET A   1      18.587  54.215  58.893  1.00 37.12           S
        ATOM      8  CE  MET A   1      18.466  54.571  57.161  1.00 36.47           C
        ATOM      0  H1  MET A   1      21.727  50.063  57.832  1.00 22.87           H   new
        ATOM      0  H2  MET A   1      20.861  49.994  56.449  1.00 22.87           H   new
        ATOM      0  H3  MET A   1      21.276  51.441  57.082  1.00 22.87           H   new
        ATOM      0  HA  MET A   1      20.015  50.335  59.118  1.00 21.80           H   new
        ATOM      0  HB2 MET A   1      19.033  52.131  56.879  1.00 24.73           H   new
        ATOM      0  HB3 MET A   1      18.055  51.740  58.300  1.00 24.73           H   new
        ATOM      0  HG2 MET A   1      20.130  52.528  59.702  1.00 30.68           H   new
        ATOM      0  HG3 MET A   1      20.701  53.255  58.194  1.00 30.68           H   new
        ATOM      0  HE1 MET A   1      17.743  55.385  57.003  1.00 36.47           H   new
        ATOM      0  HE2 MET A   1      19.451  54.877  56.779  1.00 36.47           H   new
        ATOM      0  HE3 MET A   1      18.128  53.672  56.624  1.00 36.47           H   new
        ATOM      9  N   PHE A   2      18.035  48.845  58.235  1.00 15.21           N
        ATOM     10  CA  PHE A   2      16.979  48.002  57.683  1.00 13.93           C
        ATOM     11  C   PHE A   2      15.953  48.976  57.090  1.00 13.20           C
        ATOM     12  O   PHE A   2      15.835  50.110  57.589  1.00 13.99           O
        ATOM     13  CB  PHE A   2      16.279  47.215  58.769  1.00 12.79           C
        ATOM     14  CG  PHE A   2      17.064  46.021  59.267  1.00 13.36           C
        ATOM     15  CD1 PHE A   2      18.148  46.194  60.112  1.00 13.47           C
        ATOM     16  CD2 PHE A   2      16.675  44.753  58.867  1.00 13.54           C
        ATOM     17  CE1 PHE A   2      18.834  45.073  60.549  1.00 15.48           C
        ATOM     18  CE2 PHE A   2      17.363  43.637  59.304  1.00 12.97           C
        ATOM     19  CZ  PHE A   2      18.442  43.800  60.147  1.00 14.30           C
        ATOM      0  H   PHE A   2      18.047  48.899  59.233  1.00 15.21           H   new
        ATOM      0  HA  PHE A   2      17.398  47.292  56.955  1.00 13.93           H   new
        ATOM      0  HB2 PHE A   2      16.072  47.884  59.617  1.00 12.79           H   new
        ATOM      0  HB3 PHE A   2      15.307  46.869  58.388  1.00 12.79           H   new
        ATOM      0  HD1 PHE A   2      18.457  47.201  60.429  1.00 13.47           H   new
        ATOM      0  HD2 PHE A   2      15.811  44.634  58.196  1.00 13.54           H   new
        ATOM      0  HE1 PHE A   2      19.697  45.190  61.221  1.00 15.48           H   new
        ATOM      0  HE2 PHE A   2      17.054  42.631  58.983  1.00 12.97           H   new
        ATOM      0  HZ  PHE A   2      18.996  42.919  60.504  1.00 14.30           H   new
        ATOM     20  N   LYS A   3      15.214  48.597  56.058  1.00 12.06           N
        ATOM     21  CA  LYS A   3      14.164  49.432  55.491  1.00 13.07           C
        ATOM     22  C   LYS A   3      12.869  48.669  55.651  1.00 12.04           C
        ATOM     23  O   LYS A   3      12.774  47.509  55.249  1.00 12.30           O
        ATOM     24  CB  LYS A   3      14.335  49.676  54.017  1.00 16.20           C
        ATOM     25  CG  LYS A   3      15.413  50.674  53.741  1.00 21.66           C
        ATOM     26  CD  LYS A   3      15.209  51.107  52.296  1.00 27.06           C
        ATOM     27  CE  LYS A   3      16.234  52.178  51.986  1.00 30.71           C
        ATOM     28  NZ  LYS A   3      15.769  52.985  50.874  1.00 35.99           N
        ATOM      0  H   LYS A   3      15.323  47.716  55.597  1.00 12.06           H   new
        ATOM      0  HA  LYS A   3      14.187  50.404  56.006  1.00 13.07           H   new
        ATOM      0  HB2 LYS A   3      14.575  48.728  53.514  1.00 16.20           H   new
        ATOM      0  HB3 LYS A   3      13.386  50.034  53.592  1.00 16.20           H   new
        ATOM      0  HG2 LYS A   3      15.342  51.532  54.425  1.00 21.66           H   new
        ATOM      0  HG3 LYS A   3      16.409  50.230  53.884  1.00 21.66           H   new
        ATOM      0  HD2 LYS A   3      15.328  50.251  51.616  1.00 27.06           H   new
        ATOM      0  HD3 LYS A   3      14.190  51.494  52.150  1.00 27.06           H   new
        ATOM      0  HE2 LYS A   3      16.397  52.812  52.870  1.00 30.71           H   new
        ATOM      0  HE3 LYS A   3      17.201  51.715  51.739  1.00 30.71           H   new
        ATOM      0  HZ1 LYS A   3      16.446  53.692  50.669  1.00 35.99           H   new
        ATOM      0  HZ2 LYS A   3      15.641  52.404  50.070  1.00 35.99           H   new
        ATOM      0  HZ3 LYS A   3      14.899  53.416  51.115  1.00 35.99           H   new
        """)

    def test_angle_extractor(self):
        with tempfile.NamedTemporaryFile("w") as tmp_pdb:
            print(self.PDB_STRING, file=tmp_pdb)
            tmp_pdb.flush()

            angle_extractor = AngleExtractor(FINAL_CHI_RANGE)
            data = {"pdbs": [tmp_pdb.name]}
            results = angle_extractor.run(data)

            self.assertEqual(
                results["residues"][0]["residue"].get_resname(),
                "PHE",
                "Examined correct residue")
            self.assertEqual(
                results["residues"][0]["source"],
                tmp_pdb.name,
                "Source file is correct")
            self.assertEqual(
                results["residues"][0]["chain"],
                "A",
                "Chain ID is correct")

            self.assertAlmostEqual(
                results["residues"][0]["torsion"][0],
                - 73.44,
                1, "Phi angle")
            self.assertAlmostEqual(
                results["residues"][0]["torsion"][1],
                150.9,
                1, "Psi angle")
            self.assertAlmostEqual(
                results["residues"][0]["sidechain"].angles[0],
                283.0,
                1, "chi-1 angle")
            self.assertAlmostEqual(
                results["residues"][0]["sidechain"].angles[1],
                75.6,
                1, "chi-2 angle")

class TestRotamerAssignment(unittest.TestCase):
    """Check that we assign rotamers correctly"""

    def test_molprobity(self):
        """Check assignment using MolProbity definitions."""
        residues = [{
            "sidechain": Sidechain(
                "ASN",
                (113.8016194939627, 184.31024615711995))
        }]
        component = AssignRotamers(molprobity.ROTAMERS)
        results = component.run({"residues": residues})
        self.assertEqual(
            results["residues"][0]["rotamer"].name,
            "p-10",
            "Correctly assigned p-10 rotamer")

    def test_dunbrack(self):
        """Check assignment using MolProbity definitions."""
        residues = [{
            "sidechain": Sidechain(
                "ASN",
                (113.8016194939627, 84.31024615711995))
        }]
        component = AssignRotamers(dunbrack.ROTAMERS)
        results = component.run({"residues": residues})
        self.assertTupleEqual(
            results["residues"][0]["rotamer"].name,
            (1, 1),
            "Correctly assigned (1,1) rotamer")

class TestCsvParser(unittest.TestCase):
    """Test that we can parse CSV files correctly."""

    CSV_STRING = textwrap.dedent("""\
        NONSENSE_FIELD Phi Psi Chi0 Chi1 Chi2 Chi3 PDB_FILE CHAIN RES_INDEX AA ICODE
        foo -20 30 10 0 0 0 a.pdb A 10 SER X
        bar -15.0 0 180.0 0 0 0 b.pdb B 12 VAL
        baz -180 50 300 -200 102 300 c.pdb A 33A ARG
        """)
    CORRECT_TORSION = [(-20, 30), (-15.0, 0), (-180, 50)]
    CORRECT_CHI = [(10,), (180.0,), (300, -200, 102, 300)]
    CORRECT_RESIDUES = [
        Residue((' ', 10, 'X'), 'SER', 1),
        Residue((' ', 12, ' '), 'VAL', 1),
        Residue((' ', 33, 'A'), 'ARG', 1),
    ]
    CORRECT_CHAINS = ['A', 'B', 'A']
    CORRECT_SOURCES = ['a.pdb', 'b.pdb', 'c.pdb']

    def assertTupleAlmostequal(self, a, b, places, msg):
        for val1, val2 in zip(a, b):
            self.assertAlmostEqual(val1, val2, places, msg)
    def assertResiduesEqual(self, a, b, msg=""):
        self.assertEqual(
            a.get_full_id(),
            b.get_full_id(),
            "{}: Full ID".format(msg))
        self.assertEqual(
            a.get_resname(),
            b.get_resname(),
            "{} Residue name".format(msg))

    def test_parser(self):
        """Parse a dummy CSV file."""
        mapping = {
            "phi": "Phi",
            "psi": "Psi",
            "chi": ("Chi0", "Chi1", "Chi2", "Chi3"),
            "source": "PDB_FILE",
            "chain": "CHAIN",
            "residue": {"id": (None, "RES_INDEX", "ICODE"), "resname": "AA"}
        }
        parser = CsvParser(mapping, {}, delimiter=' ')
        with tempfile.NamedTemporaryFile("w") as csv_file:
            csv_file.write(self.CSV_STRING)
            csv_file.flush()

            results = parser.run({"residue_csv": csv_file.name})
            for i, residue in enumerate(results["residues"]):
                self.assertTupleAlmostequal(
                    residue["torsion"],
                    self.CORRECT_TORSION[i],
                    1,
                    "Residue {} torsion".format(i))
                self.assertTupleAlmostequal(
                    residue["sidechain"].angles,
                    self.CORRECT_CHI[i],
                    1,
                    "Residue {} chi".format(i))
                self.assertEqual(
                    residue["source"],
                    self.CORRECT_SOURCES[i],
                    "Residue {} source".format(i))
                self.assertResiduesEqual(
                    residue["residue"],
                    self.CORRECT_RESIDUES[i],
                    "Residue {}".format(i))
                self.assertEqual(
                    residue["chain"],
                    self.CORRECT_CHAINS[i],
                    "Residue {} chain".format(i))

class TestGroupRotamers(unittest.TestCase):
    """Test the GroupRotamers component"""
    maxDiff = None
    ARG_RES = Residue((" ", 1, " "), "ARG", 1)
    SER_RES = Residue((" ", 1, " "), "SER", 1)

    TEST_RESIDUES = [
        {"residue": ARG_RES, "rotamer": Rotamer("ARG", (1, 1, 1, 1), (0, 360))},
        {"residue": ARG_RES, "rotamer": Rotamer("ARG", (1, 1, 1, 1), (0, 360))},
        {"residue": ARG_RES, "rotamer": Rotamer("ARG", (1, 1, 1, 2), (0, 360))},
        {"residue": SER_RES, "rotamer": Rotamer("SER", (1,), (0, 360))},
        {"residue": SER_RES, "rotamer": Rotamer("SER", (1,), (0, 360))},
        {"residue": SER_RES, "rotamer": Rotamer("SER", (2,), (0, 360))},
    ]

    EXPECTED_RESULT = {
        "ARG": {
            (1, 1, 1, 1): TEST_RESIDUES[0:2],
            (1, 1, 1, 2): TEST_RESIDUES[2:3],
        },
        "SER": {
            (1,): TEST_RESIDUES[3:5],
            (2,): TEST_RESIDUES[5:],
        },
    }

    def test_grouping(self):
        cpt = GroupRotamers()
        result = cpt.run({"residues": self.TEST_RESIDUES})
        self.assertDictEqual(
            result["rotamer_dict"],
            self.EXPECTED_RESULT,
            "Correctly grouped residues")
