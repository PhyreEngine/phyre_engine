import io
import os
import unittest
from pathlib import Path
import tempfile
import Bio.PDB

import phyre_engine.test
import phyre_engine.tools.rotamer as rot
from phyre_engine.tools.rotamer.angle_range import AngleRange
from phyre_engine.tools.rotamer.data.molprobity import ROTAMERS, SYMMETRIC_FINAL_CHI

DATA_DIR = Path(phyre_engine.test.__file__).parent / "data"
DATA_DIR_2 = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')

class TestSidechain(unittest.TestCase):
    """Test parsing of a single side-chain."""

    VALINE = """
ATOM     94  N   VAL A   1      22.212  26.467  28.766  1.00  9.92           N
ATOM     95  CA  VAL A   1      22.047  25.268  27.956  1.00  6.81           C
ATOM     96  C   VAL A   1      23.313  24.909  27.197  1.00 10.24           C
ATOM     97  O   VAL A   1      23.823  23.796  27.329  1.00 23.59           O
ATOM     98  CB  VAL A   1      20.888  25.446  26.949  1.00 10.72           C
ATOM     99  CG1 VAL A   1      20.769  24.220  26.049  1.00  8.13           C
ATOM    100  CG2 VAL A   1      19.597  25.667  27.701  1.00  2.66           C
"""

    LYSINE = """
ATOM     31  N   LYS A   1      15.825  35.211  28.535  1.00  9.25           N
ATOM     32  CA  LYS A   1      16.137  34.313  27.425  1.00  4.96           C
ATOM     33  C   LYS A   1      17.626  34.197  27.104  1.00 11.71           C
ATOM     34  O   LYS A   1      18.101  33.107  26.798  1.00  2.00           O
ATOM     35  CB  LYS A   1      15.384  34.759  26.166  1.00  9.89           C
ATOM     36  CG  LYS A   1      14.966  33.618  25.258  1.00 12.74           C
ATOM     37  CD  LYS A   1      14.121  34.127  24.100  1.00  9.17           C
ATOM     38  CE  LYS A   1      14.967  34.939  23.141  1.00  6.88           C
ATOM     39  NZ  LYS A   1      14.138  35.863  22.307  1.00 17.64           N
"""
    def test_glycine(self):
        res = Bio.PDB.Residue.Residue(1, "GLY", 0)
        self.assertIsNone(
            rot.Sidechain.calculate(res, SYMMETRIC_FINAL_CHI),
            "GLY has no SC")
        self.assertTupleEqual(
            rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI),
            tuple(),
            "GLY has no chi angles")

    def test_alanine(self):
        res = Bio.PDB.Residue.Residue(1, "ALA", 0)
        self.assertIsNone(
            rot.Sidechain.calculate(res, SYMMETRIC_FINAL_CHI),
            "ALA has no SC")
        self.assertTupleEqual(
            rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI),
            tuple(),
            "ALA has no chi angles")

    def test_unknown(self):
        res = Bio.PDB.Residue.Residue(1, "HOH", 0)
        with self.assertRaises(rot.UnknownResidueType) as err:
            rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI)
            self.assertEqual(
                err.exception.type, "HOH", "Exception residue type")

        with self.assertRaises(rot.UnknownResidueType) as err:
            rot.Sidechain.calculate(res, SYMMETRIC_FINAL_CHI)

    def test_valine(self):
        """Parse an atom with a single χ angle."""
        with io.StringIO(self.VALINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("val", string_fh)
            res = list(pdb.get_residues())[0]
            angles = rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI)

            # Instantiate objects and compare them
            sc1 = rot.Sidechain("VAL", angles)
            sc2 = rot.Sidechain.calculate(res, SYMMETRIC_FINAL_CHI)
            self.assertTrue(sc1.isclose(sc2), "Instantiated side-chains match")

            # Test angle parsing
            self.assertEqual(len(angles), 1, "Valine has 1 χ angle")
            self.assertAlmostEqual(angles[0], 176.8, 1, "Parsed χ1 of Valine")

    def test_lysine(self):
        """Parse an atom with four χ angles."""
        with io.StringIO(self.LYSINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("lys", string_fh)
            res = list(pdb.get_residues())[0]
            angles = rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI)

            # Instantiate objects and compare them
            sc1 = rot.Sidechain("LYS", angles)
            sc2 = rot.Sidechain.calculate(res, SYMMETRIC_FINAL_CHI)

            # Check angles are good against known values
            self.assertTrue(sc1.isclose(sc2), "Instantiated side-chains match")
            self.assertEqual(len(angles), 4, "Lysine has 4 χ angles")
            self.assertAlmostEqual(angles[0], 213.1, 1, "Parsed χ1 of Lysine")
            self.assertAlmostEqual(angles[1], 175.0, 1, "Parsed χ2 of Lysine")
            self.assertAlmostEqual(angles[2], 69.7, 1, "Parsed χ3 of Lysine")
            self.assertAlmostEqual(angles[3], 201.5, 1, "Parsed χ4 of Lysine")

    def test_missing_valine(self):
        """Check an error is thrown when parsing small
        residues (VAL) with missing atoms."""
        with io.StringIO(self.VALINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("val", string_fh)
            res = list(pdb.get_residues())[0]
            del res["CB"]
        with self.assertRaises(rot.MissingAtomError) as err:
            rot.Sidechain.calculate_angles(res, SYMMETRIC_FINAL_CHI)
            self.check_exception(err.exception, res, "CB", 0)

    def test_missing_lysine(self):
        """Check an error is thrown when parsing large residues (LYS) with
        missing atoms."""
        with io.StringIO(self.LYSINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("yls", string_fh)
            res = list(pdb.get_residues())[0]

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = res.copy()
            del copy["CB"]
            rot.Sidechain.calculate_angles(copy, SYMMETRIC_FINAL_CHI)
            self.check_exception(err.exception, copy, "CB", 0)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = res.copy()
            del copy["CD"]
            rot.Sidechain.calculate_angles(copy, SYMMETRIC_FINAL_CHI)
            self.check_exception(err.exception, copy, "CD", 1)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = res.copy()
            del copy["CE"]
            rot.Sidechain.calculate_angles(copy, SYMMETRIC_FINAL_CHI)
            self.check_exception(err.exception, copy, "CE", 2)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = res.copy()
            del copy["NZ"]
            rot.Sidechain.calculate_angles(copy, SYMMETRIC_FINAL_CHI)
            self.check_exception(err.exception, copy, "NZ", 3)


    def check_exception(self, err, res, atom, chi):
        self.assertIs(err.residue, res, "Error references residue")
        self.assertEqual(err.atom, atom, "Error specifies missing atom")
        self.assertEqual(err.chi, chi, "Error specifies correct chi angle")

class TestRotamer(unittest.TestCase):
    """Test rotameric state lookup."""
    test_angles = {
        'ASN': {'p-10': [113.8016194939627, 184.31024615711995]},
        'ASP': {'m-20': [273.8578642811458, 99.9528171316109],
                'p30': [12.358138517160553, 58.38652535571098]},
        'GLU': {'mm-40': [301.9582353204332, 308.24225036375265, 49.574469863477276]},
        'HIS': {'m170': [292.52864130315095, 164.28087964338357]},
        'ILE': {'tt': [186.7757499110814, 203.0559073373177]},
        'LEU': {'mt': [304.7679950634703, 209.51610172332454]},
        'MET': {'mmt': [309.79765829339635, 251.07803634775843, 146.83305424774474]},
        'SER': {'t': [138.10591036525867]},
        'TRP': {'m0': [335.30908697832575, 34.39320381502648],
                't-105': [239.46934245004968, 213.60986822403243]},
        'TYR': {'m-85': [259.57791656933944, 85.09418235484242]}
    }

    def test_valid_angles(self):
        """Test that each of our test angles is valid for that rotamer."""
        for aa, rots in self.test_angles.items():
            for rotamer, angles in rots.items():
                ranges = ROTAMERS[aa][rotamer]
                self.assertTrue(
                    rot.Rotamer(aa, rotamer, ranges).valid_chi(angles),
                    "Rotamer {} of AA {} is valid".format(rotamer, aa))
        self.assertFalse(
            rot.Rotamer("MET", "mmt", ROTAMERS["MET"]["mmt"]).valid_chi(
                (10, 10, 10)), "Invalid rotamer")

    def test_find_rotamers(self):
        """Test that we find the correct rotamer for each angle."""
        for aa, rots in self.test_angles.items():
            for rotamer, angles in rots.items():
                self.assertEqual(
                    rotamer,
                    rot.Rotamer.find(rot.Sidechain(aa, angles), ROTAMERS).name,
                    "Found rotamer {} of {}".format(rotamer, aa))

class TestAngleRange(unittest.TestCase):
    """Test AngleRange class."""

    def test_angrange(self):
        ang_range = AngleRange((0, 15), (345, 360))
        self.assertIn(0, ang_range, "0 in 0-15")
        self.assertNotIn(15, ang_range, "15 not in 0-15")
        self.assertNotIn(45.5, ang_range, "45.5 not in 0-15 or 345-360")
        self.assertIn(350, ang_range, "350 in 340-360")
        self.assertIn(360, ang_range, "360 in 0-15")
        self.assertIn(361, ang_range, "361 in 0-15")
