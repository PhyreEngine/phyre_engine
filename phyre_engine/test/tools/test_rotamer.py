import io
import math
import os
import unittest
from pathlib import Path
import tempfile
import Bio.PDB

import numpy as np

import phyre_engine.test
import phyre_engine.tools.rotamer.rotamer as rot
from phyre_engine.tools.rotamer.angle_range import AngleRange
from phyre_engine.tools.rotamer.data.molprobity import ROTAMERS, FINAL_CHI_RANGE
from phyre_engine.tools.rotamer.kdtree import PeriodicKDTree

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

    def setUp(self):
        """Create residue objects for valine and lysine."""
        with io.StringIO(self.VALINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("val", string_fh)
            self.valine = list(pdb.get_residues())[0]
        with io.StringIO(self.LYSINE) as string_fh:
            pdb = Bio.PDB.PDBParser().get_structure("lys", string_fh)
            self.lysine = list(pdb.get_residues())[0]

    def test_glycine(self):
        res = Bio.PDB.Residue.Residue(1, "GLY", 0)
        self.assertIsNone(
            rot.Sidechain.calculate(res, FINAL_CHI_RANGE),
            "GLY has no SC")
        self.assertTupleEqual(
            rot.Sidechain.calculate_angles(res, FINAL_CHI_RANGE),
            tuple(),
            "GLY has no chi angles")

    def test_alanine(self):
        res = Bio.PDB.Residue.Residue(1, "ALA", 0)
        self.assertIsNone(
            rot.Sidechain.calculate(res, FINAL_CHI_RANGE),
            "ALA has no SC")
        self.assertTupleEqual(
            rot.Sidechain.calculate_angles(res, FINAL_CHI_RANGE),
            tuple(),
            "ALA has no chi angles")

    def test_unknown(self):
        res = Bio.PDB.Residue.Residue(1, "HOH", 0)
        with self.assertRaises(rot.UnknownResidueType) as err:
            rot.Sidechain.calculate_angles(res, FINAL_CHI_RANGE)
            self.assertEqual(
                err.exception.type, "HOH", "Exception residue type")

        with self.assertRaises(rot.UnknownResidueType) as err:
            rot.Sidechain.calculate(res, FINAL_CHI_RANGE)

    def test_valine(self):
        """Parse an atom with a single χ angle."""
        angles = rot.Sidechain.calculate_angles(self.valine, FINAL_CHI_RANGE)

        # Instantiate objects and compare them
        sc1 = rot.Sidechain("VAL", angles)
        sc2 = rot.Sidechain.calculate(self.valine, FINAL_CHI_RANGE)
        self.assertTrue(sc1.isclose(sc2), "Instantiated side-chains match")

        # Test angle parsing
        self.assertEqual(len(angles), 1, "Valine has 1 χ angle")
        self.assertAlmostEqual(angles[0], 176.8, 1, "Parsed χ1 of Valine")

    def test_lysine(self):
        """Parse an atom with four χ angles."""
        angles = rot.Sidechain.calculate_angles(self.lysine, FINAL_CHI_RANGE)

        # Instantiate objects and compare them
        sc1 = rot.Sidechain("LYS", angles)
        sc2 = rot.Sidechain.calculate(self.lysine, FINAL_CHI_RANGE)

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
        del self.valine["CB"]
        with self.assertRaises(rot.MissingAtomError) as err:
            rot.Sidechain.calculate_angles(self.valine, FINAL_CHI_RANGE)
            self.check_exception(err.exception, self.valine, "CB", 0)

    def test_missing_lysine(self):
        """Check an error is thrown when parsing large residues (LYS) with
        missing atoms."""
        with self.assertRaises(rot.MissingAtomError) as err:
            copy = self.lysine.copy()
            del copy["CB"]
            rot.Sidechain.calculate_angles(copy, FINAL_CHI_RANGE)
            self.check_exception(err.exception, copy, "CB", 0)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = self.lysine.copy()
            del copy["CD"]
            rot.Sidechain.calculate_angles(copy, FINAL_CHI_RANGE)
            self.check_exception(err.exception, copy, "CD", 1)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = self.lysine.copy()
            del copy["CE"]
            rot.Sidechain.calculate_angles(copy, FINAL_CHI_RANGE)
            self.check_exception(err.exception, copy, "CE", 2)

        with self.assertRaises(rot.MissingAtomError) as err:
            copy = self.lysine.copy()
            del copy["NZ"]
            rot.Sidechain.calculate_angles(copy, FINAL_CHI_RANGE)
            self.check_exception(err.exception, copy, "NZ", 3)

    def test_rotamer_flip(self):
        """Check that the final χ angles are placed in the correct range."""
        final_chi_angles = {"VAL": AngleRange((0, 90), (270, 360))}
        angles = rot.Sidechain.calculate_angles(self.valine, final_chi_angles)
        # Angle should have been flipped by 180°.
        self.assertAlmostEqual(
            angles[0],
            176.8 + 180,
            1, "Flipped χ1 of Valine")

        final_chi_angles = {"VAL": AngleRange((0, 90), (0, 180))}
        angles = rot.Sidechain.calculate_angles(self.valine, final_chi_angles)
        # Angle should have been flipped by 180°.
        self.assertAlmostEqual(
            angles[0],
            176.8,
            1, "Did not flip χ1 of Valine")


    def check_exception(self, err, res, atom, angle_index):
        self.assertIs(err.residue, res, "Error references residue")
        self.assertEqual(err.atom, atom, "Error specifies missing atom")
        self.assertEqual(
            err.angle_index, angle_index, "Error specifies correct chi angle")
        self.assertEqual(
            err.feature, rot.ResidueFeature.CHI, "Error references chi")

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
        self.assertIn(-1, ang_range, "361 in 0-15")

class TestPeriodiKDTree(unittest.TestCase):
    """Test PeriodicKDTree for sanity."""

    def setUp(self):
        self.points = np.random.rand(1000, 2) * 360
        self.tree = PeriodicKDTree(self.points, 360, 360)

    def distance(self, a, b):
        """Get periodic distance between x and y."""
        dx = min(abs(a[0] - b[0]), 360 - abs(a[0] - b[0]))
        dy = min(abs(a[1] - b[1]), 360 - abs(a[1] - b[1]))
        return math.sqrt(dx ** 2 + dy ** 2)

    def test_query_k_points(self):
        """Get 10 points from the tree and verify them."""
        query_point = np.array([0, 0])
        distances, indices = self.tree.query(query_point, 10)

        self.assertEqual(
            len(distances), 10,
            "Got 10 distances from PeriodicKDTree")
        self.assertEqual(
            len(indices), 10,
            "Got 10 indices from PeriodicKDTree")

        all_distances = []
        for i, point in enumerate(self.points):
            calculated_distance = self.distance(query_point, point)
            all_distances.append((i, calculated_distance))
            if i in indices:
                kd_distance = distances[np.where(indices == i)]
                self.assertAlmostEqual(
                    calculated_distance,
                    kd_distance,
                    places=5,
                    msg="PeriodicKDTree and simple distance calculation match")

        # Go through the manually-counted list of distances and check that they
        # are the same points as found by the tree
        all_distances.sort(key=lambda x: x[1])
        self.assertListEqual(
            [x[0] for x in all_distances[0:10]],
            list(indices),
            "PeriodicKDTree returns same points as brute force.")

    def test_query_by_distance(self):
        """Query points from tree by maximum distance."""
        query_point = np.array([0, 0])
        distances, indices = self.tree.query(
            query_point,
            len(self.points),
            distance_upper_bound=10)

        all_distances = []
        for i, point in enumerate(self.points):
            calculated_distance = self.distance(query_point, point)
            all_distances.append((i, calculated_distance))
            if i in indices:
                # Should be closer than 10 units
                kd_distance = distances[np.where(indices == i)]
                self.assertAlmostEqual(
                    calculated_distance,
                    kd_distance,
                    places=5,
                    msg="PeriodicKDTree and simple distance calculation match")
                self.assertLessEqual(
                    calculated_distance, 10,
                    "Distances meets threshold")
            else:
                self.assertGreater(
                    calculated_distance, 10,
                    "Point is correctly marked as being further than 10 units")

        # Go through the manually-counted list of distances and check that they
        # are the same points as found by the tree. We slice indices from 0:-1
        # because scipy adds a final infinite distance to signify exceeding the
        # upper bound.
        close_distances = sorted(
            [x for x in all_distances if x[1] < 10],
            key=lambda x: x[1])
        self.assertListEqual(
            [x[0] for x in close_distances],
            list(indices[0:-1]),
            "PeriodicKDTree returns same points as brute force.")
