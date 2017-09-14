"""Test components in the :py:mod:`phyre_engine.component.secstruc` module."""
from pathlib import Path
import unittest
import phyre_engine.component.secstruc as secstruc
import phyre_engine.test
import phyre_engine.test.data as test_data

PDB_FILE = Path(test_data.__file__).parent / "pdb_chains/2a/12as_A.pdb"

class TestDSSP(unittest.TestCase):
    """Test DSSP component."""

    _DUMMY_DSSP = """\
==== Secondary Structure Definition by the program DSSP, CMBI version by M.L. Hekkelman/2010-10-21 ==== DATE=2017-09-01        .
REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .
                                                                                                                               .
COMPND                                                                                                                         .
SOURCE                                                                                                                         .
AUTHOR                                                                                                                         .
  328  2  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
 14376.9   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .
  208 63.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .
    9  2.7   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
   49 14.9   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
    1  0.3   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
    2  0.6   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2), SAME NUMBER PER 100 RESIDUES                              .
    1  0.3   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1), SAME NUMBER PER 100 RESIDUES                              .
   22  6.7   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2), SAME NUMBER PER 100 RESIDUES                              .
   23  7.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3), SAME NUMBER PER 100 RESIDUES                              .
   90 27.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
    5  1.5   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           .
  0  0  0  0  0  1  2  1  1  0  1  1  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0    RESIDUES PER ALPHA HELIX         .
  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    PARALLEL BRIDGES PER LADDER      .
  2  2  1  0  3  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    ANTIPARALLEL BRIDGES PER LADDER  .
  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    1    1 A A              0   0  107      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0-157.7   12.4   38.7   28.2
    2    2 A Y     >  +     0   0  152      3,-0.1     4,-1.4     2,-0.1     5,-0.0   0.608 360.0  40.8-104.2 -31.9   15.5   39.5   26.2
    3    3 A I  H  > S+     0   0  117      2,-0.2     4,-2.1     3,-0.1     5,-0.2   0.878 114.6  51.3 -83.9 -35.7   17.8   39.1   29.2
    4    4 A A  H  > S+     0   0   43      2,-0.2     4,-1.1     1,-0.2     3,-0.5   0.985 112.5  47.3 -59.5 -53.7   15.9   36.0   30.5
    5    5 A K  H  > S+     0   0   20      1,-0.2     4,-2.4     2,-0.2     3,-0.3   0.851 112.7  50.6 -51.6 -42.4   16.2   34.4   27.1
    6    6 A Q  H  X S+     0   0   55     -4,-1.4     4,-1.4     1,-0.2    -1,-0.2   0.830 109.1  47.9 -70.2 -36.2   19.9   35.3   26.9
    7    7 A R  H  X S+     0   0  194     -4,-2.1     4,-0.7    -3,-0.5    -1,-0.2   0.582 113.9  51.3 -83.1  -7.6   20.8   33.9   30.3
    8    8 A Q  H  X S+     0   0   47     -4,-1.1     4,-2.1    -3,-0.3    -2,-0.2   0.871 104.2  53.5 -91.4 -44.6   18.9   30.7   29.3
    9    9 A I  H  X S+     0   0    8     -4,-2.4     4,-3.0     1,-0.2    -2,-0.2   0.923 111.1  48.8 -53.5 -42.9   20.7   30.2   26.0
   10   10 A S  H  X S+     0   0   63     -4,-1.4     4,-2.6     2,-0.2     5,-0.3   0.852 106.3  53.8 -71.4 -40.0   24.0   30.4   27.8
   11   11 A F  H  X S+     0   0   56     -4,-0.7     4,-2.6     2,-0.2     5,-0.2   0.920 114.3  43.3 -59.0 -44.9   23.1   27.9   30.5
    """

    _DUMMY_STATES = [
        {"assigned": "C", "confidence": {"C": 1.0, "H": 0.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "C", "confidence": {"C": 1.0, "H": 0.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
        {"assigned": "H", "confidence": {"C": 0.0, "H": 1.0, "B": 0.0, "E": 0.0, "G": 0.0, "I": 0.0, "T": 0.0, "S": 0.0}},
    ]

    def test_parsing(self):
        """Test DSSP parsing with some dummy data."""
        dssp_mapping = secstruc.DSSP.parse_dssp(self._DUMMY_DSSP.split("\n"))
        self.assertListEqual(dssp_mapping, self._DUMMY_STATES)

    @phyre_engine.test.requireFields("dssp", ["tools"])
    def test_run(self):
        """Try and run mkdssp on a PDB file."""
        dssp = secstruc.DSSP(**phyre_engine.test.config["tools"]["dssp"])
        results = dssp.run({"structure": str(PDB_FILE)})
        self.assertGreater(len(results["secondary_structure"]["dssp"]), 0)

    @phyre_engine.test.requireFields("dssp", ["tools"])
    def test_length_mismatch(self):
        """Mismatch in sequence/ss length raises exception."""
        dssp = secstruc.DSSP(**phyre_engine.test.config["tools"]["dssp"])
        with self.assertRaises(secstruc.LengthMismatchError):
            dssp.run({"structure": str(PDB_FILE), "sequence": "AAA"})
