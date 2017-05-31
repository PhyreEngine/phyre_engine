import unittest
from Bio.PDB.PDBParser import PDBParser
from io import StringIO
from phyre_engine.component.rotamer.accuracy import ConditionalAccuracy, \
    AbsoluteAccuracy, PerResidueRMSD

# Structures with known chi angles. Indexed by chi angles.
STRUCTURES = {}

STRUCTURES[(-63, 130, 101, 163)] = """\
ATOM      1  N   ARG B   7      27.813  60.399  23.912  1.00  8.49           N  
ATOM      2  CA  ARG B   7      26.950  60.884  24.993  1.00  8.58           C  
ATOM      3  C   ARG B   7      25.739  59.980  25.031  1.00  8.21           C  
ATOM      4  O   ARG B   7      25.115  59.712  24.002  1.00 11.30           O  
ATOM      5  CB  ARG B   7      26.493  62.323  24.783  1.00  9.10           C  
ATOM      6  CG  ARG B   7      27.618  63.347  24.769  1.00 10.66           C  
ATOM      7  CD  ARG B   7      27.431  64.191  23.527  1.00 16.50           C  
ATOM      8  NE  ARG B   7      28.345  63.734  22.529  1.00 15.14           N  
ATOM      9  CZ  ARG B   7      28.320  63.965  21.220  1.00 14.14           C  
ATOM     10  NH1 ARG B   7      27.380  64.671  20.583  1.00 14.22           N1+
ATOM     11  NH2 ARG B   7      29.365  63.544  20.553  1.00 14.84           N
"""
STRUCTURES[(63, 130, 101, 163)] = """\
ATOM      1  N   ARG B   7      27.813  60.399  23.912  1.00  8.49           N  
ATOM      2  CA  ARG B   7      26.950  60.884  24.993  1.00  8.58           C  
ATOM      3  C   ARG B   7      25.739  59.980  25.031  1.00  8.21           C  
ATOM      4  O   ARG B   7      25.115  59.712  24.002  1.00 11.30           O  
ATOM      5  CB  ARG B   7      26.493  62.323  24.783  1.00  9.10           C  
ATOM      6  CG  ARG B   7      25.640  62.539  23.542  1.00 10.66           C  
ATOM      7  CD  ARG B   7      24.406  63.296  23.982  1.00 16.50           C  
ATOM      8  NE  ARG B   7      23.326  62.369  24.100  1.00 15.14           N  
ATOM      9  CZ  ARG B   7      22.170  62.525  24.739  1.00 14.14           C  
ATOM     10  NH1 ARG B   7      21.818  63.612  25.433  1.00 14.22           N1+
ATOM     11  NH2 ARG B   7      21.296  61.566  24.565  1.00 14.84           N
"""
STRUCTURES[(63, -110, 101, 163)] = """\
ATOM      1  N   ARG B   7      27.813  60.399  23.912  1.00  8.49           N  
ATOM      2  CA  ARG B   7      26.950  60.884  24.993  1.00  8.58           C  
ATOM      3  C   ARG B   7      25.739  59.980  25.031  1.00  8.21           C  
ATOM      4  O   ARG B   7      25.115  59.712  24.002  1.00 11.30           O  
ATOM      5  CB  ARG B   7      26.493  62.323  24.783  1.00  9.10           C  
ATOM      6  CG  ARG B   7      25.640  62.539  23.542  1.00 10.66           C  
ATOM      7  CD  ARG B   7      26.479  63.330  22.563  1.00 16.50           C  
ATOM      8  NE  ARG B   7      26.070  64.698  22.619  1.00 15.14           N  
ATOM      9  CZ  ARG B   7      26.718  65.771  22.177  1.00 14.14           C  
ATOM     10  NH1 ARG B   7      27.927  65.757  21.607  1.00 14.22           N1+
ATOM     11  NH2 ARG B   7      26.048  66.895  22.227  1.00 14.84           N
"""
STRUCTURES[(63, -110, 43, 163)] = """\
ATOM      1  N   ARG B   7      27.813  60.399  23.912  1.00  8.49           N
ATOM      2  CA  ARG B   7      26.950  60.884  24.993  1.00  8.58           C
ATOM      3  C   ARG B   7      25.739  59.980  25.031  1.00  8.21           C
ATOM      4  O   ARG B   7      25.115  59.712  24.002  1.00 11.30           O
ATOM      5  CB  ARG B   7      26.494  62.323  24.783  1.00  9.10           C
ATOM      6  CG  ARG B   7      25.641  62.539  23.542  1.00 10.66           C
ATOM      7  CD  ARG B   7      26.479  63.330  22.563  1.00 16.50           C  
ATOM      8  NE  ARG B   7      27.799  62.784  22.555  1.00 15.14           N
ATOM      9  CZ  ARG B   7      28.756  62.954  21.648  1.00 14.14           C
ATOM     10  NH1 ARG B   7      28.633  63.670  20.526  1.00 14.22           N1+
ATOM     11  NH2 ARG B   7      29.925  62.453  21.961  1.00 14.84           N
"""
STRUCTURES[(63, -110, 43, -44)] = """\
ATOM      1  N   ARG B   7      27.813  60.399  23.912  1.00  8.49           N
ATOM      2  CA  ARG B   7      26.950  60.884  24.993  1.00  8.58           C
ATOM      3  C   ARG B   7      25.739  59.980  25.031  1.00  8.21           C
ATOM      4  O   ARG B   7      25.115  59.712  24.002  1.00 11.30           O
ATOM      5  CB  ARG B   7      26.493  62.324  24.782  1.00  9.10           C
ATOM      6  CG  ARG B   7      25.640  62.539  23.541  1.00 10.66           C
ATOM      7  CD  ARG B   7      26.479  63.330  22.561  1.00 16.50           C
ATOM      8  NE  ARG B   7      27.799  62.784  22.554  1.00 15.14           N
ATOM      9  CZ  ARG B   7      28.550  62.427  23.591  1.00 14.14           C
ATOM     10  NH1 ARG B   7      28.194  62.550  24.874  1.00 14.22           N1+
ATOM     11  NH2 ARG B   7      29.671  61.819  23.294  1.00 14.84           N
"""


class TestAccuracyMetrics(unittest.TestCase):
    def setUp(self):
        self.structures = {}
        parser = PDBParser()
        for angles, structure in STRUCTURES.items():
            with StringIO(structure) as struc_fh:
                self.structures[angles] = parser.get_structure(
                    repr(angles),
                    struc_fh)

    def assertListAlmostEqual(self, list_a, list_b, msg, places=10):
        self.assertEqual(
            len(list_a),
            len(list_b),
            "{}: lengths are equal".format(msg))
        for i, (a, b) in enumerate(zip(list_a, list_b)):
            self.assertAlmostEqual(
                a, b, places,
                "{}: at {} almost equal".format(msg, i))

    def pair(self, native, model):
        """Get a native/model pair representing the whole pipeline state."""
        return {"models": [{
            "native": self.structures[native],
            "model": self.structures[model]
        }]}
    def accuracy(self, data, acc_type):
        """Get the accuracy dict of the first model."""
        return data["models"][0]["sidechain_accuracy"][acc_type]

    def test_conditional(self):
        """Check conditional accuracy of a few structures with known χ angles."""
        cond = ConditionalAccuracy()

        # All angles equal
        accuracy = self.accuracy(
            cond.run(self.pair((-63, 130, 101, 163), (-63, 130, 101, 163))),
            "conditional")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0] * 4,
            "native vs. native")

        # First chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((-63, 130, 101, 163), (63, 130, 101, 163))),
            "conditional")
        self.assertListAlmostEqual(
            accuracy["ARG"], [0.0, None, None, None],
            "bad χ1")

        # Second chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, 130, 101, 163), (63, -110, 101, 163))),
            "conditional")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 0.0, None, None],
            "bad χ2")

        # Third chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, -110, 101, 163), (63, -110, 43, 163))),
            "conditional")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 1.0, 0.0, None],
            "bad χ3")

        # Fourth chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, -110, 43, 163), (63, -110, 43, -44))),
            "conditional")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 1.0, 1.0, 0.0],
            "bad χ4")

    def test_absolute(self):
        """Check absolute accuracy structures with known χ angles."""
        cond = AbsoluteAccuracy()

        # All angles equal
        accuracy = self.accuracy(
            cond.run(self.pair((-63, 130, 101, 163), (-63, 130, 101, 163))),
            "absolute")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0] * 4,
            "native vs. native")

        # First chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((-63, 130, 101, 163), (63, 130, 101, 163))),
            "absolute")
        self.assertListAlmostEqual(
            accuracy["ARG"], [0.0] * 4,
            "bad χ1")

        # Second chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, 130, 101, 163), (63, -110, 101, 163))),
            "absolute")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 0.0, 0.0, 0.0],
            "bad χ2")

        # Third chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, -110, 101, 163), (63, -110, 43, 163))),
            "absolute")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 1.0, 0.0, 0.0],
            "bad χ3")

        # Fourth chi angle differs
        accuracy = self.accuracy(
            cond.run(self.pair((63, -110, 43, 163), (63, -110, 43, -44))),
            "absolute")
        self.assertListAlmostEqual(
            accuracy["ARG"], [1.0, 1.0, 1.0, 0.0],
            "bad χ4")

    def test_rmsd(self):
        """Check RMSDs against known RMSDs."""
        rmsd = PerResidueRMSD()

        accuracy = self.accuracy(
            rmsd.run(self.pair((-63, 130, 101, 163), (-63, 130, 101, 163))),
            "rmsd")
        self.assertAlmostEqual(
            accuracy["ARG"], 0.000,
            msg="native vs. native", places=3)

        accuracy = self.accuracy(
            rmsd.run(self.pair((-63, 130, 101, 163), (63, 130, 101, 163))),
            "rmsd")
        self.assertAlmostEqual(
            accuracy["ARG"], 2.20725995798,
            msg="native vs. bad chi1", places=5)

        accuracy = self.accuracy(
            rmsd.run(self.pair((63, 130, 101, 163), (63, -110, 101, 163))),
            "rmsd")
        self.assertAlmostEqual(
            accuracy["ARG"], 1.87490952651,
            msg="native vs. bad chi2", places=5)
