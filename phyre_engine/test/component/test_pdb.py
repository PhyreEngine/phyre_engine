"""Test the components in the :py:mod:`phyre_engine.component.pdb` module."""
import io
import tempfile
import unittest
import Bio.PDB
import phyre_engine.component.pdb
import phyre_engine.tools.pdb
import phyre_engine.test.data.minimal_template as minimal
from ast import parse

class TestPDBSeq(unittest.TestCase):
    """Test conversion of PDB structure to sequence using PDBSeq."""
    def setUp(self):
        """Create PDB structure to read."""
        self.structure = Bio.PDB.PDBParser(QUIET=True).get_structure(
            "", io.StringIO(minimal.MINIMAL_PDB))

    def test_pdbseq(self):
        """Read sequence from structure."""
        pdbseq_cmpt = phyre_engine.component.pdb.PDBSeq()
        results = pdbseq_cmpt.run({"structure_obj": self.structure})
        self.assertEqual(results["sequence"], minimal.CANONICAL_SEQ)

class TestReadStructure(unittest.TestCase):
    """Test ReadStructure component."""

    def _test_structure(self, structure_string, structure_type):
        """Write a structure and then parse it."""
        struc_file = io.StringIO(structure_string)
        reader = phyre_engine.component.pdb.ReadStructure()
        results = reader.run({"structure": struc_file})
        self.assertEqual(results["structure_type"], structure_type)
        atom_seq, _ = phyre_engine.tools.pdb.atom_seq(
            results["structure_obj"].get_residues())
        self.assertEqual(atom_seq, minimal.CANONICAL_SEQ)

    def test_pdb(self):
        """Read a PDB file."""
        self._test_structure(minimal.MINIMAL_PDB, "PDB")

    def test_mmcif(self):
        """Read an mmCIF file."""
        self._test_structure(minimal.MINIMAL_MMCIF, "MMCIF")

    def test_garbage(self):
        """Exception raised when reading garbage."""
        with self.assertRaises(phyre_engine.component.pdb.UnknownStructureTypeError):
            self._test_structure("FOO BAR", "PDB")

class TestConvertToMonomer(unittest.TestCase):
    """Test ConvertToMonomer class."""

    @staticmethod
    def structure():
        """Get a minimal structure."""
        struc_file = io.StringIO(minimal.MINIMAL_PDB)
        parser = Bio.PDB.PDBParser(QUIET=True)
        return parser.get_structure("", struc_file)

    def test_single_monomer(self):
        """A single-monomer structure is converted to a chain."""
        conv = phyre_engine.component.pdb.ConvertToMonomer()
        structure = self.structure()
        del structure[0]["D"]
        results = conv.run({"structure_obj": structure})
        self.assertIsInstance(
            results["structure_obj"],
            Bio.PDB.Chain.Chain)

    def test_multi_monomer(self):
        """Multiple monomers should cause an exception."""
        conv = phyre_engine.component.pdb.ConvertToMonomer()
        structure = self.structure()
        structure[0].add(Bio.PDB.Chain.Chain("B"))
        exception_class = phyre_engine.component.pdb.TooManyChainsError
        with self.assertRaises(exception_class):
            results = conv.run({"structure_obj": structure})
