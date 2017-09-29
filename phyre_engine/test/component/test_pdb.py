"""Test the components in the :py:mod:`phyre_engine.component.pdb` module."""
import io
import tempfile
import unittest
import Bio.PDB
import phyre_engine.component.pdb
import phyre_engine.tools.pdb
import phyre_engine.test.data.minimal_template as minimal
from unittest.mock import MagicMock, sentinel

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

class TestConvertToTemplate(unittest.TestCase):
    """Test ConvertToTemplate component."""

    @staticmethod
    def structure():
        """Get a minimal structure as a monomer."""
        struc_file = io.StringIO(minimal.MINIMAL_PDB)
        parser = Bio.PDB.PDBParser(QUIET=True)
        return parser.get_structure("", struc_file)[0]["A"]

    def test_converter(self):
        """Convert minimal PDB to template."""
        # pylint: disable=protected-access
        conv = phyre_engine.component.pdb.ConvertToTemplate()
        fh_mock = MagicMock()
        conv._open_structure = MagicMock(return_value=(fh_mock, sentinel.name))

        input_structure = self.structure()
        results = conv.run({"structure_obj": input_structure})
        self.assertIn("structure", results, "structure member added")
        self.assertEqual(results["structure"], sentinel.name)
        self.assertIsNot(results["structure_obj"], input_structure,
                         "Altered structure object returned")

        # Something was written to the new file.
        fh_mock.write.assert_called()

class TestTemplateMapping(unittest.TestCase):
    """Test TemplateMapping component."""

    def test_template_mapping(self):
        """Read template mapping from template."""
        template_in = io.StringIO(minimal.MINIMAL_PDB)
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        template_struc = pdb_parser.get_structure("", template_in)[0]["A"]
        template_buf = io.StringIO()
        template = phyre_engine.tools.template.Template.build(template_struc)
        template.write(template_buf)
        template_buf.seek(0)

        mapper = phyre_engine.component.pdb.TemplateMapping()
        results = mapper.run({"structure": template_buf})
        self.assertEqual(results["residue_mapping"], template.mapping)
