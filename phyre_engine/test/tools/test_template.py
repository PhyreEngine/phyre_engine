"""Test the :py:mod:`phyre_engine.tools.template` module."""
import datetime
import io
from pathlib import Path
import shutil
import tempfile
import unittest
import Bio.PDB
import phyre_engine.test.data.minimal_template as minimal
from phyre_engine.tools.template import Template, TemplateDatabase

class TestTemplateDatabase(unittest.TestCase):
    """Test TemplateDatabase class."""

    _METADATA = {
        "deposition_date": datetime.date(1997, 12, 2),
        "release_date": datetime.date(1998, 12, 30),
        "last_update_date": datetime.date(2011, 7, 13),
        "method": "X-RAY DIFFRACTION",
        "resolution": 2.2,
        "organism_name": "Escherichia coli",
        "organism_id": 562,
        "title": (
            "ASPARAGINE SYNTHETASE MUTANT C51A, C315A COMPLEXED "
            "WITH L-ASPARAGINE AND AMP"),
        "descriptor": (
            "ASPARAGINE SYNTHETASE, L-ASPARAGINE, "
            "ADENOSINE MONOPHOSPHATE"),
    }

    @classmethod
    def setUpClass(cls):
        """
        Create an empty temporary database, load and write a template.
        """
        cls.file_root = Path(tempfile.mkdtemp("-templatedb", "phyreengine-"))
        cls.database = cls.file_root / "test.db"
        TemplateDatabase.create(str(cls.database))

        (cls.file_root / "2a" / "12as").mkdir(parents=True)
        mmcif_buf = io.StringIO(minimal.MINIMAL_MMCIF)
        mmcif_struc = Bio.PDB.MMCIFParser().get_structure("", mmcif_buf)
        cls.sample_template = Template.build("1MIN", "A", mmcif_struc[0]["A"])
        with (cls.file_root / "2a" / "12as" / "12as_A.pdb").open("w") as tout:
            cls.sample_template.write(tout)

    @classmethod
    def tearDownClass(cls):
        """Remove temporary database."""
        shutil.rmtree(str(cls.file_root))

    def test_add_retrieve_del_template(self):
        """Add PDB and template to database."""
        template_db = TemplateDatabase(str(self.database), str(self.file_root))
        template_db.add_pdb("12AS", self._METADATA)
        template_db.add_template("12AS", "A", self.sample_template)
        template_db.commit()

        template = template_db.get_template("12AS", "A")
        self.assertEqual(
            template.canonical_seq,
            self.sample_template.canonical_seq)
        self.assertEqual(
            template.canonical_indices,
            self.sample_template.canonical_indices)
        self.assertEqual(
            template.mapping,
            self.sample_template.mapping)

        template_db.del_template("12as", "A")
        template_db.commit()
        with self.assertRaises(template_db.TemplateNotFoundException):
            template_db.get_template("12as", "A")

    def test_retrieve_pdb(self):
        """Retrieve PDB metadata from database."""
        template_db = TemplateDatabase(str(self.database), str(self.file_root))
        self.assertEqual(template_db.get_pdb("12as"), self._METADATA)

    def test_add_del_pdb(self):
        """Add and delete a PDB entry."""
        template_db = TemplateDatabase(str(self.database), str(self.file_root))
        template_db.add_pdb("1del", self._METADATA)
        template_db.commit()
        self.assertEqual(template_db.get_pdb("1del"), self._METADATA)
        template_db.del_pdb("1del")
        template_db.commit()
        with self.assertRaises(template_db.PdbNotFoundException):
            template_db.get_pdb("1del")

    def test_get_nonexist(self):
        """
        TemplateNotFoundException raised when looking up nonexistent template.
        """
        template_db = TemplateDatabase(str(self.database), str(self.file_root))
        with self.assertRaises(TemplateDatabase.TemplateNotFoundException):
            template_db.get_template("1FOO", "A")


class TestTemplate(unittest.TestCase):
    """Test Template class"""

    def setUp(self):
        """Create a temporary directory for tests to use."""
        mmcif_buf = io.StringIO(minimal.MINIMAL_MMCIF)
        mmcif_struc = Bio.PDB.MMCIFParser().get_structure("", mmcif_buf)
        self.raw_chain = mmcif_struc[0]["A"]

    def _validate(self, template):
        """Validate 'template' against known data."""
        self.assertListEqual(
            template.mapping,
            minimal.ORIG_MAPPING)
        self.assertListEqual(
            template.canonical_indices,
            minimal.CANONICAL_SEQ_INDICES)
        self.assertEqual(
            template.canonical_seq,
            minimal.CANONICAL_SEQ)

    def test_build(self):
        """Ensure that a mapping of renumbered -> original AAs was written."""
        self._validate(Template.build(None, None, self.raw_chain))
