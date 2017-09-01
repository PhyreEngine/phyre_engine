"""Test the :py:mod:`phyre_engine.tools.template` module."""
import io
import unittest
import Bio.PDB
import phyre_engine.test.data.minimal_template as minimal
from phyre_engine.tools.template import Template

class TestTemplate(unittest.TestCase):
    """Test Template class"""

    def setUp(self):
        """Create a temporary directory for tests to use."""
        mmcif_buf = io.StringIO(minimal.MINIMAL_MMCIF)
        mmcif_struc = Bio.PDB.MMCIFParser().get_structure("", mmcif_buf)
        self.raw_chain = mmcif_struc[0]["A"]

    def _load(self):
        """Write and then load a Template."""
        template = Template.build(self.raw_chain)
        buffer = io.StringIO()
        template.write(buffer)
        buffer.seek(0)
        return Template.load(buffer)

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
        self._validate(Template.build(self.raw_chain))

    def test_load(self):
        """Test that we can write and load a template."""
        self._validate(self._load())
