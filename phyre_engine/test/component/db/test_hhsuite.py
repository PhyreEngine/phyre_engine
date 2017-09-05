"""
Test components in the :py:mod:`phyre_engine.component.db.hhsuite` module.
"""
import io
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import Bio.PDB
import phyre_engine.component.db.hhsuite as hhsuite
import phyre_engine.test.data.minimal_template as minimal
from phyre_engine.tools.template import Template

class TestAddDSSP(unittest.TestCase):
    """Test AddDSSP component."""

    BEFORE_A3M = textwrap.dedent("""\
    >Query
    {query}
    """).format(query=minimal.CANONICAL_SEQ)

    AFTER_A3M = textwrap.dedent("""\
    >aa_dssp
    {query}
    >ss_dssp
    HHEECC
    >Query
    {query}
    """).format(query=minimal.CANONICAL_SEQ)

    SEC_STRUC = [
        (1, "H"),
        (2, "H"),
        (3, "E"),
        (4, "E"),
        (5, "C"),
        (6, "C"),
    ]

    def setUp(self):
        """Set up a simple pipeline in a temporary directory."""
        self.dir = pathlib.Path(tempfile.mkdtemp(prefix="test-hhsuite-"))
        a3m = self.dir / "query.a3m"
        structure = self.dir / "query.pdb"

        # Write template to PDB file
        with io.StringIO(minimal.MINIMAL_MMCIF) as mmcif_buf:
            parser = Bio.PDB.MMCIFParser()
            template_structure = parser.get_structure("", mmcif_buf)
            template = Template.build(template_structure[0]["A"])
            with structure.open("w") as template_out:
                template.write(template_out)

        # Write a3m to file.
        with a3m.open("w") as a3m_out:
            a3m_out.write(self.BEFORE_A3M)

        self.pipeline = {
            "structure": str(structure),
            "a3m": str(a3m),
            "secondary_structure": self.SEC_STRUC
        }

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.dir))

    def test_add_dssp(self):
        """Add DSSP info to an a3m file."""
        add_dssp = hhsuite.AddDSSP()
        results = add_dssp.run(self.pipeline)
        with open(results["a3m"], "r") as a3m_in:
            self.assertEqual(a3m_in.read(), self.AFTER_A3M)

        # Run again: should not alter a3m contents.
        results = add_dssp.run(results)
        with open(results["a3m"], "r") as a3m_in:
            self.assertEqual(a3m_in.read(), self.AFTER_A3M)