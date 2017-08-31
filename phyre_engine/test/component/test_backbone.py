"""
Test backbone construction tools in :py:mod:`phyre_engine.component.backbone.`
"""
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import Bio.PDB.PDBParser
import phyre_engine.component.backbone as backbone
import phyre_engine.test
import phyre_engine.tools.pdb as pdb


class BaseBackboneTest(unittest.TestCase):
    """Base test class for writing temporary files and checking results."""

    _CA_TRACE = textwrap.dedent("""\
        ATOM      2  CA  ALA A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      7  CA  TYR A   2      15.552  39.410  26.282  1.00  8.51
        ATOM     19  CA  ILE A   6      17.791  39.281  29.375  1.00 23.27
        ATOM     27  CA  ALA A   7      16.004  36.186  30.742  1.00  5.20
        ATOM     32  CA  LYS A   8      16.137  34.313  27.425  1.00  4.96
        ATOM     41  CA  GLN A   9      19.794  35.327  26.885  1.00  8.15
        ATOM     50  CA  ARG A  10      20.706  33.656  30.141  1.00 19.12
    """)

    def setUp(self):
        """Create a temporary directory for input/output files."""
        self.tmpdir = pathlib.Path(tempfile.mkdtemp())
        self.ca_trace = self.tmpdir / "ca_trace.pdb"
        with self.ca_trace.open("w") as pdb_out:
            pdb_out.write(self._CA_TRACE)
            pdb_out.flush()

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(str(self.tmpdir))

    def has_backbone(self, pdb_file):
        """Check that a given PDB file has a backbone."""
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("bb", pdb_file)
        for res in structure.get_residues():
            with self.subTest("Contains non-CA atoms", residue=res.get_id()):
                self.assertGreater(
                    len(res), 1,
                    "More than 1 atom in residue.")

@phyre_engine.test.requireFields("pd2_ca2main", ["tools"])
class TestPD2CA2main(BaseBackboneTest):
    """Construct a backbone with pd2_ca2main."""

    def has_remark(self, template):
        """Ensure REMARK 999 was written containing backbone provenance."""
        with open(template, "r") as template_in:
            self.assertGreater(len(pdb.read_remark(template_in, 999)), 1)

    def test_reconstruction(self):
        """Reconstruct a backbone."""
        config = phyre_engine.test.config["tools"]["pd2_ca2main"]
        constructor = backbone.PD2CA2main(
            flags=["ca2main:new_fixed_ca"],
            options={"ca2main:bb_min_steps": 1000},
            **config)
        result = constructor.run({"structure": str(self.ca_trace)})
        self.has_backbone(result["structure"])
        self.has_remark(result["structure"])
