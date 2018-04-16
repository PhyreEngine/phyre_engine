import copy
import io
import json
import tempfile
import unittest
import unittest.mock
import phyre_engine.component
import phyre_engine.pipeline
import phyre_engine.component.db.db as db
import phyre_engine.tools.pdb as pdb
from phyre_engine.tools.template import Template
import phyre_engine.test
import phyre_engine.test.data
import phyre_engine.test.data.minimal_template as minimal
import phyre_engine.test.tools.test_pdb
from pathlib import Path
import shutil
import textwrap
import Bio.SeqIO
import Bio.PDB.PDBParser


@unittest.mock.patch("requests.Session")
class TestHTTPMap(unittest.TestCase):
    """Test HTTPMap component."""
    PIPELINE_STATE = {"list": [{"a": 1}, {"a": 2}]}

    class Dummy(phyre_engine.component.Component):
        ADDS = []
        REMOVES = []
        REQUIRED = ["a"]

        def __init__(self, tester):
            self.tester = tester

        def run(self, data, config=None, pipeline=None):
            self.tester.assertIn("session", data)
            data["a"] *= 2
            return data

    def test_map(self, session_mock):
        """Map operates as normal."""
        components = [self.Dummy(self)]
        pipeline = phyre_engine.pipeline.Pipeline(components=components)

        map_cpt = db.HTTPMap(field="list", pipeline=pipeline)
        result = map_cpt.run(self.PIPELINE_STATE)
        self.assertEqual(result, {"list": [{"a": 2}, {"a": 4}]})
        session_mock.assert_called_once()


@phyre_engine.test.requireFields("net_tests")
class TestStructureRetriever(unittest.TestCase):
    """Test STructureRetriever component."""
    _PIPE_STATE = {"PDB": "12as"}

    def test_retrieve(self):
        """Try and download some files."""
        types = ("pdb", "cif", db.StructureType.PDB, db.StructureType.MMCIF)

        for struc_type in types:
            with self.subTest("Getting structure", type=struc_type):
                with tempfile.TemporaryDirectory() as tmpdir:
                    retriever = db.StructureRetriever(struc_type, tmpdir)
                    retriever.run(self._PIPE_STATE)

                    if isinstance(struc_type, str):
                        suffix = struc_type
                    else:
                        suffix = struc_type.value

                    pdb_12as = Path(tmpdir, "2a/12as.{}.gz".format(suffix))
                    self.assertTrue(
                        pdb_12as.exists(),
                        "{!s} should exist".format(pdb_12as))

class TestPDBList(unittest.TestCase):
    """Test PDBList component."""

    _XML_DATA = r"""<?xml version='1.0' standalone='no' ?>
    <current>
      <PDB structureId="100D" />
      <PDB structureId="101D" />
    </current>
    """

    def test_file_parsing(self):
        """Read PDB IDs from a file."""
        with tempfile.NamedTemporaryFile("w") as xml_file:
            xml_file.write(self._XML_DATA)
            xml_file.flush()
            pdb_list = db.PDBList(xml_file.name)
            results = pdb_list.run({})
            self.assertListEqual(
                results["templates"],
                [{"PDB": "100D"}, {"PDB": "101D"}])

    @phyre_engine.test.requireFields("net_tests")
    def test_download(self):
        """Read PDB IDs from the RCSB website."""
        pdb_list = db.PDBList()
        results = pdb_list.run({})
        self.assertGreater(len(results["templates"]), 0)

class TestChainPDBBuilder(unittest.TestCase):
    """Test ChainPDBBuilder class"""

    def setUp(self):
        """Create a temporary directory for tests to use."""
        self.tmpdir = Path(tempfile.mkdtemp("-test", "chain-"))
        self.mmcif_dir = self.tmpdir / "mmcif"
        self.pdb_dir = self.tmpdir / "pdb"

        (self.mmcif_dir / "2a").mkdir(parents=True)
        self.pdb_dir.mkdir()

        with (self.mmcif_dir / "2a" / "12as.cif").open("w") as mmcif_out:
            mmcif_out.write(minimal.MINIMAL_MMCIF)

    def _build(self):
        builder = db.ChainPDBBuilder(str(self.mmcif_dir), str(self.pdb_dir))
        results = builder.run({"PDB": "12as"})
        return results

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.tmpdir))

    def test_create_chain(self):
        """Create a PDB file from a minimal MMCIF."""
        results = self._build()
        pdb_12as_A_path = Path(results[0]["structure"])
        self.assertTrue(
            pdb_12as_A_path.exists(),
            "PDB file containing chain was created.")

    def test_template(self):
        """Ensure that template sequence and mapping are correct."""
        results = self._build()
        self.assertEqual(
            results[0]["sequence"],
            minimal.CANONICAL_SEQ)
        self.assertEqual(
            results[0]["original_residues"],
            minimal.ORIG_MAPPING)
        self.assertEqual(
            results[0]["canonical_indices"],
            minimal.CANONICAL_SEQ_INDICES)


class TestBuildTemplate(unittest.TestCase):
    """Test BuildTemplate component."""

    def test_template(self):
        """Ensure that template metadata is correct."""
        with io.StringIO(minimal.MINIMAL_PDB) as structure_file:
            builder = db.BuildTemplate()
            results = builder.run({
                "PDB": "1abc",
                "chain": "A",
                "structure": structure_file,
                "original_residues": minimal.ORIG_MAPPING,
                "canonical_indices": minimal.CANONICAL_SEQ_INDICES,
                "sequence": minimal.CANONICAL_SEQ,
            })
            template = results["template_obj"]
            self.assertListEqual(
                template.mapping,
                minimal.ORIG_MAPPING)
            self.assertListEqual(
                template.canonical_indices,
                minimal.CANONICAL_SEQ_INDICES)
            self.assertEqual(
                template.canonical_seq,
                minimal.CANONICAL_SEQ)


class TestPDBSequence(unittest.TestCase):
    """Tests for the PDBSequence component."""

    def test_sequence(self):
        """Sequence should match atom_seq function."""
        seq_parser = db.PDBSequence()
        results = seq_parser.run({"template_obj": minimal.MINIMAL_TEMPLATE})
        self.assertEqual(
            results["sequence"], minimal.CANONICAL_SEQ,
            "Sequence read correctly")

class TestReduceExpand(unittest.TestCase):
    """Test the Reduce and Expand components."""

    _FULL_LIST = [
        {"a": 1, "b": 1, "x": 123},
        {"a": 1, "b": 1, "x": 456},

        {"a": 1, "b": 2, "x": 789},

        {"a": 2, "b": 1, "x": 321},
        {"a": 2, "b": 1, "x": 654},

        {"a": 2, "b": 2, "x": 987},
    ]

    _REDUCTION = [
        [
            {"a": 1, "b": 1, "x": 123},
            {"a": 1, "b": 1, "x": 456},
        ],
        [
            {"a": 1, "b": 2, "x": 789},
        ],
        [
            {"a": 2, "b": 1, "x": 321},
            {"a": 2, "b": 1, "x": 654},
        ],
        [
            {"a": 2, "b": 2, "x": 987},
        ]
    ]
    _REDUCED_LIST = [
        {"a": 1, "b": 1, "x": 123},
        {"a": 1, "b": 2, "x": 789},
        {"a": 2, "b": 1, "x": 321},
        {"a": 2, "b": 2, "x": 987},
    ]

    def test_reduce(self):
        """Reduce according to a tuple of fields."""
        reducer = db.Reduce("list", ("a", "b"))
        result = reducer.run({"foo": "bar", "list": self._FULL_LIST.copy()})
        self.assertEqual(result["foo"], "bar", "Untouched extra data")
        self.assertListEqual(result["list"], self._REDUCED_LIST)
        self.assertListEqual(result["reduction"], self._REDUCTION)

    def test_expand(self):
        """Expand according to a tuple of fields."""
        expander = db.Expand("list", ("a", "b"))
        result = expander.run({
            "list": self._REDUCED_LIST.copy(),
            "reduction": self._REDUCTION.copy()
        })
        self.assertListEqual(
            sorted(result["list"], key=lambda e: e["x"]),
            sorted(self._FULL_LIST, key=lambda e: e["x"]))

class TestAnnotateCATrace(unittest.TestCase):
    """Test AnnotateCATrace component."""

    _CA_TRACE = textwrap.dedent("""\
        REMARK 150 AA
        REMARK 156 [1, 2]
        REMARK 161 [[" ", 4, " "], [" ", 5, " "], [" ", 6, " "]]
        ATOM      2  CA  ALA A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      7  CA  ALA A   2      15.552  39.410  26.282  1.00  8.51
        ATOM     32  CB  ALA A   3      16.137  34.313  27.425  1.00  4.96
    """)

    _NON_CA_TRACE = textwrap.dedent("""\
        REMARK 150 AA
        REMARK 156 [1, 2]
        REMARK 161 [[" ", 4, " "], [" ", 5, " "], [" ", 6, " "]]
        ATOM      2  CA  ALA A   1      12.501  39.048  28.539  1.00 30.68
        ATOM      5  CB  ALA A   1      12.902  39.919  29.730  1.00 16.77
        ATOM      7  CA  ALA A   2      15.552  39.410  26.282  1.00  8.51
        ATOM     32  CB  ALA A   3      16.137  34.313  27.425  1.00  4.96
    """)

    @classmethod
    def setUpClass(cls):
        """Write mock structures to PDB files."""
        parser = Bio.PDB.PDBParser(QUIET=True)
        ca_pdb = parser.get_structure("CA", io.StringIO(cls._CA_TRACE))
        bb_pdb = parser.get_structure("BB", io.StringIO(cls._NON_CA_TRACE))
        cls.ca_template = unittest.mock.MagicMock(
            chain=ca_pdb[0]["A"],
            canonical_indices=[1, 2])
        cls.bb_template = unittest.mock.MagicMock(
            chain=bb_pdb[0]["A"],
            canonical_indices=[1, 2])

    def test_annotate_ca_trace(self):
        """Check that a CA trace is marked with ``ca_trace = True``."""
        annotator = db.AnnotateCATrace()
        result = annotator.run({"template_obj": self.ca_template})
        self.assertTrue(result["ca_trace"], "Annotated as CA trace")

    def test_annotate_non_ca_trace(self):
        """Check that a PDB with backbone is marked ``ca_trace = False.``"""
        annotator = db.AnnotateCATrace()
        result = annotator.run({"template_obj": self.bb_template})
        self.assertFalse(result["ca_trace"], "Annotated as not a CA trace")
