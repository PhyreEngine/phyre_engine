"""Test components in the :py:mod:`phyre_engine.component.hhsuite` module."""
import copy
import io
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import unittest.mock

import Bio.PDB

import phyre_engine.test.data.minimal_template as minimal
import phyre_engine.component.hhsuite as hhsuite
import phyre_engine.tools.template as template
import collections

class TestAlignmentToFasta(unittest.TestCase):
    """Test the AlignmentToFasta component."""

    def test_template_seq(self):
        """Build template_obj part of the alignment from i,j pairs."""
        template_obj = unittest.mock.MagicMock(spec=template.Template)
        template_obj.canonical_seq = "ABCDEFGH"
        query_seq = "JKLMNOP"
        aln = [(1, 1), (3, 5)]

        template_seq = hhsuite.AlignmentToFasta.build_template_aln(
            query_seq, aln, template_obj)
        self.assertEqual(template_seq, "A-E----")

        # Test formatting of the fasta using different fields.
        pipe_state = {"sequence": query_seq, "name": "queryname"}
        hit = {"name": "templatename"}

        aln2fasta = hhsuite.AlignmentToFasta()
        self.assertEqual(
            aln2fasta.fasta(pipe_state, hit, template_seq),
            ">Query\nJKLMNOP\n>Template\nA-E----\n")

        aln2fasta = hhsuite.AlignmentToFasta(q_name="name", t_name="name")
        self.assertEqual(
            aln2fasta.fasta(pipe_state, hit, template_seq),
            ">queryname\nJKLMNOP\n>templatename\nA-E----\n")

class TestFastaParser(unittest.TestCase):
    """Test the FastaParser component."""

    FASTA = textwrap.dedent("""\
    No 1
    >Query
    AB-CDE
    >Template
    -A-A--

    No 2
    >Query
    AB
    >foo
    -A
    >Template
    GG
    """)

    PIPELINE = {
        "sequence": "XXABCDEFXX",
        "name": "Query",
        "templates": [
            {"query_range": range(3, 7)},
            {"query_range": range(3, 4)},
        ]
    }

    Pair = collections.namedtuple("Pair", "i j probab")
    ALIGNMENTS = [
        [Pair(4, 1, 0.5), Pair(5, 2, 0.4)],
        [Pair(5, 4, 0.3)],
    ]

    def setUp(self):
        """Create a copy of the pipeline."""
        self.pipeline = copy.deepcopy(self.PIPELINE)
        self.pipeline["pairwise_fasta"] = io.StringIO(self.FASTA)

    def test_parser(self):
        """Parse a FASTA alignment file into the fasta_alignment field."""
        parser = hhsuite.FastaParser()
        results = parser.run(self.pipeline)
        self.assertEqual(
            results["templates"][0]["sequence_alignments"], {
                "sequence": "---A-A-----",
                "query": "XXAB-CDEFXX"
        })

        self.assertEqual(
            results["templates"][1]["sequence_alignments"], {
                "foo": "---A------",
                "sequence": "--GG------",
                "query": "XXABCDEFXX"
        })

    def test_ignore(self):
        """Ignored sequences should be removed."""
        parser = hhsuite.FastaParser(ignore={"foo"})
        results = parser.run(self.pipeline)
        self.assertNotIn("foo", results["templates"][1]["sequence_alignments"])

    def test_confidences(self):
        """Test generation of confidence string from alignment pairs."""

        # Add alignments to pipeline
        for hit, aln in zip(self.pipeline["templates"], self.ALIGNMENTS):
            hit["alignment"] = aln

        parser = hhsuite.FastaParser()
        results = parser.run(self.pipeline)
        self.assertEqual(
            results["templates"][0]["sequence_alignments"]["confidence"],
            "---5-4-----")
        self.assertEqual(
            results["templates"][1]["sequence_alignments"]["confidence"],
            "----3-----")

class TestA3MSSParser(unittest.TestCase):
    """Test the A3MSSParser."""

    SAMPLE_A3M = textwrap.dedent("""\
    >ss_dssp DSSP Secondary structure
    CCCTHH
    >ss_pred Secondary structure from PSIPRED
    CCCCHH
    >ss_conf SS Confidence
    888889
    >Query
    AAAAGA
    """)

    def test_parser(self):
        """Test A3M parser."""
        parser = hhsuite.A3MSSParser()

        mock_open = unittest.mock.mock_open(read_data=self.SAMPLE_A3M)
        # Bit of a hack to work around a bug: mock_open doesn't provide an
        # iterator on its returned MagicMock.
        mock_open.return_value.__iter__ = lambda self: iter(self.readline, '')
        open_function = ".".join([parser.__module__, "open"])

        with unittest.mock.patch(open_function, mock_open):
            results = parser.run({"a3m": "query.a3m"})
            self.assertEqual(
                results["secondary_structure_sequence"], {
                    "ss_dssp": "CCCTHH",
                    "ss_pred": "CCCCHH",
                    "ss_conf": "888889"
                })


class TestAddDssp(unittest.TestCase):
    """Test AddDssp component."""

    BEFORE_A3M = textwrap.dedent("""\
    >Query
    {query}
    """).format(query=minimal.CANONICAL_SEQ)

    # Note the gaps in aa_dssp and ss_dssp. This is so we can test the gapping.
    AFTER_A3M = textwrap.dedent("""\
    >aa_dssp
    {aa_dssp}
    >ss_dssp
    HH-ECC
    >Query
    {query}
    """).format(
        aa_dssp=minimal.CANONICAL_SEQ[:2] + "-" + minimal.CANONICAL_SEQ[3:],
        query=minimal.CANONICAL_SEQ)

    SEC_STRUC = [
        {"res_id": 1, "assigned": "H"},
        {"res_id": 2, "assigned": "H"},
        # Note absence of residue 3. It should be gapped.
        {"res_id": 4, "assigned": "E"},
        {"res_id": 5, "assigned": "C"},
        {"res_id": 6, "assigned": "C"},
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
            template_obj = template.Template.build(
                "1MIN", "A", template_structure[0]["A"])
            with structure.open("w") as template_out:
                template_obj.write(template_out)

        # Write a3m to file.
        with a3m.open("w") as a3m_out:
            a3m_out.write(self.BEFORE_A3M)

        self.pipeline = {
            "template_obj": template_obj,
            "a3m": str(a3m),
            "secondary_structure": {"dssp": self.SEC_STRUC}
        }

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(str(self.dir))

    def test_add_dssp(self):
        """Add DSSP info to an a3m file."""
        add_dssp = hhsuite.AddDssp()
        results = add_dssp.run(self.pipeline)
        with open(results["a3m"], "r") as a3m_in:
            self.assertEqual(a3m_in.read(), self.AFTER_A3M)

        # Run again: should not alter a3m contents.
        results = add_dssp.run(results)
        with open(results["a3m"], "r") as a3m_in:
            self.assertEqual(a3m_in.read(), self.AFTER_A3M)


class HHToolTest(unittest.TestCase):
    """
    Test subclasses of HHSuiteTool.

    This class contains some fields to use for testing, as well as a `tool`
    method to create an instance of a tool using those fields, and a
    `verify_common` method to verify those fields.
    """

    BIN_DIR = "hh/bin"
    FLAGS = ["x", "y", "z"]
    HHLIB = "hhlib"
    OPTIONS = {"foo": "bar", "qux": "baz"}
    DATABASE = "test_db"

    def tool(self, tool_cls, *args, options={}, **kwargs):
        """
        Create an instance of `tool_cls` using some default flags and options.
        The `options` hash will be merged into `self.OPTIONS`.
        """
        tool_options = self.OPTIONS.copy()
        tool_options.update(options)

        tool_instance = tool_cls(
            *args,
            flags=self.FLAGS,
            bin_dir=self.BIN_DIR,
            HHLIB=self.HHLIB,
            options=tool_options,
            **kwargs)
        tool_instance.tool = unittest.mock.MagicMock()
        return tool_instance

    def verify_common(self, tool_name, tool_instance):
        """
        Verify that the arguments `tool_args` (passed to a mocked
        `phyre_engine.tools.external.ExternalTool.__call__` method) match the
        expected testing data.
        """
        pos_args, kw_args = tool_instance.tool.call_args
        self.assertEqual(pos_args[0], (self.BIN_DIR, tool_name))
        self.assertEqual(kw_args["flags"], self.FLAGS)

        # Existing options were retained
        for flag, expected in self.OPTIONS.items():
            self.assertEqual(kw_args["options"][flag], expected)


@unittest.mock.patch("subprocess.run")
class TestHHBlits(HHToolTest):
    """Test HHBlits component."""

    DATABASE = "test_db"
    SEQ_NAME = "seq_name"
    SEQUENCE = "AAA"
    IN_HHM_FILE = "in.hhm"
    IN_A3M_FILE = "in.a3m"
    OUT_A3M_FILE = "out.a3m"
    OUT_REPORT = "report.hhr"

    def verify_common(self, tool_name, tool_instance):
        """Verify database and output files."""
        super().verify_common(tool_name, tool_instance)
        pos_args, kw_args = tool_instance.tool.call_args

        self.assertEqual(kw_args["options"]["oa3m"], self.OUT_A3M_FILE)
        self.assertEqual(kw_args["options"]["output"], self.OUT_REPORT)
        self.assertEqual(kw_args["options"]["database"], self.DATABASE)

    def tool(self, *args, **kwargs):
        """
        Add database to constructor, and 'oa3m' and 'output' keys to options
        and call superclass.
        """
        options = {"oa3m": self.OUT_A3M_FILE, "output": self.OUT_REPORT}
        return super().tool(hhsuite.HHBlits, self.DATABASE, *args,
                            options=options, **kwargs)

    def test_input_type_seq(self, _run_mock):
        """Check command line built by HHBlits when given a sequence."""
        hhblits = self.tool(input_type=hhsuite.QueryType.SEQUENCE)
        self.assertEqual(set(hhblits.REQUIRED), set(["name", "sequence"]))
        hhblits.run({"sequence": self.SEQUENCE, "name": self.SEQ_NAME})
        self.verify_common("hhblits", hhblits)

        _, kw_args = hhblits.tool.call_args
        self.assertIn("input", kw_args["options"])

    def test_input_type_a3m(self, _run_mock):
        """Check command line built by HHBlits when given an a3m file."""
        hhblits = self.tool(input_type=hhsuite.QueryType.A3M)
        self.assertEqual(hhblits.REQUIRED, ["a3m"])
        hhblits.run({"a3m": self.IN_A3M_FILE})
        self.verify_common("hhblits", hhblits)
        _, kw_args = hhblits.tool.call_args
        self.assertEqual(kw_args["options"]["input"], self.IN_A3M_FILE)

    def test_input_type_hhm(self, _run_mock):
        """Check command line built by HHBlits when given an hhm file."""
        hhblits = self.tool(input_type=hhsuite.QueryType.HMM)
        self.assertEqual(hhblits.REQUIRED, ["hhm"])
        hhblits.run({"hhm": self.IN_HHM_FILE})
        self.verify_common("hhblits", hhblits)
        _, kw_args = hhblits.tool.call_args
        self.assertEqual(kw_args["options"]["input"], self.IN_HHM_FILE)


@unittest.mock.patch("subprocess.run")
class TestHHMake(HHToolTest):
    """Test HHMake component."""
    IN_A3M_FILE = "in.a3m"
    OUT_HHM_FILE = "out.hhm"

    def test_run(self, _run_mock):
        """Test arguments given to HHMake tool."""
        hhmake = self.tool(hhsuite.HHMake,
                           options={"output": self.OUT_HHM_FILE})
        hhmake.run({"a3m": self.IN_A3M_FILE})
        self.verify_common("hhmake", hhmake)

        _, kw_args = hhmake.tool.call_args
        self.assertEqual(kw_args["options"]["input"], self.IN_A3M_FILE)


@unittest.mock.patch("subprocess.run")
class TestCSTranslate(HHToolTest):
    """Test CSTranslate component."""
    IN_A3M_FILE = "in.a3m"
    OUT_CS219_FILE = "out.cs219"

    def test_run(self, _run_mock):
        """Test arguments given to CSTranslate tool."""
        cstranslate = self.tool(hhsuite.CSTranslate,
                           options={"outfile": self.OUT_CS219_FILE})
        cstranslate.run({"a3m": self.IN_A3M_FILE})
        self.verify_common("cstranslate", cstranslate)

        _, kw_args = cstranslate.tool.call_args
        self.assertEqual(kw_args["options"]["infile"], self.IN_A3M_FILE)
