"""Test components in the :py:mod:`phyre_engine.component.disorder` module."""
import io
import os
import pathlib
import shutil
import tempfile
import textwrap
import unittest
import phyre_engine.test
import phyre_engine.component.disorder as disorder

DATA_DIR = pathlib.Path(phyre_engine.test.__file__).parent / "data"

class TestMobiDBLite(unittest.TestCase):
    """Test MobiDBLite component."""

    _MOBIDB_OUTPUT = textwrap.dedent(r"""
    {
        "acc": "sp|P49137|MAPK2_HUMAN",
        "consensus": "DDDSSS",
        "n_pred": 8,
        "p": [
            0.88,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0
        ],
        "pred": "mobidb-lite",
        "regions": [
          [
            1,
            3
          ]
        ]
    }
    """)

    _SAMPLE_SEQ = (
        "MLSNSQGQSPPVPFPAPAPPPQPPTPALPHPPAQPPPPPPQQFPQFHVKSGLQIKKNAII"
        "DDYKVTSQVLGLGINGKVLQIFNKRTQEKFALKMLQDCPKARREVELHWRASQCPHIVRI"
        "VDVYENLYAGRKCLLIVMECLDGGELFSRIQDRGDQAFTEREASEIMKSIGEAIQYLHSI"
        "NIAHRDVKPENLLYTSKRPNAILKLTDFGFAKETTSHNSLTTPCYTPYYVAPEVLGPEKY"
        "DKSCDMWSLGVIMYILLCGYPPFYSNHGLAISPGMKTRIRMGQYEFPNPEWSEVSEEVKM"
        "LIRNLLKTEPTQRMTITEFMNHPWIMQSTKVPQTPLHTSRVLKEDKERWEDVKEEMTSAL"
        "ATMRVDYEQIKIKKIEDASNPLLLKRRKKARALEAAALAH")

    def test_parser(self):
        """Parse known output."""
        results = disorder.MobiDBLite.parse_results(self._MOBIDB_OUTPUT)
        self.assertListEqual(
            results, [
                {"assigned": "D", "confidence": {"D": 0.88, "S": 0.12}},
                {"assigned": "D", "confidence": {"D": 1.00, "S": 0.00}},
                {"assigned": "D", "confidence": {"D": 1.00, "S": 0.00}},
                {"assigned": "S", "confidence": {"D": 0.00, "S": 1.00}},
                {"assigned": "S", "confidence": {"D": 0.00, "S": 1.00}},
                {"assigned": "S", "confidence": {"D": 0.00, "S": 1.00}}])

    @phyre_engine.test.requireFields("mobidb-lite", ["tools"])
    def test_mobidb_lite(self):
        """Run MobiDB lite and parse output."""
        mdb_conf = phyre_engine.test.config["tools"]["mobidb-lite"]
        mdb = disorder.MobiDBLite(**mdb_conf)
        results = mdb.run({"sequence": self._SAMPLE_SEQ})

        self.assertListEqual(
            results["disorder"]["mobidb-lite"][0:8], [
                {"assigned": "D", "confidence": {"D": 0.88, "S": 0.12}},
                {"assigned": "D", "confidence": {"D": 1.0, "S": 0.0}},
                {"assigned": "D", "confidence": {"D": 1.0, "S": 0.0}},
                {"assigned": "D", "confidence": {"D": 1.0, "S": 0.0}},
                {"assigned": "D", "confidence": {"D": 0.88, "S": 0.12}},
                {"assigned": "D", "confidence": {"D": 1.0, "S": 0.0}},
                {"assigned": "D", "confidence": {"D": 0.88, "S": 0.12}},
                {"assigned": "D", "confidence": {"D": 1.0, "S": 0.0}}])


class TestDisopred(unittest.TestCase):
    """Test Disopred component."""

    SAMPLE_OUTPUT = textwrap.dedent("""\
    #         ----- DISOPRED version 3.1 -----
    # Disordered residues are marked with asterisks (*)
    #    Ordered residues are marked with dots (.)
        1 M * 0.78
        2 K * 0.62
        3 T . 0.45
        4 A . 0.37
        5 Y . 0.20
    """)

    EXPECTED_OUTPUT = [
        {"assigned": "D", "confidence": {"S": 0.22, "D": 0.78}},
        {"assigned": "D", "confidence": {"S": 0.38, "D": 0.62}},
        {"assigned": "S", "confidence": {"S": 0.55, "D": 0.45}},
        {"assigned": "S", "confidence": {"S": 0.63, "D": 0.37}},
        {"assigned": "S", "confidence": {"S": 0.80, "D": 0.20}},
    ]

    def setUp(self):
        """Move to a temporary directory."""
        self.orig_dir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

    def tearDown(self):
        """Move back to original directory and remove tempdir."""
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir)

    def test_parser(self):
        """Test disopred output parser."""
        with io.StringIO(self.SAMPLE_OUTPUT) as result_stream:
            got_output = disorder.Disopred.parse_results(result_stream)
            for i, (expected, got) in enumerate(zip(self.EXPECTED_OUTPUT,
                                                    got_output)):
                with self.subTest("Residue {}".format(i+1)):
                    self.assertEqual(expected["assigned"], got["assigned"])
                    self.assertAlmostEqual(
                        expected["confidence"]["S"],
                        got["confidence"]["S"],
                        places=2)
                    self.assertAlmostEqual(
                        expected["confidence"]["D"],
                        got["confidence"]["D"],
                        places=2)

    @phyre_engine.test.requireFields("disopred", ["tools"])
    def test_disopred(self):
        """Run disopred on a prebuilt MTX file."""

        args = phyre_engine.test.config["tools"]["disopred"]
        mtx_file = str(DATA_DIR / "pssm" / "12as_A.mtx")
        disopred = disorder.Disopred(**args)
        results = disopred.run({"pssm": {"mtx": mtx_file}})
        self.assertGreater(len(results["disorder"]["disopred"]), 0)

