"""Test components in the :py:mod:`phyre_engine.component.disorder` module."""
import textwrap
import unittest
import phyre_engine.test
import phyre_engine.component.disorder as disorder

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
