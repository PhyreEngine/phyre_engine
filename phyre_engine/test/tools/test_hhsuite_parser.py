import unittest
from math import log10
import phyre_engine.test.data
import pathlib
from phyre_engine.tools.hhsuite.parser import Report, Tabular

DATA_DIR = pathlib.Path(phyre_engine.test.data.__file__).parent / "hhsuite"

class TestReportParser(unittest.TestCase):
    """Test a report parser."""

    _EXPECTED_HIT_16 = Report.Hit(
        name=("tr|S9CF61|S9CF61_STRAG Uncharacterized protein "
              "OS=Streptococcus agalactiae FSL S3-105 GN=SAG0023_02635 "
              "PE=3 SV=1"),
        rank=17,
        prob=100.0,
        evalue=1.1e-44,
        pvalue=5.4e-50,
        score=307.6,
        ss=0.0,
        cols=158,
        query_range=range(170, 327),
        template_range=range(1, 158),
        template_matches=158,
        identities=80,
         similarity=1.253,
        sum_probs=155.3,
        template_neff=1.935)

    _EXPECTED_HIT_147 = Report.Hit(
        name=("tr|W0U5P1|W0U5P1_9FIRM Aspartyl/glutamyl-tRNA(Asn/Gln) "
              "amidotransferase,B subunit OS=Ruminococcus bicirculans "
              "GN=gatB PE=3 SV=1"),
        rank=148,
        prob=74.7,
        evalue=5.0,
        pvalue=2.6e-5,
        score=48.2,
        ss=0.0,
        cols=261,
        query_range=range(23, 318),
        template_range=range(675, 951),
        template_matches=957,
        identities=18,
        similarity=0.221,
        sum_probs=146.3,
        template_neff=2.645)

    def _test_hits_equal(self, hit1, hit2):
        self.assertEqual(hit1.rank, hit2.rank)
        self.assertAlmostEqual(hit1.prob, hit2.prob, 1)
        self.assertAlmostEqual(log10(hit1.evalue), log10(hit2.evalue), 3)
        self.assertAlmostEqual(log10(hit1.pvalue), log10(hit2.pvalue), 3)
        self.assertAlmostEqual(hit1.score, hit2.score, 1)
        self.assertAlmostEqual(hit1.ss, hit2.ss, 1)
        self.assertEqual(hit1.cols, hit2.cols)
        self.assertEqual(hit1.query_range.start, hit2.query_range.start)
        self.assertEqual(hit1.query_range.stop, hit2.query_range.stop)
        self.assertEqual(hit1.template_range.start, hit2.template_range.start)
        self.assertEqual(hit1.template_range.stop, hit2.template_range.stop)
        self.assertEqual(hit1.template_matches, hit2.template_matches)
        self.assertEqual(hit1.name, hit2.name)
        self.assertEqual(hit1.identities, hit2.identities)
        self.assertAlmostEqual(hit1.similarity, hit2.similarity, 1)
        self.assertAlmostEqual(hit1.sum_probs, hit2.sum_probs, 1)
        self.assertAlmostEqual(hit1.template_neff, hit2.template_neff, 3)

    def test_report(self):
        """Parse a sample report"""
        report = Report(str(DATA_DIR / "test.hhr"))
        self.assertEqual(report.summary.query, "d12asa_")
        self.assertEqual(report.summary.match_cols, 327)

        self.assertEqual(report.summary.neff, 4.17492)
        self.assertEqual(report.summary.num_searched, 428)
        self.assertEqual(report.summary.date, "Tue Oct 25 11:04:48 2016")
        self.assertEqual(report.summary.command, (
            "hhblits -d "
            "/data/databases/uniprot20_2016_02/uniprot20_2016_02 "
            "-oa3m query.a3m -o report.hhr "
            "-cpu 20 -i /tmp/tmpqa0qxl3b.fasta"))
        self.assertTupleEqual(report.summary.num_seqs, (177, 221))

        self.assertEqual(len(report.hits), 273)

        #Pick a couple of arbitrary hits
        self._test_hits_equal(report.hits[16], self._EXPECTED_HIT_16)
        self._test_hits_equal(report.hits[147], self._EXPECTED_HIT_147)

class TestTabularParser(unittest.TestCase):
    """Test a tabular parser."""

    # Correct pairs for hit 0
    _EXPECTED_HIT_0 = Tabular.Hit(
        ("3m4p_A Ehasnrs, asparaginyl-tRNA synthetase, putative; "
         "aminoacyl-tRNA synthetase, tRNA ligase, AARS, translation, "
         "ATP-binding, nucleotide-binding; HET: 4AD; 2.83A "
         "{Entamoeba histolytica} PDB: 3m4q_A"),
        [
            Tabular.ResiduePair(9, 159, 0.10, 0.00, 0.3565, "H"),
            Tabular.ResiduePair(10, 160, 0.35, 0.00, 0.3835, "H"),
            Tabular.ResiduePair(11, 161, 0.15, 0.00, 0.4059, "H"),
            Tabular.ResiduePair(12, 162, 1.39, 0.00, 0.4281, "H"),
            Tabular.ResiduePair(13, 163, -0.19, 0.00, 0.4436, "H")
        ]
    )

    # Correct pairs for hit 1
    _EXPECTED_HIT_1 = Tabular.Hit(
        ("1h4v_B Histidyl-tRNA synthetase; class IIA aminoacyl-tRNA synthetase,"
         " ATP + L-histidine tRNA(His)-> AMP + PPI + L-histidyl-tRNA(His); 2.4A"
         " {Thermus thermophilus} SCOP: c.51.1.1 d.104.1.1 PDB: 1ady_A* 1adj_A "
         "4rdx_A*"), [
            Tabular.ResiduePair(6, 17, 0.23, 0.00, 0.4388, None),
            Tabular.ResiduePair(7, 18, 0.55, 0.00, 0.5639, None),
            Tabular.ResiduePair(8, 19, 0.37, 0.00, 0.6517, None),
            Tabular.ResiduePair(9, 20, 0.64, 0.00, 0.7241, None),
            Tabular.ResiduePair(10, 21, 0.43, 0.00, 0.7749, None),
            Tabular.ResiduePair(11, 22, 0.37, 0.00, 0.8221, None)
        ]
    )
    _EXPECTED_HITS = (_EXPECTED_HIT_0, _EXPECTED_HIT_1)

    def _test_pairs_equal(self, pair1, pair2):
        self.assertEqual(pair1.i, pair2.i, "Query index identical")
        self.assertEqual(pair1.j, pair2.j, "Template index identical")
        self.assertAlmostEqual(pair1.score, pair2.score, 2, "Scores equal")
        self.assertAlmostEqual(pair1.SS, pair2.SS, 2, "SS equal")
        self.assertAlmostEqual(pair1.probab, pair2.probab, 4, "probab equal")
        self.assertEqual(pair1.dssp, pair2.dssp, "DSSP state equal")

    def _test_hits_equal(self, hit1, hit2):
        self.assertEqual(hit1.name, hit2.name)
        for i, pairs in enumerate(zip(hit1.alignment, hit2.alignment)):
            with self.subTest("Testing residue pair", i=i):
                self._test_pairs_equal(pairs[0], pairs[1])

    def test_parser(self):
        """Parse a dummy tabular file from hhsearch"""
        atab = Tabular(str(DATA_DIR / "test.atab"))
        self.assertEqual(len(atab.hits), 2, 'Read 2 hits')

        for i, hits in enumerate(zip(atab.hits, self._EXPECTED_HITS)):
            with self.subTest("Checking hit against known data", hit=i):
                self._test_hits_equal(hits[0], hits[1])

if __name__ == '__main__':
    unittest.main()
