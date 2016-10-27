import unittest
from tempfile import NamedTemporaryFile

from phyre_engine.tools.hhsuite.parser import Report, Tabular

class TestReportParser(unittest.TestCase):
    """Test a report parser."""

    def test_report(self):
        """Parse a sample report"""
        report_file = NamedTemporaryFile("w")
        report_file.write(report_text)
        report_file.flush()

        report = Report(report_file.name)
        self.assertEqual(report.summary["query"], "d12asa_")
        self.assertEqual(report.summary["match_cols"], 327)

        self.assertEqual(report.summary["neff"], 4.17492)
        self.assertEqual(report.summary["num_searched"], 428)
        self.assertEqual(report.summary["date"], "Tue Oct 25 11:04:48 2016")
        self.assertEqual(report.summary["command"], (
            "hhblits -d "
            "/data/databases/uniprot20_2016_02/uniprot20_2016_02 "
            "-oa3m query.a3m -o report.hhr "
            "-cpu 20 -i /tmp/tmpqa0qxl3b.fasta"))
        self.assertDictEqual(report.summary["num_seqs"],
                {"filtered": 177, "in": 221})

        self.assertEqual(len(report.hits), 273)

        #Pick a couple of arbitrary hits

        # 17 tr|S9CF61|S9CF61_STRAG Unchara 100.0 1.1E-44 5.4E-50  307.6   0.0  158  170-327     1-158 (158)
        #>tr|S9CF61|S9CF61_STRAG Uncharacterized protein OS=Streptococcus agalactiae FSL S3-105 GN=SAG0023_02635 PE=3 SV=1
        #Probab=100.00  E-value=1.1e-44  Score=307.55  Aligned_cols=158  Identities=80%  Similarity=1.253  Sum_probs=155.3  Template_Neff=1.935
        hit_a = report.hits[16].info
        self.assertEqual(hit_a["rank"], 17)
        self.assertEqual(hit_a["prob"], 100.0)
        self.assertEqual(hit_a["evalue"], 1.1e-44)
        self.assertEqual(hit_a["pvalue"], 5.4e-50)
        self.assertEqual(hit_a["score"], 307.6)
        self.assertEqual(hit_a["ss"], 0.0)
        self.assertEqual(hit_a["cols"], 158)
        self.assertEqual(hit_a["query_range"], range(170, 327))
        self.assertEqual(hit_a["template_range"], range(1, 158))
        self.assertEqual(hit_a["template_matches"], 158)
        self.assertEqual(hit_a["name"], (
            "tr|S9CF61|S9CF61_STRAG Uncharacterized protein "
            "OS=Streptococcus agalactiae FSL S3-105 GN=SAG0023_02635 "
            "PE=3 SV=1"))
        self.assertEqual(hit_a["identities"], 80)
        self.assertEqual(hit_a["similarity"], 1.253)
        self.assertEqual(hit_a["sum_probs"], 155.3)
        self.assertEqual(hit_a["template_neff"], 1.935)


        #148 tr|W0U5P1|W0U5P1_9FIRM Asparty  74.7       5 2.6E-05   48.2   0.0  261   23-318   675-951 (957)
        #>tr|W0U5P1|W0U5P1_9FIRM Aspartyl/glutamyl-tRNA(Asn/Gln) amidotransferase,B subunit OS=Ruminococcus bicirculans GN=gatB PE=3 SV=1
        #Probab=74.74  E-value=5  Score=48.23  Aligned_cols=261  Identities=18%  Similarity=0.221  Sum_probs=146.3  Template_Neff=2.645
        hit_b = report.hits[147].info
        self.assertEqual(hit_b["rank"], 148)
        self.assertEqual(hit_b["prob"], 74.7)
        self.assertEqual(hit_b["evalue"], 5.0)
        self.assertEqual(hit_b["pvalue"], 2.6e-5)
        self.assertEqual(hit_b["score"], 48.2)
        self.assertEqual(hit_b["ss"], 0.0)
        self.assertEqual(hit_b["cols"], 261)
        self.assertEqual(hit_b["query_range"], range(23, 318))
        self.assertEqual(hit_b["template_range"], range(675, 951))
        self.assertEqual(hit_b["template_matches"], 957)
        self.assertEqual(hit_b["name"], (
            "tr|W0U5P1|W0U5P1_9FIRM Aspartyl/glutamyl-tRNA(Asn/Gln) "
            "amidotransferase,B subunit OS=Ruminococcus bicirculans "
            "GN=gatB PE=3 SV=1"))
        self.assertEqual(hit_b["identities"], 18)
        self.assertEqual(hit_b["similarity"], 0.221)
        self.assertEqual(hit_b["sum_probs"], 146.3)
        self.assertEqual(hit_b["template_neff"], 2.645)

class TestTabularParser(unittest.TestCase):
    """Test a tabular parser."""

    def test_parser(self):
        """Parse a dummy tabular file from hhsearch"""
        atab_file = NamedTemporaryFile("w")
        atab_file.write(atab_text)
        atab_file.flush()

        atab = Tabular(atab_file.name)
        self.assertEqual(len(atab.hits), 2, 'Read 2 hits')

        hit_a = atab.hits[0]
        self.assertEqual(hit_a.info["name"], (
                "3m4p_A Ehasnrs, asparaginyl-tRNA synthetase, putative; "
                "aminoacyl-tRNA synthetase, tRNA ligase, AARS, translation, "
                "ATP-binding, nucleotide-binding; HET: 4AD; 2.83A "
                "{Entamoeba histolytica} PDB: 3m4q_A"
            ), "Read name from atab file for hit 1")
        self.assertListEqual(hit_a.aln, [
            {"i":  9, "j":  159, "score":  0.10, "SS": 0.00, "probab": 0.3565, "dssp": "H"},
            {"i": 10, "j":  160, "score":  0.35, "SS": 0.00, "probab": 0.3835, "dssp": "H"},
            {"i": 11, "j":  161, "score":  0.15, "SS": 0.00, "probab": 0.4059, "dssp": "H"},
            {"i": 12, "j":  162, "score":  1.39, "SS": 0.00, "probab": 0.4281, "dssp": "H"},
            {"i": 13, "j":  163, "score": -0.19, "SS": 0.00, "probab": 0.4436, "dssp": "H"},
        ], "Alignment for hit A")

report_text = """\
Query         d12asa_
Match_columns 327
No_of_seqs    177 out of 221
Neff          4.17492
Searched_HMMs 428
Date          Tue Oct 25 11:04:48 2016
Command       hhblits -d /data/databases/uniprot20_2016_02/uniprot20_2016_02 -oa3m query.a3m -o report.hhr -cpu 20 -i /tmp/tmpqa0qxl3b.fasta 

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 sp|Q63D03|ASNA_BACCZ Aspartate 100.0  5E-146  3E-151 1043.2   0.0  311    3-324     9-320 (327)
  2 tr|A0A090P0W1|A0A090P0W1_9VIBR 100.0  1E-115  5E-121  793.5   0.0  222   93-324     1-223 (230)
  3 tr|S9W7Y1|S9W7Y1_9TRYP Asparta 100.0  1E-101  7E-107  703.3   0.0  240    1-241     6-253 (253)
  4 tr|R5AR36|R5AR36_9FIRM Asparta 100.0 8.6E-93 3.8E-98  643.5   0.0  261    5-275    18-281 (355)
  5 tr|A0A0R2EZP5|A0A0R2EZP5_LACFE 100.0 3.7E-80 1.6E-85  550.5   0.0  229   26-264     3-233 (293)
  6 tr|W1BHY9|W1BHY9_KLEPN Asparta 100.0 5.9E-80 3.1E-85  544.2   0.0  194    1-194     4-197 (197)
  7 tr|T0E9S1|T0E9S1_CLOSO Asparta 100.0 4.5E-66 2.3E-71  427.3   0.0  118  197-323     2-119 (127)
  8 tr|B7XN60|B7XN60_ENTBH Asparta 100.0 4.9E-66 2.6E-71  446.2   0.0  145    6-153    11-159 (170)
  9 tr|B7XPB1|B7XPB1_ENTBH Asparta 100.0 1.2E-59 5.5E-65  413.5   0.0  237   40-288     3-246 (255)
 10 tr|B7XQP4|B7XQP4_ENTBH Asparta 100.0 4.4E-56 1.9E-61  380.7   0.0  190   40-241     3-199 (203)
 11 tr|A0A074ZZK9|A0A074ZZK9_9TREM 100.0 8.2E-55 3.7E-60  415.1   0.0  304    7-322   270-574 (583)
 12 tr|K1U8F8|K1U8F8_9ZZZZ Asparta 100.0 7.9E-53 3.5E-58  354.8   0.0  141   38-180     3-145 (178)
 13 tr|K1TDC4|K1TDC4_9ZZZZ Asparta 100.0 3.4E-52 1.5E-57  331.6   0.0  118  174-300     2-119 (119)
 14 tr|A0A0P8SA41|A0A0P8SA41_ENTCL 100.0 5.3E-49 2.4E-54  317.6   0.0  124  105-228     1-124 (124)
 15 tr|J9G8M4|J9G8M4_9ZZZZ Asparta 100.0   4E-48 1.8E-53  311.1   0.0  112    8-120     9-120 (122)
 16 tr|X4JHK5|X4JHK5_SALEN Asparag 100.0 8.3E-47 3.8E-52  311.1   0.0  137    1-137     4-140 (140)
 17 tr|S9CF61|S9CF61_STRAG Unchara 100.0 1.1E-44 5.4E-50  307.6   0.0  158  170-327     1-158 (158)
 18 tr|A0A0E9FL15|A0A0E9FL15_CHLTH 100.0 1.9E-43 1.1E-48  276.7   0.0   78  249-326     1-78  (83)
 19 tr|Q3EXH4|Q3EXH4_BACTI Asparta 100.0 1.2E-40 5.3E-46  274.6   0.0  106    6-113    12-117 (147)
 20 tr|B7XP24|B7XP24_ENTBH Asparta 100.0   5E-40 2.2E-45  273.7   0.0  149   40-191     3-158 (158)
 21 tr|G4ELX7|G4ELX7_MYCIO Asparag 100.0 4.5E-33 2.2E-38  221.3   0.0   92    1-93      7-98  (98)
 22 tr|S7XF60|S7XF60_SPRLO Putativ  99.9 4.6E-31 2.1E-36  222.5   0.0  107  152-258    29-161 (167)
 23 tr|B7XQM8|B7XQM8_ENTBH Asparta  99.9 1.2E-28 5.4E-34  200.4   0.0  104   40-146     3-110 (128)
 24 tr|K6E9S2|K6E9S2_9BACI Asparag  99.9 1.5E-28 6.8E-34  177.3   0.0   55   64-119     1-55  (62)
 25 tr|W1XGP1|W1XGP1_9ZZZZ Asparta  99.9 5.3E-28 2.4E-33  178.4   0.0   69   49-117     2-70  (70)
 26 tr|A0A0J9TKM5|A0A0J9TKM5_PLAVI  99.9 1.7E-27   8E-33  188.2   0.0   93  125-217     5-98  (98)
 27 tr|K1TZK1|K1TZK1_9ZZZZ Asparta  99.8   4E-25 1.8E-30  158.0   0.0   57  197-262     2-58  (58)
 28 tr|D8FGH8|D8FGH8_9FIRM Putativ  99.8 4.3E-25   2E-30  165.9   0.0   69    8-78      3-71  (73)
 29 tr|K1Y4G7|K1Y4G7_9BACT Asparta  99.8 5.3E-25 2.5E-30  167.8   0.0   73    3-75      6-78  (78)
 30 tr|S7W4K5|S7W4K5_SPRLO Unchara  99.8 1.4E-24 6.3E-30  177.5   0.0   97    6-105    26-122 (130)
 31 tr|A0A0H2REP6|A0A0H2REP6_MORMO  99.8 1.6E-23 7.1E-29  153.9   0.0   66  261-326     1-66  (67)
 32 tr|A4GX00|A4GX00_9LACT AsnA (F  99.8 9.4E-23 4.2E-28  141.5   0.0   36  283-318     2-37  (49)
 33 tr|A0A090P1X0|A0A090P1X0_9VIBR  99.5 9.2E-18 4.3E-23  126.2   0.0   55    2-56     15-69  (69)
 34 tr|B7XNE7|B7XNE7_ENTBH Asparta  99.3   2E-15 8.8E-21  108.6   0.0   54   39-95      2-55  (56)
 35 tr|K1TB71|K1TB71_9ZZZZ Asparta  99.1 2.4E-13 1.1E-18  118.1   0.0   59    2-60     94-152 (176)
 36 tr|W1VJV4|W1VJV4_STRPA Asparta  98.6 2.1E-10 9.8E-16   84.5   0.0   54    1-54      4-57  (57)
 37 tr|A0A090PG64|A0A090PG64_9VIBR  98.5 6.1E-10 2.7E-15   74.9   0.0   28   60-87      5-32  (37)
 38 tr|K0NLJ8|K0NLJ8_9LACO Unchara  98.5 8.3E-10 3.7E-15   72.6   0.0   28  299-326     1-28  (33)
 39 tr|A0A084U2P1|A0A084U2P1_MYCIO  98.2   8E-09 3.8E-14   99.8   0.0  291    1-312    11-308 (325)
 40 tr|B7XQ79|B7XQ79_ENTBH Asparta  98.2 8.6E-09 3.9E-14   76.3   0.0   48  106-153     1-48  (59)
 41 tr|A0A0G1B9N3|A0A0G1B9N3_9BACT  97.7 2.7E-07 1.8E-12   94.1   0.0  255   23-312    53-328 (344)
 42 tr|X1IS21|X1IS21_9ZZZZ Unchara  97.5 1.5E-06 6.8E-12   59.1   0.0   28  298-325     1-28  (35)
 43 tr|R5PEQ4|R5PEQ4_9BACT Asparta  97.2   1E-05 4.6E-11   60.2   0.0   36    2-37     16-51  (56)
 44 tr|B7XNP0|B7XNP0_ENTBH Unchara  96.9 3.8E-05 1.7E-10   59.5   0.0   47   39-97      2-48  (68)
 45 tr|K2C8Q9|K2C8Q9_9BACT Unchara  96.7   9E-05   4E-10   47.9   0.0   21  305-325     1-21  (27)
 46 tr|R1C3V4|R1C3V4_EMIHU Lysyl-t  95.3   0.005 2.2E-08   46.2   0.0   29  283-311     5-33  (52)
 47 tr|A0A079G0B3|A0A079G0B3_ECOLX  95.2  0.0056 2.9E-08   47.0   0.0   38  273-311     5-42  (51)
 48 tr|U1F973|U1F973_9ACTN Lysyl-t  95.1  0.0056   3E-08   49.7   0.0   40  270-311    19-58  (68)
 49 tr|F3VDV0|F3VDV0_SHIDY Lysyl-t  95.1  0.0062   3E-08   43.5   0.0   28  284-311     1-28  (37)
 50 tr|T0YGT2|T0YGT2_9ZZZZ Lysyl-t  95.1  0.0055 3.3E-08   57.9   0.0   55  255-311   107-163 (165)
 51 tr|D9PH19|D9PH19_9ZZZZ Lysyl-t  94.9  0.0078 4.6E-08   63.6   0.0  274    7-312    52-360 (375)
 52 tr|W7QHE3|W7QHE3_9GAMM Unchara  94.4   0.017 9.2E-08   53.0   0.0   56  255-311    67-126 (133)
 53 tr|X1CIB2|X1CIB2_9ZZZZ Unchara  94.3   0.022   1E-07   51.0   0.0   56  255-311     8-68  (144)
 54 tr|A0A0G9GN21|A0A0G9GN21_LACPN  94.3   0.018 1.1E-07   50.7   0.0   56  256-312    34-93  (103)
 55 tr|G4NE53|G4NE53_MAGO7 Lysine-  94.0   0.025 1.4E-07   63.4   0.0   42  270-312   571-612 (638)
 56 tr|Q48986|Q48986_MYCCA Lys-tRN  93.6   0.046   2E-07   42.3   0.0   32  280-311    18-49  (58)
 57 tr|S2EWF1|S2EWF1_9ARCH Unchara  93.6   0.048 2.2E-07   45.2   0.0   68  243-311     4-72  (83)
 58 tr|A0A0C2D2D0|A0A0C2D2D0_9BILA  93.6   0.047 2.2E-07   43.5   0.0   36  275-311     6-41  (65)
 59 tr|A0A0G1WAG4|A0A0G1WAG4_9BACT  93.3   0.042 2.7E-07   60.1   0.0  257   17-312   152-432 (463)
 60 tr|K1T7I3|K1T7I3_9ZZZZ Lysyl-t  93.3   0.056 2.7E-07   53.8   0.0   56  255-311   118-177 (260)
 61 tr|A0A0G0L971|A0A0G0L971_9BACT  93.3   0.062 2.8E-07   42.0   0.0   28    3-30     33-60  (61)
 62 tr|V2PR40|V2PR40_SALET Asparag  93.3   0.063 2.8E-07   37.0   0.0   29    1-29      4-32  (32)
 63 tr|A0A0G0WIP2|A0A0G0WIP2_9BACT  92.8   0.076 3.9E-07   51.2   0.0  127  173-320    42-175 (185)
 64 tr|A0A0M1J520|A0A0M1J520_9GAMM  92.8    0.09   4E-07   41.1   0.0   35  277-312    23-57  (60)
 65 tr|A0A090PT73|A0A090PT73_9FLAO  92.7   0.095 4.2E-07   45.8   0.0   38  273-311     5-42  (116)
 66 tr|A0A0G7ZNA4|A0A0G7ZNA4_9MOLU  92.6   0.098 4.5E-07   53.0   0.0  290    3-321    15-314 (325)
 67 tr|U1QP96|U1QP96_9EURY Asparty  92.5   0.081 4.8E-07   47.5   0.0   68  243-311    32-100 (111)
 68 tr|A0A0K0F1X5|A0A0K0F1X5_9BILA  92.5     0.1 4.9E-07   46.7   0.0   41  271-312    50-90  (115)
 69 tr|W9LCI1|W9LCI1_FUSOX Asparty  92.4   0.086   5E-07   51.3   0.0   67  244-311   100-167 (185)
 70 tr|W4YTD6|W4YTD6_STRPU Unchara  92.4    0.12 5.2E-07   42.6   0.0   38  273-311    14-51  (80)
 71 tr|X0UCQ7|X0UCQ7_9ZZZZ Unchara  92.3     0.1 5.7E-07   42.1   0.0   39  273-311    12-50  (61)
 72 sp|O67258|SYK_AQUAE Lysine--tR  92.2    0.12 5.8E-07   57.3   0.0  263    7-311   275-582 (597)
 73 tr|A0A0E0KFX6|A0A0E0KFX6_ORYPU  92.2     0.1 5.9E-07   59.1   0.0   56  255-311   556-615 (652)
 74 tr|B7BCA8|B7BCA8_9PORP tRNA li  92.1    0.14 6.3E-07   48.7   0.0   56  255-311    66-125 (208)
 75 tr|B7XRA1|B7XRA1_ENTBH Asparta  92.0    0.15 6.8E-07   43.1   0.0   80   40-122     3-86  (93)
 76 tr|A0A0G1IIJ4|A0A0G1IIJ4_9BACT  91.9    0.12   7E-07   55.5   0.0   80  227-312   291-372 (382)
 77 tr|A0A0G1D8P3|A0A0G1D8P3_9BACT  91.5    0.13 8.5E-07   52.4   0.0  203   86-311    22-245 (256)
 78 tr|M8AA95|M8AA95_TRIUA Unchara  91.3    0.18 9.4E-07   53.3   0.0   41  271-312   320-360 (368)
 79 tr|A0A078CAH0|A0A078CAH0_BRANA  91.3    0.18 9.7E-07   57.8   0.0   58  254-312   536-597 (720)
 80 tr|A0A0G1M6Q9|A0A0G1M6Q9_9BACT  91.2     0.2   1E-06   51.7   0.0  228   61-313    16-268 (299)
 81 tr|Q9FEL1|Q9FEL1_TOBAC Lysyl-t  91.0    0.24 1.1E-06   46.6   0.0   58  253-311    55-116 (169)
 82 tr|G7DZZ8|G7DZZ8_MIXOS Unchara  90.9    0.25 1.1E-06   53.6   0.0  112  191-308   515-656 (691)
 83 tr|G5DY35|G5DY35_9PIPI Putativ  90.9    0.21 1.2E-06   52.1   0.0  206   90-311    42-293 (302)
 84 tr|A0A079FAQ6|A0A079FAQ6_ECOLX  90.5    0.31 1.4E-06   39.8   0.0   39  272-311    27-65  (74)
 85 tr|R5CXA5|R5CXA5_9FIRM Asparta  90.5    0.31 1.4E-06   38.7   0.0   28    1-28     14-41  (63)
 86 tr|A0A0G0IMH4|A0A0G0IMH4_9BACT  90.4    0.26 1.4E-06   53.2   0.0  276    6-310    88-398 (404)
 87 tr|A0A095BRF7|A0A095BRF7_SCHHA  90.2    0.34 1.5E-06   46.7   0.0   58  253-311    39-100 (206)
 88 tr|X1CXE8|X1CXE8_9ZZZZ Unchara  90.1    0.35 1.6E-06   35.1   0.0   31  282-312     3-33  (39)
 89 tr|A0A0A0L4F7|A0A0A0L4F7_CUCSA  89.5    0.44 1.9E-06   39.9   0.0   40  271-311    43-82  (84)
 90 tr|A0A059WSR4|A0A059WSR4_9BACT  88.9    0.47 2.4E-06   46.0   0.0   46  275-321   131-176 (176)
 91 tr|A0A0J9VCE3|A0A0J9VCE3_PLAVI  88.9    0.49 2.4E-06   53.5   0.0   45  266-311   625-669 (677)
 92 tr|A0A0G2UQE3|A0A0G2UQE3_THAPS  88.8    0.44 2.5E-06   54.7   0.0   59  253-312   599-661 (669)
 93 tr|G1XRE9|G1XRE9_ARTOA Lysine-  88.7    0.44 2.6E-06   54.5   0.0   42  270-312   546-587 (626)
 94 tr|A0A0G4LDB2|A0A0G4LDB2_9PEZI  88.6    0.56 2.6E-06   41.6   0.0   40  273-313    39-78  (110)
 95 tr|Q05FQ1|Q05FQ1_CARRP Putativ  88.2    0.62   3E-06   40.0   0.0   25  286-310    61-85  (87)
 96 sp|P43151|GPA_LEIDO Putative g  87.9    0.73 3.2E-06   48.6   0.0   58  253-311   258-319 (464)
 97 tr|B9Y6V6|B9Y6V6_9FIRM Asparta  87.9    0.52 3.3E-06   53.6   0.0   76  235-312   468-548 (580)
 98 tr|A6DS64|A6DS64_9BACT Lysyl-t  87.6    0.53 3.5E-06   49.2   0.0  266    5-316    10-291 (292)
 99 tr|A6MK47|A6MK47_CALJA Lysyl-t  87.3    0.85 3.8E-06   43.1   0.0   40  271-311   112-151 (180)
100 tr|C5LM59|C5LM59_PERM5 Lysyl-t  87.2    0.88 3.9E-06   47.2   0.0   51  260-311   316-368 (388)
101 tr|F2TYB2|F2TYB2_SALR5 Lysine-  87.2    0.86 3.9E-06   51.1   0.0   59  252-311   489-551 (796)
102 tr|A0A0M7R9W6|A0A0M7R9W6_ECOLX  87.1    0.89 3.9E-06   34.9   0.0   29  275-304     7-35  (50)
103 tr|A0A0C3FBM0|A0A0C3FBM0_9HOMO  87.0    0.91   4E-06   50.9   0.0   71  244-318   739-810 (895)
104 tr|J8USK7|J8USK7_PSEPU Unchara  86.8    0.96 4.3E-06   38.4   0.0   18  290-307     7-24  (88)
105 tr|W7K2V9|W7K2V9_PLAFA Asparag  86.5       1 4.6E-06   51.1   0.0   58  254-312   991-1048(1058)
106 tr|A0A0N5CP41|A0A0N5CP41_THECL  86.4     1.1 4.7E-06   34.2   0.0   36  271-307    12-47  (48)
107 tr|A0A0A9R9K1|A0A0A9R9K1_ARUDO  86.3       1 4.8E-06   37.2   0.0   39  273-312    25-63  (71)
108 tr|A0A098ED06|A0A098ED06_9ZZZZ  85.9     1.2 5.2E-06   35.7   0.0   32  229-260    24-57  (62)
109 tr|W7TMK7|W7TMK7_9STRA Lysyl-t  85.8     1.2 5.4E-06   43.1   0.0   38  273-311    23-60  (190)
110 tr|T0YBN7|T0YBN7_9ZZZZ tRNA sy  85.5     1.3 5.7E-06   41.9   0.0  134  156-312    18-152 (175)
111 tr|T2JYP4|T2JYP4_CROWT Lysyl-t  85.4     1.3 5.9E-06   39.4   0.0   37  275-312    10-46  (116)
112 tr|H9JJU4|H9JJU4_BOMMO Unchara  85.3     1.3 5.9E-06   36.9   0.0   30  286-317    42-71  (78)
113 tr|S9QX59|S9QX59_9RHOB Asparag  84.8     1.5 6.6E-06   37.7   0.0   41  278-318    46-88  (92)
114 tr|A0A0D2X2W4|A0A0D2X2W4_CAPO3  84.6     1.5 6.8E-06   48.0   0.0   55  256-311   564-622 (652)
115 tr|A0A067QMI8|A0A067QMI8_9HOMO  84.4     1.6 7.1E-06   35.6   0.0   42   92-133    12-53  (68)
116 tr|M1CK15|M1CK15_SOLTU Unchara  84.1     1.7 7.5E-06   35.1   0.0   38  271-309    16-53  (64)
117 tr|A0A0G1BLE2|A0A0G1BLE2_9BACT  84.1     1.7 7.5E-06   34.0   0.0   43  268-311     7-49  (54)
118 tr|A0A0E8TK72|A0A0E8TK72_STREE  83.1       2 8.9E-06   35.3   0.0   29  283-311     4-35  (70)
119 tr|R5QIY9|R5QIY9_9FIRM Lysine-  82.9       2 9.2E-06   45.9   0.0   55  256-311   183-241 (396)
120 tr|A0A0R3TSL9|A0A0R3TSL9_HYMNN  82.8     2.1 9.3E-06   40.6   0.0   42  285-326    49-91  (176)
121 tr|D4JFP2|D4JFP2_9FIRM Lysine-  82.6     1.8 9.6E-06   50.0   0.0  201   93-311   247-485 (649)
122 tr|A0A017T1J2|A0A017T1J2_9DELT  82.3       2   1E-05   47.1   0.0  278    5-317   116-430 (431)
123 tr|H2KSF2|H2KSF2_CLOSI Lysine-  82.2     2.2   1E-05   47.4   0.0   58  253-311   462-523 (647)
124 tr|D8LM77|D8LM77_ECTSI Lysine-  81.8     2.5 1.1E-05   47.6   0.0  139  170-311   390-551 (849)
125 tr|M7SLK5|M7SLK5_EUTLA Putativ  81.6     2.1 1.1E-05   51.2   0.0  254   24-317   271-558 (914)
126 tr|A0A0F6TGI3|A0A0F6TGI3_9BETA  81.0     2.7 1.2E-05   38.8   0.0   37  204-240    95-132 (144)
127 tr|W4F9U4|W4F9U4_9STRA Unchara  80.7     2.9 1.3E-05   40.3   0.0   58  253-311    45-106 (194)
128 tr|A0A0N8EBA7|A0A0N8EBA7_9CRUS  80.5       3 1.3E-05   42.4   0.0   70  245-317   215-285 (292)
129 tr|G5ANI0|G5ANI0_HETGA Lysine-  80.3     2.8 1.4E-05   46.0   0.0   40  271-311   379-418 (447)
130 tr|A0A0G4LE24|A0A0G4LE24_9PEZI  79.8     3.3 1.5E-05   42.8   0.0   42  273-315   266-307 (337)
131 tr|V4KXH3|V4KXH3_EUTSA Unchara  79.6     2.8 1.5E-05   48.3   0.0  258   24-317   334-621 (628)
132 tr|A0A0F7UAI7|A0A0F7UAI7_NEOCL  79.5     3.3 1.5E-05   48.5   0.0   41  270-311  1016-1056(1065)
133 tr|A0A089QLH2|A0A089QLH2_9PROC  79.4     3.4 1.5E-05   33.4   0.0   22  166-187     6-27  (63)
134 tr|A0A060SS49|A0A060SS49_PYCCI  79.0     3.6 1.6E-05   48.3   0.0   40  271-311   552-591 (1431)
135 tr|A0A0G4MVL4|A0A0G4MVL4_9PEZI  78.9     3.4 1.6E-05   36.0   0.0   41  271-312    42-82  (88)
136 tr|D1NED1|D1NED1_HAEIF Lysine-  78.2       4 1.8E-05   39.1   0.0   40  271-311   136-175 (183)
137 tr|A0A078F524|A0A078F524_BRANA  78.1       4 1.8E-05   35.1   0.0   38  273-311    27-64  (89)
138 tr|A0A086M299|A0A086M299_TOXGO  77.5     4.3 1.9E-05   46.0   0.0   40  271-311   770-809 (818)
139 tr|E9ICF4|E9ICF4_SOLIN Putativ  77.3     4.4   2E-05   45.0   0.0   66  245-311   443-509 (686)
140 tr|W7W9Z0|W7W9Z0_9BURK Lysine-  76.7     4.7 2.1E-05   46.7   0.0   40  271-311   460-499 (1129)
141 tr|W4FAF2|W4FAF2_9STRA Unchara  76.7     4.8 2.1E-05   40.5   0.0   58  253-311   112-173 (261)
142 tr|L1LCI1|L1LCI1_BABEQ Lysine-  76.4     4.9 2.2E-05   46.9   0.0   38  273-311  1206-1243(1262)
143 tr|A0A0D7BRC0|A0A0D7BRC0_9HOMO  76.4     4.9 2.2E-05   33.1   0.0   24  282-305     8-31  (69)
144 tr|A0A067BUS2|A0A067BUS2_SAPPC  76.3       5 2.2E-05   36.0   0.0   44   54-100    21-66  (113)
145 tr|A0A0F8YFG0|A0A0F8YFG0_9ZZZZ  76.1     5.1 2.3E-05   42.4   0.0   49  261-310    53-103 (403)
146 tr|A0A0G2YPK2|A0A0G2YPK2_BACIU  75.7     4.2 2.4E-05   44.1   0.0  260   22-312    34-330 (347)
147 tr|A0A0D2IX76|A0A0D2IX76_9CHLO  75.2     5.6 2.5E-05   34.9   0.0   39  273-312    56-94  (99)
148 tr|W0U5P1|W0U5P1_9FIRM Asparty  74.7       5 2.6E-05   48.2   0.0  261   23-318   675-951 (957)
149 tr|W4KAS8|W4KAS8_9HOMO Unchara  74.5       6 2.7E-05   31.6   0.0   41   75-116    17-57  (57)
150 tr|A0A0A2MJI3|A0A0A2MJI3_9FLAO  73.6     6.6 2.9E-05   31.9   0.0   31    9-39     15-45  (62)
151 tr|A0A037UZF0|A0A037UZF0_9RHIZ  73.4     6.7   3E-05   37.4   0.0   37   61-97     46-87  (169)
152 tr|H9IUD4|H9IUD4_BOMMO Unchara  72.7     7.1 3.2E-05   45.4   0.0  109   95-207   909-1037(1109)
153 tr|D0MYK2|D0MYK2_PHYIT Asparty  72.4     7.3 3.3E-05   39.7   0.0   55  254-309   211-265 (284)
154 tr|V9W5Z7|V9W5Z7_9BACL Putativ  72.4     7.4 3.3E-05   34.3   0.0   19  286-304    48-66  (99)
155 tr|W4XMC8|W4XMC8_STRPU Unchara  71.4     8.1 3.6E-05   38.0   0.0   35  190-225    61-95  (210)
156 tr|G6CVH5|G6CVH5_DANPL Asparty  71.2     8.2 3.7E-05   40.0   0.0   65  244-309   168-233 (325)
157 tr|D7G3Y0|D7G3Y0_ECTSI Lysine-  71.0     8.4 3.7E-05   43.3   0.0   41  271-312   648-688 (710)
158 tr|B9SM60|B9SM60_RICCO Putativ  70.4     8.8 3.9E-05   31.7   0.0   34  274-307     2-35  (67)
159 tr|A0A0F7U2K7|A0A0F7U2K7_9EURO  70.3     8.9   4E-05   37.8   0.0   25  284-308    46-70  (212)
160 tr|A0A0G5IVB7|A0A0G5IVB7_PSEAI  70.2     8.6   4E-05   29.1   0.0   29  293-321    11-39  (40)
161 tr|B9TDH8|B9TDH8_RICCO Putativ  69.4     9.6 4.3E-05   33.8   0.0   18  281-299    78-95  (102)
162 tr|A0A060T8I3|A0A060T8I3_BLAAD  69.4     9.6 4.3E-05   38.6   0.0   50  246-297    31-87  (262)
163 tr|F0VCT3|F0VCT3_NEOCL Lysyl-t  69.0     9.9 4.4E-05   46.9   0.0   61  250-311  2478-2538(2547)
164 tr|R8NW90|R8NW90_BACCE Unchara  68.5      10 4.6E-05   30.0   0.0   34  274-307     8-41  (53)
165 tr|X1K747|X1K747_9ZZZZ Unchara  67.6     8.9   5E-05   36.6   0.0   75  236-312    20-99  (133)
166 tr|A0A093XZD4|A0A093XZD4_9PEZI  67.4     9.5   5E-05   44.6   0.0   55  256-311   511-567 (645)
167 tr|A0A093V0G4|A0A093V0G4_TALMA  67.3      11 5.1E-05   42.1   0.0  213   54-305   292-527 (659)
168 tr|K3YKL8|K3YKL8_SETIT Unchara  66.2      12 5.6E-05   31.3   0.0   16  294-309     6-21  (72)
169 tr|F6I4L4|F6I4L4_VITVI Lysine-  65.9      13 5.7E-05   44.7   0.0   58  253-311  1492-1553(1578)
170 tr|C5LAS7|C5LAS7_PERM5 Asparty  64.8      14 6.2E-05   37.5   0.0   80  229-317   170-250 (257)
171 tr|A3TSB1|A3TSB1_OCEBH Putativ  63.5      15 6.8E-05   36.4   0.0   71  232-306   102-176 (213)
172 tr|A0A068VDE3|A0A068VDE3_COFCA  63.4      14 6.9E-05   43.1   0.0   41  271-312   616-656 (664)
173 tr|A0A0C3J1A7|A0A0C3J1A7_9RHOO  63.2      16   7E-05   29.6   0.0   32   21-52     23-54  (58)
174 tr|A0A0F9MLR7|A0A0F9MLR7_9ZZZZ  62.3      17 7.4E-05   27.7   0.0   28  120-149     3-30  (42)
175 tr|U2R9P7|U2R9P7_9FIRM Unchara  61.9      17 7.6E-05   30.8   0.0   25  277-301     7-31  (71)
176 tr|A0A0F8B681|A0A0F8B681_CERFI  61.2      17   8E-05   43.1   0.0  247   24-304   443-709 (1054)
177 tr|F3GQQ4|F3GQQ4_PSESJ Nicotin  61.0      18 8.1E-05   28.1   0.0   22   23-44      9-30  (47)
178 tr|A0A0L9TWR3|A0A0L9TWR3_PHAAN  61.0      18 8.1E-05   35.9   0.0   32  271-302   139-170 (211)
179 tr|I3KV14|I3KV14_ORENI Unchara  59.3      20 9.1E-05   35.1   0.0   43  100-142    23-65  (166)
180 tr|A0A0P6D0T2|A0A0P6D0T2_9CRUS  59.3      20 9.1E-05   39.8   0.0   40  271-311   507-546 (568)
181 tr|C0NGY7|C0NGY7_AJECG Asparty  59.1      20 9.2E-05   42.8   0.0  220   54-308   734-995 (1009)
182 tr|W7TE52|W7TE52_9STRA Asparty  57.5      20  0.0001   42.0   0.0  256   25-318   312-595 (601)
183 tr|R5GV28|R5GV28_9BACT Unchara  54.7      27 0.00012   29.6   0.0   27  297-323    41-67  (73)
184 tr|L1JTG7|L1JTG7_GUITH Unchara  52.7      29 0.00014   38.0   0.0   43  273-318    90-132 (383)
185 tr|E4TU80|E4TU80_MARTH Unchara  52.5      32 0.00014   32.7   0.0   44   33-77     63-111 (145)
186 tr|A0A0G1C013|A0A0G1C013_9BACT  52.5      30 0.00014   37.9   0.0  120  187-312    78-202 (394)
187 tr|X1KLZ9|X1KLZ9_9ZZZZ Unchara  52.5      26 0.00014   32.6   0.0   72  110-183     6-80  (106)
188 tr|H1XVZ9|H1XVZ9_9BACT TonB fa  52.2      32 0.00014   37.8   0.0   70    3-74    376-445 (478)
189 tr|A0A0A2KRV9|A0A0A2KRV9_PENIT  52.1      32 0.00014   41.2   0.0   46  265-311   464-509 (1059)
190 tr|E9LR59|E9LR59_9EUKA Lysyl-t  51.9      33 0.00015   37.6   0.0   32  280-311   217-248 (463)
191 tr|J2U6W9|J2U6W9_9BURK Unchara  51.8      32 0.00015   30.2   0.0   39  218-256    41-79  (82)
192 tr|E1ZP53|E1ZP53_CHLVA Putativ  51.3      34 0.00015   38.6   0.0   40  271-311   412-451 (624)
193 tr|K0REP7|K0REP7_THAOC Unchara  51.2      34 0.00015   30.3   0.0   35  179-213    27-63  (94)
194 tr|A0A0F8Z584|A0A0F8Z584_9ZZZZ  51.1      34 0.00015   28.6   0.0   31    1-31     22-52  (66)
195 tr|X1GNU4|X1GNU4_9ZZZZ Unchara  50.6      35 0.00016   26.6   0.0   32  133-164     2-33  (46)
196 tr|A0A067TRT7|A0A067TRT7_9AGAR  50.2      36 0.00016   37.3   0.0   44   97-141     9-53  (459)
197 tr|G0AIS0|G0AIS0_COLFT Unchara  49.9      37 0.00016   34.8   0.0   25   91-115    73-97  (251)
198 tr|A0A078D0F7|A0A078D0F7_BRANA  49.9      37 0.00017   36.2   0.0   30  283-312   195-224 (353)
199 tr|T1C6F5|T1C6F5_9ZZZZ Asparag  49.6      35 0.00017   32.9   0.0   79  229-312    49-127 (138)
200 tr|T0REB3|T0REB3_9DELT WYL dom  48.9      39 0.00018   35.8   0.0   30   84-113    35-64  (336)
201 tr|V6LXH9|V6LXH9_9EUKA Unchara  48.7      40 0.00018   32.6   0.0   35    9-43     47-81  (164)
202 tr|A0A0G2FYK3|A0A0G2FYK3_9PEZI  48.4      41 0.00018   37.2   0.0  105    4-123    16-130 (490)
203 tr|M7XA29|M7XA29_9BACT Unchara  48.2      41 0.00018   29.4   0.0   32  274-305    46-77  (84)
204 tr|I7ATP4|I7ATP4_9CAUD Unchara  48.2      34 0.00018   33.9   0.0   55  190-244     7-64  (155)
205 tr|W1SLJ3|W1SLJ3_9BACI Unchara  48.1      41 0.00018   34.4   0.0   44  273-316   184-228 (246)
206 tr|A0A0C3PMV5|A0A0C3PMV5_PHLGI  48.0      41 0.00018   34.0   0.0   31   70-100   130-160 (227)
207 tr|A0A0D7X2X4|A0A0D7X2X4_9BACL  47.7      42 0.00019   35.0   0.0   23  221-243    50-72  (291)
208 tr|A0A0R2AJ50|A0A0R2AJ50_LACPA  47.4      41 0.00019   28.3   0.0   32  166-197    25-56  (64)
209 tr|D0NTV5|D0NTV5_PHYIT Asparty  46.4      46  0.0002   36.0   0.0   62  254-318   321-382 (388)
210 tr|E9EJY1|E9EJY1_METRA Methyl   46.2      42 0.00021   32.8   0.0   54  217-270    18-71  (145)
211 tr|A0A0Q9GYJ6|A0A0Q9GYJ6_9BACI  45.6      47 0.00021   32.7   0.0   36  271-306     5-40  (169)
212 tr|A0A0G0YPD7|A0A0G0YPD7_9BACT  45.5      48 0.00022   34.1   0.0   55  154-208    66-120 (255)
213 tr|A0A0Q7F9L2|A0A0Q7F9L2_9CAUL  45.0      50 0.00022   28.7   0.0   34  235-268    15-48  (81)
214 tr|B2ACG1|B2ACG1_PODAN Podospo  44.3      52 0.00023   33.3   0.0   19  203-221    97-115 (224)
215 tr|A0A0G0L7W7|A0A0G0L7W7_9BACT  43.6      52 0.00024   32.8   0.0   37  229-265    53-89  (180)
216 tr|A6DL15|A6DL15_9BACT Unchara  41.8      60 0.00027   30.9   0.0   30  169-198    44-73  (140)
217 tr|E4XWZ8|E4XWZ8_OIKDI Unchara  41.6      61 0.00027   37.1   0.0   36  263-298    46-81  (652)
218 tr|A0A067PEU0|A0A067PEU0_9HOMO  41.1      63 0.00028   31.7   0.0   30   39-68    125-154 (172)
219 tr|A0A0N0Y7P1|A0A0N0Y7P1_9ACTN  40.7      64 0.00028   29.2   0.0   48   42-89     10-57  (102)
220 tr|A0A0F9V184|A0A0F9V184_9ZZZZ  40.5      62 0.00029   31.2   0.0   43   62-104    20-65  (137)
221 tr|J9E6C3|J9E6C3_WUCBA Unchara  39.7      68  0.0003   29.9   0.0   48   45-94     50-97  (121)
222 tr|A0A0J6EHQ9|A0A0J6EHQ9_9BACI  39.1      70 0.00031   27.7   0.0   28  196-223    34-61  (78)
223 tr|A0A0P9J0K9|A0A0P9J0K9_9PSED  38.5      73 0.00033   32.7   0.0   28  172-199   200-227 (241)
224 tr|E2NUA3|E2NUA3_9FIRM Unchara  38.3      74 0.00033   26.0   0.0   22   87-108    21-42  (55)
225 tr|L0L8H6|L0L8H6_9CAUD Unchara  37.9      69 0.00034   31.5   0.0   48  133-181     3-51  (148)
226 tr|A0A0F0HS90|A0A0F0HS90_9PSEU  37.9      76 0.00034   28.7   0.0   60  226-289    29-88  (99)
227 tr|S8DGJ7|S8DGJ7_9LAMI Unchara  37.7      77 0.00034   32.6   0.0   52  116-169    60-113 (243)
228 tr|A0A0B0EA93|A0A0B0EA93_9BACT  37.3      78 0.00035   31.6   0.0   53  271-323    18-72  (193)
229 tr|A0A0R0HVZ6|A0A0R0HVZ6_SOYBN  37.3      76 0.00035   33.5   0.0   54  116-169   115-168 (272)
230 tr|G9MSI5|G9MSI5_HYPVG Unchara  37.2      79 0.00035   35.0   0.0   60  119-179    55-121 (456)
231 tr|A0A0F9Y4R0|A0A0F9Y4R0_TRIHA  37.0      80 0.00036   34.4   0.0   62  118-180    54-122 (390)
232 tr|U1X922|U1X922_9BURK Methylt  36.9      81 0.00036   29.1   0.0   43    5-47     31-73  (111)
233 tr|G8M0B5|G8M0B5_CLOCD Unchara  35.5      83 0.00039   28.7   0.0   55  125-179    29-84  (95)
234 tr|A0A0J6GXZ2|A0A0J6GXZ2_9BACI  34.9      75 0.00041   30.0   0.0   29  284-312    45-73  (109)
235 tr|C1EBJ0|C1EBJ0_MICSR Unchara  34.6      93 0.00041   33.8   0.0   56  238-295   256-312 (374)
236 tr|M2NJW0|M2NJW0_BAUCO Lysine-  34.4      94 0.00042   38.3   0.0   39  274-313  1418-1456(1487)
237 tr|A0A0D2IXC4|A0A0D2IXC4_9CHLO  34.4      94 0.00042   29.7   0.0   34  170-204    92-125 (140)
238 tr|X1J2A8|X1J2A8_9ZZZZ Unchara  34.0      96 0.00043   29.1   0.0   33    2-34      4-36  (123)
239 tr|A0A0M2WTF0|A0A0M2WTF0_9BURK  34.0      91 0.00043   35.0   0.0   77  216-304   216-292 (413)
240 tr|E8YJA2|E8YJA2_9BURK Putativ  33.8      97 0.00043   27.5   0.0   40   45-84      7-46  (87)
241 tr|I2G2J9|I2G2J9_USTH4 Unchara  33.6      99 0.00044   26.9   0.0   24   45-68     32-55  (77)
242 tr|T0CYU0|T0CYU0_9BACL Unchara  32.9   1E+02 0.00046   35.1   0.0   62  173-246   428-493 (576)
243 tr|A0A0E0LIU8|A0A0E0LIU8_ORYPU  32.8   1E+02 0.00046   33.1   0.0   66  177-256   252-318 (337)
244 tr|A0A0C3C5V9|A0A0C3C5V9_HEBCY  32.7   1E+02 0.00046   34.3   0.0   53  154-207   169-230 (463)
245 tr|A0A061IL71|A0A061IL71_CRIGR  32.0 1.1E+02 0.00048   30.7   0.0   37  203-241    77-113 (190)
246 tr|A0A0P4VRZ8|A0A0P4VRZ8_9EUCA  31.7 1.1E+02 0.00049   33.4   0.0   29   71-100   146-174 (388)
247 tr|A0A0G0CFV7|A0A0G0CFV7_9BACT  31.2 1.2E+02 0.00051   29.1   0.0   27  183-209    93-119 (136)
248 tr|T1BE53|T1BE53_9ZZZZ Asparty  31.1      86 0.00052   35.5   0.0  117  173-312   215-337 (376)
249 tr|Q14KG0|Q14KG0_SPICI Putativ  30.4 1.2E+02 0.00054   30.0   0.0   29  164-192    81-109 (161)
250 tr|R6KA16|R6KA16_9FIRM Unchara  30.4 1.2E+02 0.00054   27.4   0.0   29  221-249    58-86  (94)
251 tr|M7C583|M7C583_CHEMY Unchara  30.1 1.2E+02 0.00055   27.8   0.0   31  216-246    50-80  (104)
252 tr|V6LT09|V6LT09_9EUKA Unchara  29.6 1.3E+02 0.00057   32.5   0.0   35    9-43     47-81  (331)
253 tr|A8SD86|A8SD86_9FIRM Unchara  29.6 1.3E+02 0.00057   28.2   0.0   35  270-304    66-102 (117)
254 tr|A5BX84|A5BX84_VITVI Putativ  29.0 1.3E+02 0.00059   32.6   0.0   28  226-253   298-325 (350)
255 tr|R7E9K0|R7E9K0_9BACE Unchara  28.7 1.3E+02  0.0006   27.0   0.0   34  168-201    45-78  (81)
256 tr|K8XZW9|K8XZW9_RHOOP Unchara  27.9 1.4E+02 0.00064   28.9   0.0   67  192-261    21-98  (148)
257 tr|T1DKS2|T1DKS2_ANOAQ Putativ  27.8 1.3E+02 0.00064   28.9   0.0   31  261-294    57-87  (123)
258 tr|L1JXI5|L1JXI5_GUITH Unchara  27.4 1.5E+02 0.00066   34.6   0.0   51  226-277   615-665 (697)
259 tr|A0A0R1QJV5|A0A0R1QJV5_9LACO  27.1 1.5E+02 0.00067   25.9   0.0   37  178-215    11-47  (70)
260 tr|M9M860|M9M860_PSEA3 Unchara  26.0 1.6E+02 0.00072   36.7   0.0  170   64-259  1183-1359(1514)
261 tr|A0A0J9SNC8|A0A0J9SNC8_PLAVI  25.9 1.6E+02 0.00073   27.5   0.0   24  151-174     2-25  (115)
262 tr|F8NHN4|F8NHN4_SERL9 Putativ  25.7 1.7E+02 0.00074   31.6   0.0   29  297-325     4-32  (317)
263 tr|L8WKL1|L8WKL1_THACA CK1/CK1  25.3 1.7E+02 0.00076   34.6   0.0   37  183-224    50-86  (797)
264 tr|H9J8U5|H9J8U5_BOMMO Unchara  25.2 1.7E+02 0.00076   28.4   0.0   22  175-196    33-54  (147)
265 tr|Q22DK8|Q22DK8_TETTS Unchara  23.2   2E+02 0.00088   32.7   0.0   32    1-33    225-256 (497)
266 tr|C0JP56|C0JP56_RHYSE MAT1-1-  22.6 2.1E+02 0.00093   29.3   0.0   28  208-235   131-158 (204)
267 tr|W9DTD3|W9DTD3_METTI PDK rep  22.4 2.1E+02 0.00094   34.1   0.0   76  213-288   189-269 (822)
268 tr|A0A0C1D7Y4|A0A0C1D7Y4_9FLAO  22.2 2.1E+02 0.00095   25.8   0.0   35  119-153    37-71  (90)
269 tr|R9NAC4|R9NAC4_9FIRM Unchara  22.1 2.1E+02 0.00096   27.2   0.0   30   12-41      3-35  (124)
270 tr|G7UUD7|G7UUD7_PSEUP Elongat  22.0 2.1E+02 0.00097   33.2   0.0   29  284-312   495-523 (529)
271 tr|A0A061D193|A0A061D193_BABBI  21.2 2.3E+02   0.001   27.0   0.0   49  181-229    25-73  (125)
272 tr|A0A0C9N3C5|A0A0C9N3C5_9FUNG  20.2 2.5E+02  0.0011   33.3   0.0   26  284-309   507-532 (738)
273 tr|R5MDR2|R5MDR2_9BACE Unchara  20.1   2E+02  0.0011   31.1   0.0   34  168-201   197-230 (237)

No 1
>sp|Q63D03|ASNA_BACCZ Aspartate--ammonia ligase OS=Bacillus cereus (strain ZK / E33L) GN=asnA PE=3 SV=1
Probab=100.00  E-value=4.8e-146  Score=1043.25  Aligned_cols=311  Identities=62%  Similarity=1.034  Sum_probs=303.8  Template_Neff=4.200

Q d12asa_           3 IAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHD   82 (327)
Q Consensus         3 ~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~   82 (327)
                      ...|+.|.++|..|.+.|.++|.|++|+||++.+.+.|-+|||+|.|++|.+.++.+ +...|||||||||||..|+..+
T Consensus         9 ~eT~~aI~~iK~~F~~~L~~~LnL~rVsaPL~v~~~tGlNDdL~g~erpV~F~v~~~-~~~~EIVqSLAKWKR~aL~~y~   87 (327)
T sp|Q63D03|ASNA    9 IETQQAIKFIKDFFQRNLSKELNLIRVSAPLFVRKGSGLNDNLNGVERPVSFDIKDL-DATAEVVHSLAKWKRMALKRYG   87 (327)
Confidence            456889999999999999999999999999999999999999999999999999998 8999999999999999999999


Q d12asa_          83 FSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQ  161 (327)
Q Consensus        83 ~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~  161 (327)
                      |..|+|+||.|.|+|+|| .++++||+||||||||+|+...+|....||+||+.|+..|+.||..+.+.|. |.+.||++
T Consensus        88 ~~~geGi~TdmnaIRrDe-~ld~~HS~yVDQwDWE~vI~~~dR~l~~Lk~~V~kIy~~l~~te~~v~~~yp~l~~~Lp~~  166 (327)
T sp|Q63D03|ASNA   88 FSPGEGLYTDMNAIRRDE-ELDNIHSVYVDQWDWEKVISKEDRNLDYLKETVRKIYKAIKETEKEVSEKYPQLKPFLPEE  166 (327)
Confidence            999999999999999999 7999999999999999999999999999999999999999999999999986 67799999


Q d12asa_         162 IHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLE  241 (327)
Q Consensus       162 I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~  241 (327)
                      |+|+|+|||..+||+|++|+||++|+|..||||+.|||++|+||..||.||||||||+         |||||+||||+++
T Consensus       167 I~FIts~eL~~~YP~Lt~kERE~~i~ke~gAvFi~gIG~~L~~g~~Hd~RapDYDDW~---------LNGDilv~n~~l~  237 (327)
T sp|Q63D03|ASNA  167 ITFITSQELEDRYPDLTPKERENAIAKEYGAVFLIGIGGKLSDGKPHDGRAPDYDDWS---------LNGDILVWNPVLE  237 (327)
Confidence            9999999999999999999999999999999999999999999999999999999999         9999999999999


Q d12asa_         242 DAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRE  321 (327)
Q Consensus       242 ~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~  321 (327)
                      +|+||||||||||+++|+.||+++|.++|+++.||+.++.|++|+|||||||||||.|+||+-.|||.||+++||....+
T Consensus       238 ~a~ElSSMGIRVD~~sL~~Ql~~~~~~~r~~~~yH~~ll~~~LP~TIGGGIGQSRl~M~lL~K~HIGEVQ~SvW~~e~~~  317 (327)
T sp|Q63D03|ASNA  238 KAFELSSMGIRVDEEALKRQLKITGCEDRLELPFHKMLLNGELPLTIGGGIGQSRLCMFLLRKAHIGEVQASVWPDEMRE  317 (327)
Confidence            99999999999999999999999999999999999999999999999999999999999999999999999999998877


Q d12asa_         322 SVP  324 (327)
Q Consensus       322 ~~~  324 (327)
                      ...
T Consensus       318 ~~~  320 (327)
T sp|Q63D03|ASNA  318 ECE  320 (327)
Confidence            653


No 2
>tr|A0A090P0W1|A0A090P0W1_9VIBR Aspartate-ammonia ligase OS=Vibrio ponticus GN=JCM19238_71 PE=3 SV=1
Probab=100.00  E-value=9.5e-116  Score=793.48  Aligned_cols=222  Identities=66%  Similarity=1.106  Sum_probs=219.2  Template_Neff=3.452

Q d12asa_          93 MKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEF-GLAPFLPDQIHFVHSQELL  171 (327)
Q Consensus        93 MnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y-~l~~~Lp~~I~FI~sqeL~  171 (327)
                      |+|+||| |+++++||+||||||||+||+..+|...+||.||+.||++|++||..|.++| ++.|+||++|+|||||||+
T Consensus         1 m~air~d-~~ld~~HS~yVdQwDWEkvi~~e~Rn~d~LK~tV~~Iy~~~~~te~~v~~~yp~i~~~LP~~I~FI~SQEL~   79 (230)
T tr|A0A090P0W1|    1 MNAIRPD-DELDNIHSIYVDQWDWEKVINKEDRNIDYLKETVRKIYKAIKATEKYVYEEYPGIKPKLPDKITFIHSQELE   79 (230)
Confidence            8999999 9999999999999999999999999999999999999999999999999999 9999999999999999999


Q d12asa_         172 SRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGI  251 (327)
Q Consensus       172 ~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGi  251 (327)
                      ++||||++|+||++|+|..||||+.|||+||++|.+||.||||||||+         |||||+||||||.+|||||||||
T Consensus        80 d~YPDl~pk~RE~~i~ke~gAVFiigIG~~L~~G~~Hd~RApDYDDW~---------LNGDIlvw~~vl~~a~ElSSMGI  150 (230)
T tr|A0A090P0W1|   80 DRYPDLSPKERENAIAKEHGAVFIIGIGGKLSSGKPHDGRAPDYDDWS---------LNGDILVWNPVLDCALELSSMGI  150 (230)
Confidence            999999999999999999999999999999999999999999999999         99999999999999999999999


Q d12asa_         252 RVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVP  324 (327)
Q Consensus       252 rVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~~~  324 (327)
                      |||+++|++||..+|+++|++++||+++|.|++|+|||||||||||+|+||+-.|||.||+++||..+++...
T Consensus       151 RVD~~sL~~Ql~~~g~~~r~~l~~Hk~ll~g~lP~TIGGGIGQSRl~MflL~k~HIGEVQ~svWp~~~~~~~~  223 (230)
T tr|A0A090P0W1|  151 RVDADSLKRQLKIAGCEERLELPYHKMLLNGELPQTIGGGIGQSRLCMFLLRKAHIGEVQCSVWPEEVREECK  223 (230)
Confidence            9999999999999999999999999999999999999999999999999999999999999999999987654


No 3
>tr|S9W7Y1|S9W7Y1_9TRYP Aspartate--ammonia ligase OS=Angomonas deanei GN=AGDE_07378 PE=4 SV=1
Probab=100.00  E-value=1.3e-101  Score=703.30  Aligned_cols=240  Identities=82%  Similarity=1.251  Sum_probs=229.3  Template_Neff=2.856

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      +||.+|+||+|||++|+++|+++|+||+|+|||++|+|+||||||||.||+||++||.+|++.||||||||||||+|||+
T Consensus         6 ~~~~tq~ai~~vk~~f~~~l~~~LnL~~v~~Pl~~~~~~G~nDnL~G~ek~V~~~vk~~~~~~~evVhSlAKWKR~~L~~   85 (253)
T tr|S9W7Y1|S9W7    6 AYIETQRQISFVKSHFSRELEKRLNLIRVSAPILSRVGDGTQDNLSGCEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   85 (253)
Confidence            57899999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEF-GLAPFLP  159 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y-~l~~~Lp  159 (327)
                      |+|..|+||||||+|+|+|| +|+++||+||||||||+||+++||.+.+||+||+.||+.||+||.+|+++| +|.|+||
T Consensus        86 y~f~~~~GlyT~M~AlR~de-~Ld~~HS~YVDQWDWE~Vi~~eeRn~~~Lk~tV~~iy~~lk~te~~V~~~fp~l~~~LP  164 (253)
T tr|S9W7Y1|S9W7   86 YGFSAGEGLYTHMKALRPDE-RLDNLHSVYVDQWDWERVMGDEERNFSYLKSTVRKIYAAIKATEAEVHKKFPGLAPFLP  164 (253)
Confidence            99999999999999999999 899999999999999999999999999999999999999999999999999 9999999


Q d12asa_         160 DQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGL-------NGD  232 (327)
Q Consensus       160 ~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gL-------NGD  232 (327)
                      ++|||+|+|+|+.+||||++|+|||+|||.+|||||+||||||+||++||.|||||||||||++...+..       ...
T Consensus       165 ~~I~Fi~s~~l~~~ypdl~~KeRE~~iak~~gAVFl~~IG~kL~~G~pHD~RApDYDDWs~~~~~~~~~~~~~~~~~~~~  244 (253)
T tr|S9W7Y1|S9W7  165 DQIHFIHSQELLDRYPDLDPKERERAIAKELGAVFLMGIGGKLSDGHPHDGRAPDYDDWSSPSEMSSSDIGFPCAEPPME  244 (253)
Confidence            9999999999999999999999999999999999999999999999999999999999999998633221       111


Q d12asa_         233 ILVWNPVLE  241 (327)
Q Consensus       233 ilv~n~~l~  241 (327)
                      -.||||+|.
T Consensus       245 ~~~~~p~l~  253 (253)
T tr|S9W7Y1|S9W7  245 DSVWNPVLD  253 (253)
Confidence            239999874


No 4
>tr|R5AR36|R5AR36_9FIRM Aspartate--ammonia ligase OS=Firmicutes bacterium CAG:103 GN=BN455_00724 PE=3 SV=1
Probab=100.00  E-value=8.6e-93  Score=643.48  Aligned_cols=261  Identities=47%  Similarity=0.807  Sum_probs=239.2  Template_Neff=1.000

Q d12asa_           5 KQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS   84 (327)
Q Consensus         5 Tq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~   84 (327)
                      .|..|..+|..|..||...|.|..|.||++-.-..|-.|||.|-|+.|-..+... .-.-.||.|||||||..|-..||-
T Consensus        18 tqkaigllkrlfedqlaaklnlfrvsaplfleeasglndnlngyerpvlfdipqa-gkeaqvvqslakwkrmalhrydfy   96 (355)
T tr|R5AR36|R5AR   18 TQKAIGLLKRLFEDQLAAKLNLFRVSAPLFLEEASGLNDNLNGYERPVLFDIPQA-GKEAQVVQSLAKWKRMALHRYDFY   96 (355)
Confidence            4778999999999999999999999999999999999999999999987654321 224568999999999999999999


Q d12asa_          85 AGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATE---AAVSEEFGLAPFLPDQ  161 (327)
Q Consensus        85 ~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te---~~v~~~y~l~~~Lp~~  161 (327)
                      .|.||||.|.++|.|||-|..||||||||||||.+.....|...-|||||--|.|.+--|.   .|+--....-|-|..|
T Consensus        97 pgkglytdmnairrdedvldnlhsvyvdqwdwekiiessdrnldylkstvmdivaavcdtqrtmraiypqlqvlpelerq  176 (355)
T tr|R5AR36|R5AR   97 PGKGLYTDMNAIRRDEDVLDNLHSVYVDQWDWEKIIESSDRNLDYLKSTVMDIVAAVCDTQRTMRAIYPQLQVLPELERQ  176 (355)
Confidence            9999999999999999999999999999999999999999999999999999998876554   4555566677888899


Q d12asa_         162 IHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLE  241 (327)
Q Consensus       162 I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~  241 (327)
                      +.||..|||-.|||||..|+||.|..+.....|++||||.|..|..||-|.||||||.         .||||+.||..|.
T Consensus       177 vtfvtaqeledrypdltpkerehafvrehhttfiigiggalrsgkphdgrspdyddwt---------mngdilfwnelld  247 (355)
T tr|R5AR36|R5AR  177 VTFVTAQELEDRYPDLTPKEREHAFVREHHTTFIIGIGGALRSGKPHDGRSPDYDDWT---------MNGDILFWNELLD  247 (355)
Confidence            9999999999999999999999999999999999999999999999999999999996         6999999999999


Q d12asa_         242 DAFELSSMGIRVDADTLKHQLALTGDEDRLELEW  275 (327)
Q Consensus       242 ~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~  275 (327)
                      -||||||||||||++.|..||.+..-.||-|-.+
T Consensus       248 cafelssmgirvdpesldrqltiadcddrrerty  281 (355)
T tr|R5AR36|R5AR  248 CAFELSSMGIRVDPESLDRQLTIADCDDRRERTY  281 (355)
Confidence            9999999999999999999999988888776543


No 5
>tr|A0A0R2EZP5|A0A0R2EZP5_LACFE Aspartate--ammonia ligase OS=Lactobacillus fermentum GN=IV46_GL001416 PE=4 SV=1
Probab=100.00  E-value=3.7e-80  Score=550.48  Aligned_cols=229  Identities=50%  Similarity=0.827  Sum_probs=215.5  Template_Neff=1.000

Q d12asa_          26 LIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSP  105 (327)
Q Consensus        26 L~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld~  105 (327)
                      |-...||.+-.-..|-.|||...|..|....|.||.-..|+|||||||||-.|....|..-|||||.|.|+|.||| |..
T Consensus         3 lqrmsapmfveqstglndnlnrveapvsftmkdlpgetieivhslakwkrialkkygfglheglytnmnairkded-len   81 (293)
T tr|A0A0R2EZP5|    3 LQRMSAPMFVEQSTGLNDNLNRVEAPVSFTMKDLPGETIEIVHSLAKWKRIALKKYGFGLHEGLYTNMNAIRKDED-LEN   81 (293)
Confidence            3445688888888899999999999999999999999999999999999999999999999999999999999998 678


Q d12asa_         106 LHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAP-FLPDQIHFVHSQELLSRYPDLDAKGRER  184 (327)
Q Consensus       106 ~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l~~-~Lp~~I~FI~sqeL~~~YP~LtpkeRE~  184 (327)
                      .||+||||||||.|....||.-.|||.||..++.-||-.|-.|--.|--|- .|||.||||..|||-.|+|.|..-.||.
T Consensus        82 fhsiyvdqwdwekviakeerteatlkatvrqvfkvikhmehevwynfpqavyhlpdeihfvttqeledrwpelspmered  161 (293)
T tr|A0A0R2EZP5|   82 FHSIYVDQWDWEKVIAKEERTEATLKATVRQVFKVIKHMEHEVWYNFPQAVYHLPDEIHFVTTQELEDRWPELSPMERED  161 (293)
Confidence            999999999999999999999999999999999999999998887776654 4899999999999999999999999999


Q d12asa_         185 AIAKDLGAVFLVGIGGKLS-DGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLA  263 (327)
Q Consensus       185 ~i~ke~gAvFi~gIG~~L~-~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~  263 (327)
                      .||..||+||+.-||.||. .|..||-.||||||||         |||||+.|-.-|..-.|.||||||||...+|.||.
T Consensus       162 kiaeelgcvfvmkigdklqrsgephdgyapdyddws---------lngdiifwyeplqtklevssmgirvderamkeqle  232 (293)
T tr|A0A0R2EZP5|  162 KIAEELGCVFVMKIGDKLQRSGEPHDGYAPDYDDWS---------LNGDIIFWYEPLQTKLEVSSMGIRVDERAMKEQLE  232 (293)
Confidence            9999999999999999985 6889999999999997         79999999999999999999999999999999986


Q d12asa_         264 L  264 (327)
Q Consensus       264 ~  264 (327)
                      -
T Consensus       233 k  233 (293)
T tr|A0A0R2EZP5|  233 K  233 (293)
Confidence            4


No 6
>tr|W1BHY9|W1BHY9_KLEPN Aspartate--ammonia ligase OS=Klebsiella pneumoniae IS22 PE=3 SV=1
Probab=100.00  E-value=5.9e-80  Score=544.20  Aligned_cols=194  Identities=85%  Similarity=1.284  Sum_probs=192.2  Template_Neff=2.794

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      +||++|++|+|||++|+|+|+++|||+|||+||++|+|||+||||+|.||+|+++||.+||++||||||||||||+||++
T Consensus         4 s~~~~q~~i~~vk~~f~r~l~~~l~l~~V~~Pi~~r~~~G~~DnLnG~E~pV~~~vk~~~~~~~eVVHSLAKWKR~~L~~   83 (197)
T tr|W1BHY9|W1BH    4 SYIAKQRAISFVKSHFERQLEKKLGLIEVQAPILSRVGDGTQDNLNGCEKPVQVKVKTLPDAQFEVVHSLAKWKRQTLGQ   83 (197)
Confidence            47999999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPD  160 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l~~~Lp~  160 (327)
                      ++|++|+||||||||+|||||+++|+||+||||||||+||+|++|++.+||.||++||+|||+||.+|+++|++.|||||
T Consensus        84 ~~f~~geGl~~~M~aiR~dED~~~~~HS~yVdQwDWEkVi~~e~Rni~~LK~tV~~IY~~ik~te~~V~~~y~l~p~Lp~  163 (197)
T tr|W1BHY9|W1BH   84 HGFSAGEGLYTHMKALRPDEDRLDPIHSVYVDQWDWEKVIPDEERNIETLKETVEAIYAGIKATELAVSAEYGLEPFLPD  163 (197)
Confidence            99999999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_         161 QIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVF  194 (327)
Q Consensus       161 ~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvF  194 (327)
                      ||+|+|||||+++||||++++|+++|+|..|++|
T Consensus       164 qI~fihsqeL~~~~pdl~~k~~~~~i~k~~~~~~  197 (197)
T tr|W1BHY9|W1BH  164 QIHFIHSQELKSRYPDLDAKGRERAIPKELGAVL  197 (197)
Confidence            9999999999999999999999999999999986


No 7
>tr|T0E9S1|T0E9S1_CLOSO Aspartate-ammonia ligase family protein OS=[Clostridium] sordellii ATCC 9714 GN=H477_1191 PE=3 SV=1
Probab=100.00  E-value=4.5e-66  Score=427.32  Aligned_cols=118  Identities=65%  Similarity=1.092  Sum_probs=115.1  Template_Neff=2.537

Q d12asa_         197 GIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWH  276 (327)
Q Consensus       197 gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h  276 (327)
                      .||..++.|..||-||||||||.         |||||+||||+|..+|||||||||||+++|++||..+|.++|++++||
T Consensus         2 ~~~~~l~~g~~HdgRapDYDDW~---------LNGDilv~~~vL~~~lElSSMGIRVd~~sL~~QL~~~g~~~r~~~~fh   72 (127)
T tr|T0E9S1|T0E9    2 QIGKMLSSGEPHDGRAPDYDDWE---------LNGDILVYNPVLDIALELSSMGIRVDEKSLKRQLEIAGCEERLELDFH   72 (127)
Confidence            57889999999999999999997         999999999999999999999999999999999999999999999999


Q d12asa_         277 QALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESV  323 (327)
Q Consensus       277 ~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~~  323 (327)
                      ++++.|+.|.||||||||||+.|++|+-.|||.||+++||.++|+..
T Consensus        73 k~~l~~~lP~TiGGGIGqSRicMf~L~kaHIGEVq~S~Wp~~~~~~~  119 (127)
T tr|T0E9S1|T0E9   73 KALLNGELPLTIGGGIGQSRICMFFLRKAHIGEVQASVWPDEMREEC  119 (127)
Confidence            99999999999999999999999999999999999999999999864


No 8
>tr|B7XN60|B7XN60_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_25239 PE=4 SV=1
Probab=100.00  E-value=4.9e-66  Score=446.24  Aligned_cols=145  Identities=46%  Similarity=0.747  Sum_probs=135.4  Template_Neff=2.631

Q d12asa_           6 QRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSA   85 (327)
Q Consensus         6 q~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~   85 (327)
                      |.-|..+|..|++.|-+.|+|..|.|||.-+-..|..|||+|. +.|...++.-  ...|||||||||||-.|....|+.
T Consensus        11 q~aIk~~~d~fe~~l~~~l~L~rvsaP~fv~~~sglnD~LnG~-rpV~fd~~~~--~~~eivHSLAKWKR~aL~kygf~~   87 (170)
T tr|B7XN60|B7XN   11 QEAIKTIKDFFERELGEELNLTRVSAPIFVRPESGLNDNLNGV-RPVRFDIKDG--EWAEIVHSLAKWKRMALYKYGFQI   87 (170)
Confidence            5568889999999999999999999999999999999999999 9998888876  778999999999999999999999


Q d12asa_          86 GEGLYTHMKALRPDEDR----LSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG  153 (327)
Q Consensus        86 geGiyTdMnAIRrDE~~----ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~  153 (327)
                      ++|+||.|.|.|.|||-    -|.+||+||||||||.|....+|...+||.||+.|..++++|+..|-+++.
T Consensus        88 ~~GiytdM~aIRrdEd~D~~~~~~iHS~yVdQwDWEkvI~ke~Rnid~Lk~tV~~Iy~~~~~~~~~~~~~y~  159 (170)
T tr|B7XN60|B7XN   88 GEGIYTDMNAIRRDEDLDNICTSRIHSYYVDQWDWEKVIDKEDRNIDFLKETVRKIYSAFKKTEEYVYEKYP  159 (170)
Confidence            99999999999999985    346999999999999999999999999999999999999999998877653


No 9
>tr|B7XPB1|B7XPB1_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_24385 PE=3 SV=1
Probab=100.00  E-value=1.2e-59  Score=413.50  Aligned_cols=237  Identities=42%  Similarity=0.724  Sum_probs=210.9  Template_Neff=1.000

Q d12asa_          40 GTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR----LSPLHSVYVDQWD  115 (327)
Q Consensus        40 GlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~----ld~~HSiyVDQWD  115 (327)
                      |..|.|.|.|.   |+......--..||||||||||-.|-.-.-..-.|++..|+.+|.-||.    .|..||.-|.|+|
T Consensus         3 gvnddlngeep---vrfqtsdgewcsvvhslakwkrwmlwklrddgvsgiwcdmrgirkcedadvictsrmhsyqveqfd   79 (255)
T tr|B7XPB1|B7XP    3 GVNDDLNGEEP---VRFQTSDGEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRGIRKCEDADVICTSRMHSYQVEQFD   79 (255)
Confidence            44566777653   3344444555679999999999888665555667899999999987764    5889999999999


Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG---LAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGA  192 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~---l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gA  192 (327)
                      ||.|.....|...-||.||..|+.|+||....|..|+.   |..-|.|.+.||.||||-.-+|.|+|..||.|||+..|.
T Consensus        80 wekvipsekrnieflketvreiyrglkaakervdaeypgmflrnkladevrfvqsqeleelfpelnaqerenaiareags  159 (255)
T tr|B7XPB1|B7XP   80 WEKVIPSEKRNIEFLKETVREIYRGLKAAKERVDAEYPGMFLRNKLADEVRFVQSQELEELFPELNAQERENAIAREAGS  159 (255)
Confidence            99999999999999999999999999999988887763   556688999999999999999999999999999999999


Q d12asa_         193 VFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLE  272 (327)
Q Consensus       193 vFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~  272 (327)
                      ||++|||-.|..|..|.-|..|||||.         ||||++|||.|+--..|.||||||||-..|+.||.........|
T Consensus       160 vfiigighplksgkehgqrsadyddwd---------lngdlivwnevigcgmevssmgirvdhtslrrqleaknesrkae  230 (255)
T tr|B7XPB1|B7XP  160 VFIIGIGHPLKSGKEHGQRSADYDDWD---------LNGDLIVWNEVIGCGMEVSSMGIRVDHTSLRRQLEAKNESRKAE  230 (255)
Confidence            999999999999999999999999995         89999999999999999999999999999999999988888899


Q d12asa_         273 LEWHQALLRGEMPQTI  288 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TI  288 (327)
                      .|+|...-.|-+|.|.
T Consensus       231 yefhrgvpegvipptm  246 (255)
T tr|B7XPB1|B7XP  231 YEFHRGVPKECSPPTM  246 (255)
Confidence            9999999999888874


No 10
>tr|B7XQP4|B7XQP4_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_21779 PE=4 SV=1
Probab=100.00  E-value=4.4e-56  Score=380.73  Aligned_cols=190  Identities=43%  Similarity=0.760  Sum_probs=165.9  Template_Neff=1.000

Q d12asa_          40 GTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR----LSPLHSVYVDQWD  115 (327)
Q Consensus        40 GlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~----ld~~HSiyVDQWD  115 (327)
                      |..|.|.|.|.   |+......--..||||||||||-.|-.-.-..-.|++..|+.+|.-||.    .|..||.-|.|+|
T Consensus         3 gvnddlngeep---vrfqtsdgewcsvvhslakwkrwmlwklrddgvsgiwcdmrgirkcedadvictsrmhsyqveqfd   79 (203)
T tr|B7XQP4|B7XQ    3 GVNDDLNGEEP---VRFQTSDGEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRGIRKCEDADVICTSRMHSYQVEQFD   79 (203)
Confidence            44566777653   3344444555679999999999888665555667899999999987764    5889999999999


Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG---LAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGA  192 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~---l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gA  192 (327)
                      ||.|.....|...-||.||..|..|+||....|..|+.   |..-|.|.+.||.||||-.-+|.|++..||.|||+..|.
T Consensus        80 wekvipsekrnieflketvreiyrglkaakervdaeypgmflrnkladevrfvqsqeleelfpelnaqerenaiareags  159 (203)
T tr|B7XQP4|B7XQ   80 WEKVIPSEKRNIEFLKETVREIYRGLKAAKERVDAEYPGMFLRNKLADEVRFVQSQELEELFPELNAQERENAIAREAGS  159 (203)
Confidence            99999999999999999999999999999988887763   556688999999999999999999999999999999999


Q d12asa_         193 VFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLE  241 (327)
Q Consensus       193 vFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~  241 (327)
                      ||++|||-.|..|..|.-|..|||||.         ||||++|||.|..
T Consensus       160 vfiigighplksgkehgqrsadyddwd---------lngdlivwnevmd  199 (203)
T tr|B7XQP4|B7XQ  160 VFIIGIGHPLKSGKEHDQSHADYDDWD---------LNGDLIVWNEVMD  199 (203)
Confidence            999999999999999999999999995         8999999999875


No 11
>tr|A0A074ZZK9|A0A074ZZK9_9TREM Uncharacterized protein (Fragment) OS=Opisthorchis viverrini GN=T265_12984 PE=3 SV=1
Probab=100.00  E-value=8.2e-55  Score=415.06  Aligned_cols=304  Identities=47%  Similarity=0.811  Sum_probs=287.6  Template_Neff=1.000

Q d12asa_           7 RQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAG   86 (327)
Q Consensus         7 ~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~g   86 (327)
                      ..|.++|..|...|...|.|+.|.||++-.---|-||.|||.|+|+-..|..-..+-  +-.|||||||-.||.....--
T Consensus       270 egikyikdqfqqalatnlnllrvsaplfvpnrlglqddlsgveraiyfdvrsgeeav--inqslakwkrmalgkyqlkry  347 (583)
T tr|A0A074ZZK9|  270 EGIKYIKDQFQQALATNLNLLRVSAPLFVPNRLGLQDDLSGVERAIYFDVRSGEEAV--INQSLAKWKRMALGKYQLKRY  347 (583)
Confidence            468899999999999999999999999988888999999999999998888765553  457999999999999999999


Q d12asa_          87 EGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGL-APFLPDQIHFV  165 (327)
Q Consensus        87 eGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l-~~~Lp~~I~FI  165 (327)
                      .||||.|.|+|.|| -.|||||.||||||||-|....||.-..|...||.|+...|.|||||-..+.+ ..-||..|.|+
T Consensus       348 kglytdmnairrde-iisplhsyyvdqwdwemvierderterkleevvekiytafkqteaavlakypifskklpekitfi  426 (583)
T tr|A0A074ZZK9|  348 KGLYTDMNAIRRDE-IISPLHSYYVDQWDWEMVIERDERTERKLEEVVEKIYTAFKQTEAAVLAKYPIFSKKLPEKITFI  426 (583)
Confidence            99999999999987 58999999999999999999999999999999999999999999999877754 45699999999


Q d12asa_         166 HSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFE  245 (327)
Q Consensus       166 ~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~E  245 (327)
                      ..|||-..||||....||-...|..|||||--||.+|..|..||-||||||||.         |||||+.++.......|
T Consensus       427 rtqeledqypdlkpsereyqavkkygavflkqigkrlksgeihdsrapdyddwe---------lngdilfydanadrcie  497 (583)
T tr|A0A074ZZK9|  427 RTQELEDQYPDLKPSEREYQAVKKYGAVFLKQIGKRLKSGEIHDSRAPDYDDWE---------LNGDILFYDANADRCIE  497 (583)
Confidence            999999999999999999999999999999999999999999999999999995         79999999999999999


Q d12asa_         246 LSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRES  322 (327)
Q Consensus       246 lSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~  322 (327)
                      |||||||||.-.+..||.-+|-|.|-..|+||..|.||.|-||||||||||+.|.+|.-.|+|.||...||..||..
T Consensus       498 lssmgirvdevsmekqlresgcesrkqfeyhqkvlngeypltigggigqsricmfflekqhvgevqcsiwpehvrqe  574 (583)
T tr|A0A074ZZK9|  498 LSSMGIRVDEVSMEKQLRESGCESRKQFEYHQKVLNGEYPLTIGGGIGQSRICMFFLEKQHVGEVQCSIWPEHVRQE  574 (583)
Confidence            99999999999999999999999999999999999999999999999999999999999999999999999999864


No 12
>tr|K1U8F8|K1U8F8_9ZZZZ Aspartate--ammonia ligase OS=human gut metagenome GN=LEA_03427 PE=4 SV=1
Probab=100.00  E-value=7.9e-53  Score=354.84  Aligned_cols=141  Identities=43%  Similarity=0.725  Sum_probs=115.2  Template_Neff=1.000

Q d12asa_          38 GDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWE  117 (327)
Q Consensus        38 ~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWE  117 (327)
                      |.|-.|.|.|.|+.|...|..| |-.-|||||||||||-.|....|..|.||.|.|.|+|.||. |..|||+||||||||
T Consensus         3 gsglnddlngverpvsfdvpcl-deraevvhslakwkryalaeygfrpgqglvtdmnairrdee-ldnlhsiyvdqwdwe   80 (178)
T tr|K1U8F8|K1U8    3 GSGLNDDLNGVERPVSFDVPCL-DERAEVVHSLAKWKRYALAEYGFRPGQGLVTDMNAIRRDEE-LDNLHSIYVDQWDWE   80 (178)
Confidence            6788999999999999998877 45789999999999999999999999999999999999875 788999999999999


Q d12asa_         118 RVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG--LAPFLPDQIHFVHSQELLSRYPDLDAK  180 (327)
Q Consensus       118 kvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~--l~~~Lp~~I~FI~sqeL~~~YP~Ltpk  180 (327)
                      .|.....|...-|..||.-|...+-.....+.=.|-  -+--|...-.|+..|||-.-||||...
T Consensus        81 kvitakdrtlpflqetvrdivdavcsaadelrwkfpelkairltreptfittqeledlypdltpq  145 (178)
T tr|K1U8F8|K1U8   81 KVITAKDRTLPFLQETVRDIVDAVCSAADELRWKFPELKAIRLTREPTFITTQELEDLYPDLTPQ  145 (178)
Confidence            999999999988999888776544322222221221  012234455799999999999999764


No 13
>tr|K1TDC4|K1TDC4_9ZZZZ Aspartate/ammonia ligase (Fragment) OS=human gut metagenome GN=OBE_05291 PE=4 SV=1
Probab=100.00  E-value=3.4e-52  Score=331.58  Aligned_cols=118  Identities=63%  Similarity=1.049  Sum_probs=115.5  Template_Neff=1.000

Q d12asa_         174 YPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRV  253 (327)
Q Consensus       174 YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirV  253 (327)
                      ||||..|.||..+.+.-.||||.-||-.|..|.|||-||||||||.         +||||+|+-|||.-|+|||||||||
T Consensus         2 ypdlepkerenkllrekkavflmqigkelksgqrhdgrapdyddwe---------lngdilvyypvldmayelssmgirv   72 (119)
T tr|K1TDC4|K1TD    2 YPDLEPKERENKLLREKKAVFLMQIGKELKSGQRHDGRAPDYDDWE---------LNGDILVYYPVLDMAYELSSMGIRV   72 (119)
Confidence            8999999999999999999999999999999999999999999995         7999999999999999999999999


Q d12asa_         254 DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTML  300 (327)
Q Consensus       254 d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~  300 (327)
                      |.|.|-.|+...|-||||||.+|..|+..|.|-|||||||||||.|.
T Consensus        73 dedsllsqikiagcedrleldfhkkllnkelpytigggigqsrlcmy  119 (119)
T tr|K1TDC4|K1TD   73 DEDSLLSQIKIAGCEDRLELDFHKKLLNKELPYTIGGGIGQSRLCMY  119 (119)
Confidence            99999999999999999999999999999999999999999999873


No 14
>tr|A0A0P8SA41|A0A0P8SA41_ENTCL Aspartate--ammonia ligase (Fragment) OS=Enterobacter cloacae subsp. cloacae GN=AN693_29235 PE=4 SV=1
Probab=100.00  E-value=5.3e-49  Score=317.57  Aligned_cols=124  Identities=69%  Similarity=1.167  Sum_probs=120.6  Template_Neff=1.280

Q d12asa_         105 PLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRER  184 (327)
Q Consensus       105 ~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~  184 (327)
                      |.||||||||||||||++|+|...+||+||+||.|||++||+|||.+||||||||..|||||||+|++|+|++|+|+|||
T Consensus         1 ~~hsvyvdqwdwe~v~~~~~r~l~~lk~tv~aiya~i~~te~avsk~fgla~flp~~i~fvhs~~l~~rfp~~~~k~re~   80 (124)
T tr|A0A0P8SA41|    1 QIHSVYVDQWDWERVMSDEQRHLGYLKSTVRAIYAGIKETEEAVSKKFGLAPFLPKDIQFVHSQELVQRFPNMNDKEREN   80 (124)
Confidence            57999999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_         185 AIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAG  228 (327)
Q Consensus       185 ~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~g  228 (327)
                      ||+|.+|||||+|||||||||.+||||+||||||||..|-..||
T Consensus        81 ai~ke~gavfligiggklsdgk~hdvrapdyddwstagehdyag  124 (124)
T tr|A0A0P8SA41|   81 AICKEYGAVFLIGIGGKLSDGKPHDVRAPDYDDWSTAGEHDYAG  124 (124)
Confidence            99999999999999999999999999999999999998866554


No 15
>tr|J9G8M4|J9G8M4_9ZZZZ Aspartate-ammonia ligase OS=gut metagenome GN=EVA_08334 PE=4 SV=1
Probab=100.00  E-value=4e-48  Score=311.06  Aligned_cols=112  Identities=50%  Similarity=0.840  Sum_probs=106.8  Template_Neff=1.180

Q d12asa_           8 QISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGE   87 (327)
Q Consensus         8 aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~ge   87 (327)
                      .|..+|..|...|...|.|-.|.||+.---|-|..|.|.|.|++|...+|.+.+||-|||||||||||-||.......+-
T Consensus         9 gik~ik~ffq~nls~elrlrrvtaplfvl~g~ginddlngver~v~fpikdlg~a~aevvhslakwkrltla~y~iep~y   88 (122)
T tr|J9G8M4|J9G8    9 GIKMIKDFFQMNLSSELRLRRVTAPLFVLQGMGINDDLNGVERPVSFPIKDLGDAQAEVVHSLAKWKRLTLADYNIEPGY   88 (122)
Confidence            57778999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_          88 GLYTHMKALRPDEDRLSPLHSVYVDQWDWERVM  120 (327)
Q Consensus        88 GiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI  120 (327)
                      |.||.|.|+|.||. |..|||.|||||||+.+.
T Consensus        89 g~ytdmnairadee-l~nlhslyvdqwdw~~~~  120 (122)
T tr|J9G8M4|J9G8   89 GIYTDMNAIRADEE-LDNLHSLYVDQWDWEKVY  120 (122)
Confidence            99999999999985 789999999999999864


No 16
>tr|X4JHK5|X4JHK5_SALEN Asparagine synthetase AsnA (Fragment) OS=Salmonella enterica subsp. enterica serovar Enteritidis str. SA20090419 GN=AU95_23665 
Probab=100.00  E-value=8.3e-47  Score=311.13  Aligned_cols=137  Identities=79%  Similarity=1.205  Sum_probs=136.1  Template_Neff=1.371

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      +||++||||||||+||+|+|++|||+||||+|||||||||+||||||.|++|||||+++|+|+|||||||||||||||||
T Consensus         4 ~~i~~q~~isfvk~~f~~~l~~~l~~~evq~pils~vgdg~qdnlsg~e~~v~vkv~~~p~a~fevvhslakwkr~tl~~   83 (140)
T tr|X4JHK5|X4JH    4 SFIHQQQQISFVKNHFSQQLIDRLEIIEVQGPILSQVGDGMQDNLSGCEHPVQVKVKNIPDAQFEVVHSLAKWKRQTLGQ   83 (140)
Confidence            69999999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAI  137 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kI  137 (327)
                      |+|++|||||+||||||||||+|+|+||||||||||||||+||+|||++||+|||+|
T Consensus        84 ~~f~~gegl~~hmkalrpded~l~p~hsvyvdqwdwe~v~~dg~r~~~~lk~tve~i  140 (140)
T tr|X4JHK5|X4JH   84 FDFNEGEGLFTHMKALRPDEDSLSPTHSVYVDQWDWERVMPDGRRQFSYLKSTVEKI  140 (140)
Confidence            999999999999999999999999999999999999999999999999999999987


No 17
>tr|S9CF61|S9CF61_STRAG Uncharacterized protein OS=Streptococcus agalactiae FSL S3-105 GN=SAG0023_02635 PE=3 SV=1
Probab=100.00  E-value=1.1e-44  Score=307.55  Aligned_cols=158  Identities=80%  Similarity=1.253  Sum_probs=155.3  Template_Neff=1.935

Q d12asa_         170 LLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSM  249 (327)
Q Consensus       170 L~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSm  249 (327)
                      ++..||||++|+||+||||++|||||+||||+||||++||+||||||||+||||.||+||||||+||||||+++||+|||
T Consensus         1 ~~~~ypdl~~kerenaiake~gavfligig~~l~~g~~hd~rapdyddw~t~~~~g~~glngdilvwn~vl~~~~elssm   80 (158)
T tr|S9CF61|S9CF    1 LEEKYPDLNAKERENAIAKEAGAVFIIGIGHPLSDGKRHDQRAPDYDDWSTESSEGFAGLNGDILVWNPVLGDAFELSSM   80 (158)
Confidence            45789999999999999999999999999999999999999999999999999999999999999999999999999999


Q d12asa_         250 GIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL  327 (327)
Q Consensus       250 GirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~~~~~~  327 (327)
                      |||||+|+||+||++||||||+|+|||+++|||+||+||||||||||++|+|||++||||||++|||++|+|++++++
T Consensus        81 girvd~~~lk~ql~~~g~edr~e~e~h~~~l~g~~p~tigggigqsr~~m~ll~~~higevq~~vwp~~v~e~~~~~~  158 (158)
T tr|S9CF61|S9CF   81 GIRVDADSLKRQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMFLLQLRHIGEVQCSVWPKKVREECAALL  158 (158)
Confidence            999999999999999999999999999999999999999999999999999999999999999999999999999864


No 18
>tr|A0A0E9FL15|A0A0E9FL15_CHLTH Aspartate--ammonia ligase OS=Chlamydia trachomatis GN=asnA_2 PE=4 SV=1
Probab=100.00  E-value=1.9e-43  Score=276.74  Aligned_cols=78  Identities=81%  Similarity=1.273  Sum_probs=76.0  Template_Neff=3.190

Q d12asa_         249 MGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSL  326 (327)
Q Consensus       249 mGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~~~~~  326 (327)
                      ||||||+++|+.||+++|++||++|+||++|++||||+|||||||||||+|+|||.+|||+||+++||..+++.++..
T Consensus         1 MGIRVD~~~L~~QL~~~g~edR~~l~~Hk~Ll~gelP~TIGGGIGQSRl~MllLqk~HIGEVQ~svW~~e~~~~~~~~   78 (83)
T tr|A0A0E9FL15|    1 MGIRVDADALKRQLALTGDEDRLKLEWHQALLNGELPQTIGGGIGQSRLCMLLLQKAHIGEVQCSVWPAEVRESCPAI   78 (83)
Confidence            899999999999999999999999999999999999999999999999999999999999999999999999988753


No 19
>tr|Q3EXH4|Q3EXH4_BACTI Aspartate--ammonia ligase OS=Bacillus thuringiensis serovar israelensis ATCC 35646 GN=RBTH_05986 PE=4 SV=1
Probab=100.00  E-value=1.2e-40  Score=274.56  Aligned_cols=106  Identities=44%  Similarity=0.745  Sum_probs=94.7  Template_Neff=1.000

Q d12asa_           6 QRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSA   85 (327)
Q Consensus         6 q~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~   85 (327)
                      |-.|.-||..|..||..||.|..|.||+.-.--.|-.|.|.|.|+.+....-. ..-..|+|||||||||-.|-...+.|
T Consensus        12 qiaikevktffedqlakrlelfrvsaplfvtkksglndhlngverpiefdmlh-sgeeleivhslakwkrfalheygyea   90 (147)
T tr|Q3EXH4|Q3EX   12 QIAIKEVKTFFEDQLAKRLELFRVSAPLFVTKKSGLNDHLNGVERPIEFDMLH-SGEELEIVHSLAKWKRFALHEYGYEA   90 (147)
Confidence            44567799999999999999999999999999999999999999988765432 45678999999999999999999999


Q d12asa_          86 GEGLYTHMKALRPDEDRLSPLHSVYVDQ  113 (327)
Q Consensus        86 geGiyTdMnAIRrDE~~ld~~HSiyVDQ  113 (327)
                      ||||||.|.|+|.||. |..-||+|||.
T Consensus        91 geglytnmnairrdee-ldathsiyvdh  117 (147)
T tr|Q3EXH4|Q3EX   91 GEGLYTNMNAIRRDEE-LDATHSIYVDH  117 (147)
Confidence            9999999999999975 77889999985


No 20
>tr|B7XP24|B7XP24_ENTBH Aspartate-ammonia ligase (Fragment) OS=Enterocytozoon bieneusi (strain H348) GN=EBI_25530 PE=4 SV=1
Probab=100.00  E-value=5e-40  Score=273.66  Aligned_cols=149  Identities=40%  Similarity=0.647  Sum_probs=125.5  Template_Neff=1.000

Q d12asa_          40 GTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR----LSPLHSVYVDQWD  115 (327)
Q Consensus        40 GlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~----ld~~HSiyVDQWD  115 (327)
                      |..|.|.|.|.   |.......--..||||||||||-.|-.-.-..-.|++..|+.+|.-||.    .|..||.-|.|+|
T Consensus         3 gvnddlngeep---vrfqtsdgewcsvvhslakwkrwmlwklrddgvsgiwcdmrgirkcedadvictsrmhsyqveqfd   79 (158)
T tr|B7XP24|B7XP    3 GVNDDLNGEEP---VRFQTSDGEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRGIRKCEDADVICTSRMHSYQVEQFD   79 (158)
Confidence            45567777653   3344444555679999999999888665555667899999999987764    5889999999999


Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG---LAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLG  191 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~---l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~g  191 (327)
                      ||.|.....|...-||.||..|+.|+||....|..|+.   |..-|.|.+.||.||||-.-+|.|++..||.|||+..|
T Consensus        80 wekvipsekrnieflketvreiyrglkaakervdaeypgmflrnkladevrfvqsqeleelfpelnaqerenaiareag  158 (158)
T tr|B7XP24|B7XP   80 WEKVIPSEKRNIEFLKETVREIYRGLKAAKERVDAEYPGMFLRNKLADEVRFVQSQELEELFPELNAQERENAIAREAG  158 (158)
Confidence            99999999999999999999999999999988887763   55668899999999999999999999999999998754


No 21
>tr|G4ELX7|G4ELX7_MYCIO Asparagine synthetase AsnA (Fragment) OS=Mycoplasma iowae 695 GN=GUU_02027 PE=4 SV=1
Probab=99.95  E-value=4.5e-33  Score=221.35  Aligned_cols=92  Identities=97%  Similarity=1.345  Sum_probs=91.0  Template_Neff=2.007

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      +||++|++|||+|++|+|+|+++|||||||+||+||+|||+||||||+||||||+|| +||++||||||||||||+||++
T Consensus         7 ~~~~~q~~i~~~k~~f~~~~~~~l~l~~v~~pil~~~~~g~~dnlsg~e~~vq~~~~-~~~~~~ev~~slakwkr~~l~~   85 (98)
T tr|G4ELX7|G4EL    7 SYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGCEKAVQVKVK-LPDAQFEVVHSLAKWKRQTLGQ   85 (98)
Confidence            589999999999999999999999999999999999999999999999999999999 9999999999999999999999


Q d12asa_          81 HDFSAGEGLYTHM   93 (327)
Q Consensus        81 y~~~~geGiyTdM   93 (327)
                      +||++|||+||||
T Consensus        86 ~d~~~~eg~y~~m   98 (98)
T tr|G4ELX7|G4EL   86 YDFSAGEGLYTHM   98 (98)
Confidence            9999999999998


No 22
>tr|S7XF60|S7XF60_SPRLO Putative aspartate-ammonia ligase (Fragment) OS=Spraguea lophii (strain 42_110) GN=SLOPH_2596 PE=4 SV=1
Probab=99.94  E-value=4.6e-31  Score=222.49  Aligned_cols=107  Identities=45%  Similarity=0.705  Sum_probs=94.7  Template_Neff=1.000

Q d12asa_         152 FGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG------  225 (327)
Q Consensus       152 y~l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~------  225 (327)
                      |--..-||.||.|+.||||-.-||......||-...|--.|+|+.|||-.|..|..||.||.||||||..-|--      
T Consensus        29 fltkntlpeqitfissqeledlypnmtpsereyeyvklkkaifimgighelksgkphdlraddyddws~dwekvirkedr  108 (167)
T tr|S7XF60|S7XF   29 FLTKNTLPEQITFISSQELEDLYPNMTPSEREYEYVKLKKAIFIMGIGHELKSGKPHDLRADDYDDWSXDWEKVIRKEDR  108 (167)
Confidence            44456789999999999999999999999999999999999999999999999999999999999999754311      


Q d12asa_         226 --------------------HAGLNGDILVWNPVLEDAFELSSMGIRVDADTL  258 (327)
Q Consensus       226 --------------------~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L  258 (327)
                                          ..-.||||..+|......||.|||||||..+-|
T Consensus       109 nfeflksivrkiysavlstqkdl~ngdiilyngeidrefeissmgirv~meil  161 (167)
T tr|S7XF60|S7XF  109 NFEFLKSIVRKIYSAVLSTQKDLXNGDIILYNGEIDREFEISSMGIRVXMEIL  161 (167)
Confidence                                123599999999999999999999999987654


No 23
>tr|B7XQM8|B7XQM8_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_23995 PE=4 SV=1
Probab=99.91  E-value=1.2e-28  Score=200.41  Aligned_cols=104  Identities=37%  Similarity=0.589  Sum_probs=83.2  Template_Neff=1.000

Q d12asa_          40 GTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR----LSPLHSVYVDQWD  115 (327)
Q Consensus        40 GlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~----ld~~HSiyVDQWD  115 (327)
                      |..|.|.|.|.   |.......--..||||||||||-.|-.-.-..-.|++..|+.+|.-||.    -|..||.-|.|+|
T Consensus         3 gvnddlngeep---vrfqtsdgewcsvvhslakwkrwmlwklrddgvsgiwcdmrgirkcedadvictsrmhsyqveqfd   79 (128)
T tr|B7XQM8|B7XQ    3 GVNDDLNGEEP---VRFQTSDGEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRGIRKCEDADVICTSRMHSYQVEQFD   79 (128)
Confidence            45567777653   3344444455679999999999888665555667899999999987764    5889999999999


Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAGIKATEA  146 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~al~~te~  146 (327)
                      ||.|.....|...-||.||..|..-+|....
T Consensus        80 wekvipsekrnieflketvreiyprmkqqks  110 (128)
T tr|B7XQM8|B7XQ   80 WEKVIPSEKRNIEFLKETVREIYPRMKQQKS  110 (128)
Confidence            9999999999999999999999988776543


No 24
>tr|K6E9S2|K6E9S2_9BACI Asparagine synthetase AsnA OS=Bacillus bataviensis LMG 21833 GN=BABA_07056 PE=4 SV=1
Probab=99.91  E-value=1.5e-28  Score=177.27  Aligned_cols=55  Identities=58%  Similarity=1.059  Sum_probs=51.4  Template_Neff=1.000

Q d12asa_          64 FEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERV  119 (327)
Q Consensus        64 ~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkv  119 (327)
                      .|||.|||||||-.|....|...||||..|.|+|.||| |..-||++|||||||..
T Consensus         1 mevvqslakwkrmalarygfpvneglyidmnairrded-ldnqhsifvdqwdwetl   55 (62)
T tr|K6E9S2|K6E9    1 MEVVQSLAKWKRMALARYGFPVNEGLYIDMNAIRRDED-LDNQHSIFVDQWDWETL   55 (62)
Confidence            37999999999999999999999999999999999998 66789999999999964


No 25
>tr|W1XGP1|W1XGP1_9ZZZZ Aspartate-ammonia ligase (Fragment) OS=human gut metagenome GN=Q604_UNBC16002G0001 PE=4 SV=1
Probab=99.90  E-value=5.3e-28  Score=178.37  Aligned_cols=69  Identities=68%  Similarity=1.200  Sum_probs=66.7  Template_Neff=1.000

Q d12asa_          49 EKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWE  117 (327)
Q Consensus        49 ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWE  117 (327)
                      |..|.|||-..||-.+|||||||||||.||....|..||||+.||||||||||.|..+|||||||||||
T Consensus         2 ehpvsvkvlqipdetyevvhslakwkrhtlarfgfgegeglfvhmkalrpdedsldaihsvyvdqwdwe   70 (70)
T tr|W1XGP1|W1XG    2 EHPVSVKVLQIPDETYEVVHSLAKWKRHTLARFGFGEGEGLFVHMKALRPDEDSLDAIHSVYVDQWDWE   70 (70)
Confidence            456889999999999999999999999999999999999999999999999999999999999999997


No 26
>tr|A0A0J9TKM5|A0A0J9TKM5_PLAVI Uncharacterized protein OS=Plasmodium vivax North Korean GN=PVNG_02423 PE=4 SV=1
Probab=99.89  E-value=1.7e-27  Score=188.24  Aligned_cols=93  Identities=38%  Similarity=0.546  Sum_probs=85.2  Template_Neff=1.520

Q d12asa_         125 RQFSTLKSTVEAIWAGIKATEAAVSEE-FGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLS  203 (327)
Q Consensus       125 Rnl~~Lk~tV~kIy~al~~te~~v~~~-y~l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~  203 (327)
                      |...-|+|.|.-|+...|-||..+... -.|.--||+.|.|+.+.||+.-||++...+||-...|.-.|.|+.-||.+|+
T Consensus         5 r~~~fl~sivkki~k~~k~tek~~n~~~~~l~~klp~ni~fi~~~~l~~~y~~~s~~ere~~~vk~~ka~fiy~ig~~l~   84 (98)
T tr|A0A0J9TKM5|    5 RNEQFLFSIVKKIFKSFKETEKEFNSQYPNLNRKLPDNISFISSKDLFNMYPNKSPDEREYEYVKKHKAIFIYQIGHKLP   84 (98)
Confidence            445568999999999999999877543 3577889999999999999999999999999999999999999999999999


Q d12asa_         204 DGHRHDVRAPDYDD  217 (327)
Q Consensus       204 ~G~~Hd~RapDYDD  217 (327)
                      ||..|..||-||||
T Consensus        85 d~~~h~~ra~dydd   98 (98)
T tr|A0A0J9TKM5|   85 DKSPHQFRAFDYDD   98 (98)
Confidence            99999999999997


No 27
>tr|K1TZK1|K1TZK1_9ZZZZ Aspartate--ammonia ligase (Fragment) OS=human gut metagenome GN=LEA_03426 PE=4 SV=1
Probab=99.84  E-value=4e-25  Score=157.96  Aligned_cols=57  Identities=63%  Similarity=1.049  Sum_probs=52.8  Template_Neff=1.000

Q d12asa_         197 GIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQL  262 (327)
Q Consensus       197 gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql  262 (327)
                      -|||.|..|-|||-|+||||||.         +|-||+.|...|.-|.||||||||||+.....||
T Consensus         2 qiggqlksgirhdgrapdyddwt---------lncdilfwhkalgcalelssmgirvdpaamtrql   58 (58)
T tr|K1TZK1|K1TZ    2 QIGGQLKSGIRHDGRAPDYDDWT---------LNCDILFWHKALGCALELSSMGIRVDPAAMTRQL   58 (58)
Confidence            48999999999999999999995         7999999999999999999999999998776664


No 28
>tr|D8FGH8|D8FGH8_9FIRM Putative Aspartate--ammonia ligase (Fragment) OS=Peptoniphilus sp. oral taxon 836 str. F0141 GN=HMPREF9131_1461 PE=4 SV=1
Probab=99.84  E-value=4.3e-25  Score=165.94  Aligned_cols=69  Identities=42%  Similarity=0.697  Sum_probs=62.7  Template_Neff=1.382

Q d12asa_           8 QISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTL   78 (327)
Q Consensus         8 aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL   78 (327)
                      .|.-+|..|+|.|...|.|-.|.||..-|-..|-.||||| |++|..+.+.. |...|+|||||||||..|
T Consensus         3 ai~~lkd~f~r~l~~~lnl~rvsaplf~r~~~glnd~lsg-e~~v~f~~~~~-~~nveiv~slakwkrnal   71 (73)
T tr|D8FGH8|D8FG    3 AIKSLKDYFQRDLRSSLNLTRVSAPLFVRPSSGLNDNLSG-EKPVNFNPRNY-NINIEIIQSLAKWKRNAL   71 (73)
Confidence            4666899999999999999999999999999999999998 78999988875 578999999999999765


No 29
>tr|K1Y4G7|K1Y4G7_9BACT Aspartate-ammonia ligase (Fragment) OS=uncultured bacterium (gcode 4) GN=ACD_77C00154G0001 PE=4 SV=1
Probab=99.84  E-value=5.3e-25  Score=167.80  Aligned_cols=73  Identities=41%  Similarity=0.677  Sum_probs=68.2  Template_Neff=1.510

Q d12asa_           3 IAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR   75 (327)
Q Consensus         3 ~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR   75 (327)
                      |.....|..||..|+..|...|.|+.|.+||+-.-|.|..|||.|.|+.|....|+..|-.--||||||||||
T Consensus         6 ~~te~~i~~vk~~f~~~l~~~l~l~ris~pi~~~~g~gind~lng~er~v~f~~~~~q~~ravvvhslakwkr   78 (78)
T tr|K1Y4G7|K1Y4    6 LETEISIALVKETFSTELSRQLSLARISSPIAILDGTGINDDLNGCERPVAFPLKAMQDRRAVVVHSLAKWKR   78 (78)
Confidence            3445567889999999999999999999999999999999999999999999999999999999999999997


No 30
>tr|S7W4K5|S7W4K5_SPRLO Uncharacterized protein OS=Spraguea lophii (strain 42_110) GN=SLOPH_146 PE=4 SV=1
Probab=99.83  E-value=1.4e-24  Score=177.55  Aligned_cols=97  Identities=40%  Similarity=0.639  Sum_probs=82.9  Template_Neff=1.000

Q d12asa_           6 QRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSA   85 (327)
Q Consensus         6 q~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~   85 (327)
                      ..||.||||.|.....|.|.|-.|.||++-.-..|-.|.|.|.|||+......  .....+|||||||||..|....- .
T Consensus        26 mqqinfvkstfakkfsenlqltkvpaplfvtpasglndtltgiekaisfttag--kidctivhslakwkrealkryeg-d  102 (130)
T tr|S7W4K5|S7W4   26 MQQINFVKSTFAKKFSENLQLTKVPAPLFVTPASGLNDTLTGIEKAISFTTAG--KIDCTIVHSLAKWKREALKRYEG-D  102 (130)
Confidence            36899999999999999999999999999999999999999999998765443  34567899999999999876542 2


Q d12asa_          86 GEGLYTHMKALRPDEDRLSP  105 (327)
Q Consensus        86 geGiyTdMnAIRrDE~~ld~  105 (327)
                      -.|.||.|-|.|.|||.-.-
T Consensus       103 idgvytdmlairkdedpeti  122 (130)
T tr|S7W4K5|S7W4  103 IDGVYTDMLAIRKDEDPETI  122 (130)
Confidence            36899999999999986543


No 31
>tr|A0A0H2REP6|A0A0H2REP6_MORMO Asparagine synthetase AsnA (Fragment) OS=Morganella morganii GN=ABN09_02135 PE=4 SV=1
Probab=99.80  E-value=1.6e-23  Score=153.92  Aligned_cols=66  Identities=65%  Similarity=1.103  Sum_probs=63.1  Template_Neff=1.000

Q d12asa_         261 QLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSL  326 (327)
Q Consensus       261 Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~~~~~~  326 (327)
                      ||.|||..-|.||.||..|..||||||||||||||||.|||||-.||||||..||...|||.+..+
T Consensus         1 qldltgqqvrmeldwhkslaagempqtigggigqsrlvmlllqrshigqvqcsvwapevretitdm   66 (67)
T tr|A0A0H2REP6|    1 QLDLTGQQVRMELDWHKSLAAGEMPQTIGGGIGQSRLVMLLLQRSHIGQVQCSVWAPEVRETITDM   66 (67)
Confidence            678999999999999999999999999999999999999999999999999999999999987654


No 32
>tr|A4GX00|A4GX00_9LACT AsnA (Fragment) OS=Lactococcus garvieae GN=asnA PE=4 SV=1
Probab=99.78  E-value=9.4e-23  Score=141.47  Aligned_cols=36  Identities=61%  Similarity=1.055  Sum_probs=33.9  Template_Neff=1.000

Q d12asa_         283 EMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAA  318 (327)
Q Consensus       283 ~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~  318 (327)
                      |.|-||||||||||+.||||.--|||.||...|...
T Consensus         2 elpltigggigqsrmcmlllskvhigevqvslwdee   37 (49)
T tr|A4GX00|A4GX    2 ELPLTIGGGIGQSRMCMLLLSKVHIGEVQVSLWDEE   37 (49)
Confidence            789999999999999999999999999999999754


No 33
>tr|A0A090P1X0|A0A090P1X0_9VIBR Aspartate-ammonia ligase OS=Vibrio ponticus GN=JCM19238_73 PE=4 SV=1
Probab=99.53  E-value=9.2e-18  Score=126.23  Aligned_cols=55  Identities=33%  Similarity=0.457  Sum_probs=53.0  Template_Neff=1.464

Q d12asa_           2 YIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV   56 (327)
Q Consensus         2 ~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v   56 (327)
                      +.+||.||+.+|+.|+..|+.+|+|.|||||||+.+.+||||||||+||||-|++
T Consensus        15 ~~dt~~ai~~lk~~fe~~la~~l~l~rvsaplf~~~~~glndnlng~erpv~~~~   69 (69)
T tr|A0A090P1X0|   15 VYDTQKAIKKLKRYFEDQLAYKLYLFRVSAPLFVKPESGLNDNLNGYERPVAFIC   69 (69)
Confidence            5789999999999999999999999999999999999999999999999999874


No 34
>tr|B7XNE7|B7XNE7_ENTBH Aspartate-ammonia ligase (Fragment) OS=Enterocytozoon bieneusi (strain H348) GN=EBI_25665 PE=4 SV=1
Probab=99.33  E-value=2e-15  Score=108.63  Aligned_cols=54  Identities=35%  Similarity=0.565  Sum_probs=50.1  Template_Neff=1.000

Q d12asa_          39 DGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKA   95 (327)
Q Consensus        39 sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnA   95 (327)
                      +|+|||||| |.||.|..++  +..+.||||||||||+.|.+..-..-.||++||+.
T Consensus         2 sgvnddlng-eepvrfqtsd--gewcsvvhslakwkrwmlwklrddgvsgiwcdmrg   55 (56)
T tr|B7XNE7|B7XN    2 SGVNDDLNG-EEPVRFQTSD--GEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRG   55 (56)
Confidence            699999999 6699999886  79999999999999999999999889999999963


No 35
>tr|K1TB71|K1TB71_9ZZZZ Aspartate-ammonia ligase (Fragment) OS=human gut metagenome GN=OBE_03908 PE=4 SV=1
Probab=99.09  E-value=2.4e-13  Score=118.09  Aligned_cols=59  Identities=31%  Similarity=0.491  Sum_probs=56.4  Template_Neff=1.000

Q d12asa_           2 YIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP   60 (327)
Q Consensus         2 ~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~   60 (327)
                      +.+||+||..+|..|+.+|+.+|||.||+||||+...+||||||||.||||-|++...+
T Consensus        94 vydtqkaigllkrlfedqlaaklnlfrvsaplfleeasglndnlngyerpvlfdipqag  152 (176)
T tr|K1TB71|K1TB   94 VYDTQKAIGLLKRLFEDQLAAKLNLFRVSAPLFLEEASGLNDNLNGYERPVLFDIPQAG  152 (176)
Confidence            56899999999999999999999999999999999999999999999999999998753


No 36
>tr|W1VJV4|W1VJV4_STRPA Aspartate-ammonia ligase (Fragment) OS=Streptococcus parasanguinis DORA_23_24 GN=Q616_SPPC00677G0003 PE=4 SV=1
Probab=98.59  E-value=2.1e-10  Score=84.52  Aligned_cols=54  Identities=65%  Similarity=1.019  Sum_probs=51.5  Template_Neff=1.349

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQV   54 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F   54 (327)
                      +||++||||||||+||+|+||+|||+||||+|||||||||+||||||.|..|.|
T Consensus         4 ~~i~~q~~isfvk~~f~~~l~~~l~~ievq~pils~vgdg~qdnlsgiehpvsv   57 (57)
T tr|W1VJV4|W1VJ    4 SFIHQQREISFVKSTFTQYLKDKLGLIEVQGPILSKVGDGMQDNLSGIEHPVSV   57 (57)
Confidence            699999999999999999999999999999999999999999999999987653


No 37
>tr|A0A090PG64|A0A090PG64_9VIBR Aspartate-ammonia ligase OS=Vibrio ponticus GN=JCM19238_72 PE=4 SV=1
Probab=98.49  E-value=6.1e-10  Score=74.94  Aligned_cols=28  Identities=57%  Similarity=0.850  Sum_probs=25.3  Template_Neff=1.000

Q d12asa_          60 PDAQFEVVHSLAKWKRQTLGQHDFSAGE   87 (327)
Q Consensus        60 ~~~~~EIVhSLAKWKR~aL~ky~~~~ge   87 (327)
                      ++...|||||||||||+||.+|||-..+
T Consensus         5 ngahvevvhslakwkrmalgrygflede   32 (37)
T tr|A0A090PG64|    5 NGAHVEVVHSLAKWKRMALGRYGFLEDE   32 (37)
Confidence            4789999999999999999999997655


No 38
>tr|K0NLJ8|K0NLJ8_9LACO Uncharacterized protein OS=Lactobacillus equicursoris DSM 19284 = JCM 14600 = CIP 110162 GN=BN147_07820 PE=4 SV=1
Probab=98.46  E-value=8.3e-10  Score=72.56  Aligned_cols=28  Identities=50%  Similarity=0.873  Sum_probs=26.4  Template_Neff=1.000

Q d12asa_         299 MLLLQLPHIGQVQAGVWPAAVRESVPSL  326 (327)
Q Consensus       299 M~lL~k~HIgEVq~svW~~~~~~~~~~~  326 (327)
                      |+||.++||||||+++||+++++.|.+.
T Consensus         1 mlllgrahigevqagiwpddmieqcaea   28 (33)
T tr|K0NLJ8|K0NL    1 MLLLGRAHIGEVQAGIWPDDMIEQCAEA   28 (33)
Confidence            8999999999999999999999999863


No 39
>tr|A0A084U2P1|A0A084U2P1_MYCIO Uncharacterized protein OS=Mycoplasma iowae DK-CPA GN=P271_52 PE=4 SV=1
Probab=98.21  E-value=8e-09  Score=99.84  Aligned_cols=291  Identities=15%  Similarity=0.225  Sum_probs=247.8  Template_Neff=1.677

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      |++.|..-|+...+....-.....+|+...-|+...-.+.+|=+  .+-|||.|+--.. -.-.||+-|-.---|+...-
T Consensus        11 sf~~~v~~i~~~~~~~~~~~~~~y~~i~ldlp~~~n~~~~~n~~--~~~r~i~fd~~~d-~~iyeii~~~dn~iry~~~~   87 (325)
T tr|A0A084U2P1|   11 SFKDVVKFIDYFYNQLAETIKQKYNLINLDLPLVSNMKSDVNLL--NNNRAINFDNYND-KNIYEIIYEPDNMIRYYCWF   87 (325)
Confidence            57788888999888888888899999999999999999999876  6689999997765 67789999999999999999


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEA--AVSEEFGLAPFL  158 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~--~v~~~y~l~~~L  158 (327)
                      |......-+++-...|.||- .+.++.|+----..+|--+.+++|+-++|++.++-||+.+.+.-.  .+.+.|.+   -
T Consensus        88 ~e~~~~~~~~sky~~i~rd~-~~n~~~~~e~~~~nfe~~~~e~~~~~~~~~~~~~~~w~~~~~~~~~s~~~~~~~~---~  163 (325)
T tr|A0A084U2P1|   88 LELTNNDVVVSKYKQINRDA-IINNSSSIENNMLNFEFFILEEQKKEEYVLDLINYFWNIFLKIVKSSSLNKNYRL---E  163 (325)
Confidence            99999999999999999998 899999988888899999999999999999999999999876532  23333433   2


Q d12asa_         159 PDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHD-----VRAPDYDDWSTPSELGHAGLNGDI  233 (327)
Q Consensus       159 p~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd-----~RapDYDDW~t~~~~~~~gLNGDi  233 (327)
                      .+.+.-++..|+..+|--++-+.-..++.+-.|...+-.|..+.    .||     ..+-|..|.+         ----+
T Consensus       164 ~~~~~c~s~~~~~~~~l~~p~~~~~~~~~~~~gi~~~k~~~~kf----~h~~~~~~~~~~~s~d~~---------nt~sl  230 (325)
T tr|A0A084U2P1|  164 TKKIRCVSLKEIKKMYLVLPIKDAVDKFILNNGIHLIKDISNKF----EHDSNVYLEKSSDSHDFE---------NTYSL  230 (325)
Confidence            35677788899999999999999999999999877655555444    565     3456666666         44468


Q d12asa_         234 LVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       234 lv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ++|+...+..=||--...|.+=+++++|..+-+. .-....|-.-+.+|.---|.+--|--..|+.+||.|..|.|+-+
T Consensus       231 ~~~~~~~~~~~elit~~~rpnw~~~~kq~~i~~~-~~~~n~f~di~~~~~~~~~~s~ki~~d~l~~y~lsk~~~~eips  308 (325)
T tr|A0A084U2P1|  231 LFFDENSQQVKELITITFRPNWDTYKKQKGINGE-KILNNNFTDILKKDSEVNTCSFKINFDLLIYYFLSKTDIQEIPS  308 (325)
Confidence            8999999999999999999999999999997666 45778888889999888888889999999999999999999865


No 40
>tr|B7XQ79|B7XQ79_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_27632 PE=4 SV=1
Probab=98.21  E-value=8.6e-09  Score=76.27  Aligned_cols=48  Identities=35%  Similarity=0.578  Sum_probs=47.0  Template_Neff=1.000

Q d12asa_         106 LHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG  153 (327)
Q Consensus       106 ~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~  153 (327)
                      +||+.|.|+|=.+||..+.||++|||+||+.||..|+.+.+.+...||
T Consensus         1 mhsyqveqfdgkkvipsekrnieflketvreiyrglkaakervdaeyp   48 (59)
T tr|B7XQ79|B7XQ    1 MHSYQVEQFDGKKVIPSEKRNIEFLKETVREIYRGLKAAKERVDAEYP   48 (59)
Confidence            599999999999999999999999999999999999999999999998


No 41
>tr|A0A0G1B9N3|A0A0G1B9N3_9BACT Aspartyl-tRNA synthetase (Fragment) OS=Microgenomates (Daviesbacteria) bacterium GW2011_GWA2_42_7 GN=UV41_C0041G0001 PE=3 SV=1
Probab=97.71  E-value=2.7e-07  Score=94.08  Aligned_cols=255  Identities=20%  Similarity=0.314  Sum_probs=179.3  Template_Neff=5.200

Q d12asa_          23 RLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR  102 (327)
Q Consensus        23 ~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~  102 (327)
                      .-|.++|.+|.+      |+-+-.|....  |.+.-- +..+-..||..=.|-+|+..  |   +-+|+=--.+|..- .
T Consensus        53 ~~gF~ei~~P~l------l~~~~Egg~~~--f~v~Yf-g~~a~LaQS~qlykQ~~i~~--~---~kVF~i~P~fRaE~-s  117 (344)
T tr|A0A0G1B9N3|   53 DRGFVEIHPPKL------LWPASEGGATL--FEVDYF-GKKAYLAQSPQLYKQMAIAA--F---EKVFEIGPNFRAEK-S  117 (344)
Confidence            348999999999      44444656666  555444 67899999999999999987  3   45777777888776 6


Q d12asa_         103 LSPLHSVYVDQWDWERVMGDG-ERQFSTLKSTVEAIWAGIKAT-EAAV---SEEFG-LAPFLPDQIHFVHSQELLSR---  173 (327)
Q Consensus       103 ld~~HSiyVDQWDWEkvI~~~-dRnl~~Lk~tV~kIy~al~~t-e~~v---~~~y~-l~~~Lp~~I~FI~sqeL~~~---  173 (327)
                      -++-|=.-.-|+|.|.-+... +.-++.+.++|..|++.+++- ...+   ...+| .+..+ ..++|=.+-++++.   
T Consensus       118 ~T~RHL~EFt~lD~Em~~~~~~~dvl~~iE~li~~i~~~l~~~~~~~l~~~~~~~p~~~~pf-~R~~y~eai~~L~~~~~  196 (344)
T tr|A0A0G1B9N3|  118 NTHRHLTEFTQLDFEMAFADHYDDVMDLIEELIKYIFKRLQEHCKEELETLERQLPPPKKPF-PRITYDEAIEMLREAGV  196 (344)
Confidence            789999999999999999886 889999999999998887763 2222   22233 22111 23333222222222   


Q d12asa_         174 -----YPDLDAKGRERAIAK------DLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLED  242 (327)
Q Consensus       174 -----YP~LtpkeRE~~i~k------e~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~  242 (327)
                           ..|+++. -|+.+++      -...+||+.+-..+  -.--+.+.|+....+         .|=|++.      .
T Consensus       197 ~~~~~~~D~~~~-~E~~l~~~~~~~~~~~p~~i~~~P~~~--r~FY~~~~~~~~~~~---------~s~DL~~------g  258 (344)
T tr|A0A0G1B9N3|  197 EKIGEEEDLGTE-HERLLGKLVKEKYGTDPFFITNYPAAI--RPFYMMPDPDDPGYS---------NSYDLIM------G  258 (344)
Confidence                 2355443 3455555      34556666655333  344566667655566         7778776      6


Q d12asa_         243 AFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       243 a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      --|++|-|-|+ |.+.|.++.+..+-. .....|+-.+..--+|.+-|+|||-=|+.||||...||+||+.
T Consensus       259 ~gEi~sGsqR~~d~e~L~~r~~~~gid-~~~~~~Yld~~k~g~~PhaG~GiGlERlv~~~~gl~~Ire~~l  328 (344)
T tr|A0A0G1B9N3|  259 YGEILSGAQRIHDYEELLERMKEKGID-PESFKWYLDLFKYGCPPHAGFGIGLERLVMWLLGLKHIREVTL  328 (344)
Confidence            78999999998 578888888888884 4444444444445799999999999999999999999999874


No 42
>tr|X1IS21|X1IS21_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S03H2_55971 PE=4 SV=1
Probab=97.49  E-value=1.5e-06  Score=59.09  Aligned_cols=28  Identities=36%  Similarity=0.713  Sum_probs=26.5  Template_Neff=1.285

Q d12asa_         298 TMLLLQLPHIGQVQAGVWPAAVRESVPS  325 (327)
Q Consensus       298 ~M~lL~k~HIgEVq~svW~~~~~~~~~~  325 (327)
                      +|+||+|+|||||.+.+||.-..+.|++
T Consensus         1 fm~ll~~ah~gevs~tvwp~ilk~mc~k   28 (35)
T tr|X1IS21|X1IS    1 FMYLLRKAHIGEVSVTVWPQILKDMCRK   28 (35)
Confidence            6999999999999999999999999975


No 43
>tr|R5PEQ4|R5PEQ4_9BACT Aspartate--ammonia ligase OS=Prevotella sp. CAG:487 GN=BN679_01862 PE=4 SV=1
Probab=97.16  E-value=1e-05  Score=60.22  Aligned_cols=36  Identities=25%  Similarity=0.458  Sum_probs=33.1  Template_Neff=1.000

Q d12asa_           2 YIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV   37 (327)
Q Consensus         2 ~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~   37 (327)
                      .+.|+++|+.||++|+++|+.+|.|.||+||||...
T Consensus        16 mkrteqgiklikdffqenlstelrlsrvtaplf~~~   51 (56)
T tr|R5PEQ4|R5PE   16 MKRTEQGIKLIKDFFQENLSTELRLSRVTAPLFXXX   51 (56)
Confidence            357899999999999999999999999999999754


No 44
>tr|B7XNP0|B7XNP0_ENTBH Uncharacterized protein OS=Enterocytozoon bieneusi (strain H348) GN=EBI_25947 PE=4 SV=1
Probab=96.89  E-value=3.8e-05  Score=59.47  Aligned_cols=47  Identities=36%  Similarity=0.522  Sum_probs=36.1  Template_Neff=1.000

Q d12asa_          39 DGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALR   97 (327)
Q Consensus        39 sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIR   97 (327)
                      +|+|||||| |.|-+-+-.   -...-|||||||=|+-||.        -+|..||.|-
T Consensus         2 sgvnddlng-espsdsdhp---revvhvvhslakgkkdale--------tlyegmngig   48 (68)
T tr|B7XNP0|B7XN    2 SGVNDDLNG-ESPSDSDHP---REVVHVVHSLAKGKKDALE--------TLYEGMNGIG   48 (68)
Confidence            699999999 777554432   5678899999999999985        4666777664


No 45
>tr|K2C8Q9|K2C8Q9_9BACT Uncharacterized protein (Fragment) OS=uncultured bacterium (gcode 4) GN=ACD_39C01328G0001 PE=4 SV=1
Probab=96.69  E-value=9e-05  Score=47.89  Aligned_cols=21  Identities=43%  Similarity=0.729  Sum_probs=18.3  Template_Neff=1.000

Q d12asa_         305 PHIGQVQAGVWPAAVRESVPS  325 (327)
Q Consensus       305 ~HIgEVq~svW~~~~~~~~~~  325 (327)
                      .||||||+|+||.+..-.|+.
T Consensus         1 rhigevqvsvwpadevvrcea   21 (27)
T tr|K2C8Q9|K2C8    1 RHIGEVQVSVWPADEVVRCEA   21 (27)
Confidence            499999999999998877764


No 46
>tr|R1C3V4|R1C3V4_EMIHU Lysyl-tRNA synthetase (Fragment) OS=Emiliania huxleyi GN=EMIHUDRAFT_244259 PE=3 SV=1
Probab=95.29  E-value=0.005  Score=46.16  Aligned_cols=29  Identities=41%  Similarity=0.584  Sum_probs=27.3  Template_Neff=1.000

Q d12asa_         283 EMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       283 ~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      -||.|.|=|.|..||||+|-..+.|.||-
T Consensus         5 glpptagwgmgidrlcmlltdnesikevl   33 (52)
T tr|R1C3V4|R1C3    5 GLPPTAGWGMGIDRLCMLLTDNESIKEVL   33 (52)
Confidence            48999999999999999999999999985


No 47
>tr|A0A079G0B3|A0A079G0B3_ECOLX tRNA synthetases class II family protein OS=Escherichia coli 2-460-02_S3_C3 GN=AC49_1465 PE=3 SV=1
Probab=95.15  E-value=0.0056  Score=47.01  Aligned_cols=38  Identities=45%  Similarity=0.538  Sum_probs=32.7  Template_Neff=2.717

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .+|-.++.-| +|.|-|=|||..||.|+|-....|.||-
T Consensus         5 eDf~~ALEYG-mPPtgG~GiGIDRLvMllT~~~sIrdVi   42 (51)
T tr|A0A079G0B3|    5 EDFVTALEYG-MPPTGGLGIGIDRLVMLLTNSATIRDVI   42 (51)
Confidence            3566666666 9999999999999999999999999985


No 48
>tr|U1F973|U1F973_9ACTN Lysyl-tRNA ligase (Fragment) OS=Propionibacterium granulosum DSM 20700 GN=lysS PE=4 SV=1
Probab=95.14  E-value=0.0056  Score=49.71  Aligned_cols=40  Identities=45%  Similarity=0.662  Sum_probs=34.8  Template_Neff=2.908

Q d12asa_         270 RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       270 r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....+|-.++.-| +|.|-|-|||..||+|+|-. .-|.+|-
T Consensus        19 ~~Dedfv~ALeyG-mPPt~G~GiGiDRLvmllt~-~SIRdvi   58 (68)
T tr|U1F973|U1F9   19 MMDEDFVEALEYG-MPPTAGFGIGIDRLVMLLTD-PSIRDVI   58 (68)
Confidence            4566788888877 99999999999999999999 8888874


No 49
>tr|F3VDV0|F3VDV0_SHIDY Lysyl-tRNA synthetase domain protein OS=Shigella dysenteriae 155-74 GN=SD15574_4901 PE=3 SV=1
Probab=95.13  E-value=0.0062  Score=43.48  Aligned_cols=28  Identities=54%  Similarity=0.697  Sum_probs=25.5  Template_Neff=1.893

Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +|.|-|-|||..|+.|+|-++.-|..|-
T Consensus         1 mPPtgG~GiGIDR~iMLlT~~~sIRDvl   28 (37)
T tr|F3VDV0|F3VD    1 MPPTGGLGIGIDRLIMLLTNSPSIRDVL   28 (37)
Confidence            6999999999999999999999888763


No 50
>tr|T0YGT2|T0YGT2_9ZZZZ Lysyl-tRNA synthetase (Fragment) OS=mine drainage metagenome GN=B1B_17746 PE=4 SV=1
Probab=95.08  E-value=0.0055  Score=57.92  Aligned_cols=55  Identities=38%  Similarity=0.598  Sum_probs=46.6  Template_Neff=4.139

Q d12asa_         255 ADTLKHQLALTGDED--RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~~~~--r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+.+|++ .+..+  -...+|-.++.-| +|.|-|=|||..||+|+|....-|.+|-
T Consensus       107 r~RFe~Q~~-~Gd~EA~~~Dedfi~ALEyG-mPPTgG~GiGIDRlvMlltn~~sIrdVl  163 (165)
T tr|T0YGT2|T0YG  107 RQRFEEQAK-AGDDEAQEIDEDFIEALEYG-MPPTGGLGIGIDRLVMLLTNSPSIRDVL  163 (165)
Confidence            456788888 44444  3778899999999 9999999999999999999999999873


No 51
>tr|D9PH19|D9PH19_9ZZZZ Lysyl-tRNA synthetase OS=sediment metagenome GN=lysS PE=3 SV=1
Probab=94.87  E-value=0.0078  Score=63.59  Aligned_cols=274  Identities=20%  Similarity=0.263  Sum_probs=151.3  Template_Neff=4.100

Q d12asa_           7 RQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP-DAQFEVVHSLAKWKRQTLGQHDFSA   85 (327)
Q Consensus         7 ~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~-~~~~EIVhSLAKWKR~aL~ky~~~~   85 (327)
                      +.|+.|+++|.+     .+.++|..|+.....-|=      .-||  |...... +...-.--|..=    .||+.=...
T Consensus        52 ~ii~~iR~ff~~-----~gFlEVeTP~L~~i~GGA------~ArP--F~Th~nald~dlyLRIAPEL----yLKrLiVGG  114 (375)
T tr|D9PH19|D9PH   52 KIISSIRRFFDN-----RGFLEVETPMLQPIPGGA------AARP--FVTHHNALDMDLYLRIAPEL----YLKRLVVGG  114 (375)
Confidence            445666666654     577888888776554221      1233  4432210 111111111111    133332322


Q d12asa_          86 GEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQIHF  164 (327)
Q Consensus        86 geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~F  164 (327)
                      -+-+|-=-..- |+| -.|..|+=-+--+-|=..-.....-++..++.+..+..++..+...   .|+ -.-.+......
T Consensus       115 ~ekVfEI~r~F-RNE-Gid~~HnPEFtmlE~Y~AYady~d~m~~te~Li~~~a~~~~g~~~~---~~~~~~id~~~p~~r  189 (375)
T tr|D9PH19|D9PH  115 FERVFEIGRNF-RNE-GIDTTHNPEFTMLEFYQAYADYNDLMDLTEELISHLAQEVLGTTKV---PYGGPEIDFSPPWRR  189 (375)
Confidence            23344333334 455 5899999888888888877777777777777777776666544322   133 23345666666


Q d12asa_         165 VHSQELLSRYPDL------DAKGRERAIAKDLGAVF--LVGIGGKLS--------DGHRHDVRAPDYDDWSTPSELGHAG  228 (327)
Q Consensus       165 I~sqeL~~~YP~L------tpkeRE~~i~ke~gAvF--i~gIG~~L~--------~G~~Hd~RapDYDDW~t~~~~~~~g  228 (327)
                      |+..++.+.| ++      +..+...++|++.+.-+  -.+.|..+.        ....+-..--||..-.+|-..    
T Consensus       190 i~~~e~i~~~-g~d~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~l~d~l~~~~ve~~l~~PTFI~d~P~~~SPLAK----  264 (375)
T tr|D9PH19|D9PH  190 ISMVEALKKY-GIDIPPDLDDPEELKAIAEKIGVKVDKPWTAGRLLDELFEHFVEPKLIQPTFITDYPVEMSPLAK----  264 (375)
Confidence            7777777776 43      33445555666444100  111111111        111111111222222211111    


Q d12asa_         229 LNGDILVWNPVLEDAFELSSMGIRVD------------ADTLKHQLALTG-D--ED--RLELEWHQALLRGEMPQTIGGG  291 (327)
Q Consensus       229 LNGDilv~n~~l~~a~ElSSmGirVd------------~~~L~~Ql~~~~-~--~~--r~~~~~h~~ll~~~lP~TIgGG  291 (327)
                      .+.+    ||.+-..||+==.|+-+-            ++.+..|++... .  .+  ....+|-+++.-| +|.|-|=|
T Consensus       265 ~~~~----~p~~~eRFElfi~G~Ei~NaysELNDP~~Qr~RF~~Q~~~~~~~gd~Ea~~~DedFi~ALEyG-mPPt~G~G  339 (375)
T tr|D9PH19|D9PH  265 RHRS----NPGLTERFELFIGGKEIANAYSELNDPVDQRERFEEQAKDKEVAGDDEAMPIDEDFIEALEYG-MPPTGGIG  339 (375)
Confidence            1111    344555666666666542            356777777755 2  22  3778899999999 99999999


Q d12asa_         292 IGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       292 IgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ||-.||+|+|....+|.||..
T Consensus       340 iGIDRLvMlLT~~~~Irevil  360 (375)
T tr|D9PH19|D9PH  340 IGIDRLVMLLTNSHSIREVIL  360 (375)
Confidence            999999999999999999964


No 52
>tr|W7QHE3|W7QHE3_9GAMM Uncharacterized protein OS=Halomonas sp. BC04 GN=Q427_17735 PE=4 SV=1
Probab=94.37  E-value=0.017  Score=52.99  Aligned_cols=56  Identities=36%  Similarity=0.464  Sum_probs=46.3  Template_Neff=3.300

Q d12asa_         255 ADTLKHQLALTG----DEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~----~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+.+|.+.+.    .+.....+|-.+|.-| +|.|-|=|||-.||+|+|..-+-|.||-
T Consensus        67 r~RF~~q~~~r~~~Gd~e~~iDedfl~ALE~G-mPPt~G~giGiDRLvmllT~~~sIrevl  126 (133)
T tr|W7QHE3|W7QH   67 RERFEEQQRARAAAGDDEMPIDEDFLEALEYG-MPPTGGLGIGIDRLVMLLTGAQSIREVI  126 (133)
Confidence            355677777665    3345777888898888 9999999999999999999999999984


No 53
>tr|X1CIB2|X1CIB2_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S01H4_36895 PE=4 SV=1
Probab=94.31  E-value=0.022  Score=50.95  Aligned_cols=56  Identities=36%  Similarity=0.524  Sum_probs=43.4  Template_Neff=1.000

Q d12asa_         255 ADTLKHQLALTGDED-----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~~~~-----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+..||+.+...+     ....+|-+++.-| +|.|-|=|||-.||.|||-..+.|.||-
T Consensus         8 rerfeaqlklaergddeatasidedflraleyg-mpptsgmgigmdrlimfltnnqsiqevl   68 (144)
T tr|X1CIB2|X1CI    8 RERFEAQLKLAERGDDEATASIDEDFLRALEYG-MPPTSGMGIGMDRLIMFLTNNQSIQEVL   68 (144)
Confidence            345566776554322     2456777887777 8999999999999999999999999985


No 54
>tr|A0A0G9GN21|A0A0G9GN21_LACPN Lysyl-tRNA synthetase (Fragment) OS=Lactobacillus plantarum GN=WP50_37500 PE=3 SV=1
Probab=94.25  E-value=0.018  Score=50.69  Aligned_cols=56  Identities=34%  Similarity=0.442  Sum_probs=43.1  Template_Neff=3.979

Q d12asa_         256 DTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       256 ~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +.+..|...+..-+    ....+|-.++.-| +|.|-|=|||-.||+|+|..-.-|.+|-.
T Consensus        34 ~Rf~~q~~~r~~g~~ea~~~D~d~l~ALeyG-mPp~~G~giGiDRLvmllt~~~sIrdVil   93 (103)
T tr|A0A0G9GN21|   34 ERFEAQERMREAGDEEAMPMDEDFIEALEYG-MPPTGGLGIGIDRLVMLLTGADSIRDVIL   93 (103)
Confidence            44556655553222    3566777888777 99999999999999999999999999853


No 55
>tr|G4NE53|G4NE53_MAGO7 Lysine--tRNA ligase OS=Magnaporthe oryzae (strain 70-15 / ATCC MYA-4617 / FGSC 8958) GN=MGG_00161 PE=3 SV=1
Probab=93.99  E-value=0.025  Score=63.39  Aligned_cols=42  Identities=36%  Similarity=0.524  Sum_probs=37.4  Template_Neff=3.500

Q d12asa_         270 RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       270 r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      .....|-.++.-| +|.|-|=|||.-||||+|+.+.-|.||-.
T Consensus       571 ~iD~~Yv~Amk~G-MPPtGG~G~GIDRL~MLLt~~~rI~dVL~  612 (638)
T tr|G4NE53|G4NE  571 SIDEQYVEAMKYG-MPPTGGWGLGIDRLCMLLTGKKRIEDVLS  612 (638)
Confidence            3667888888777 99999999999999999999999999853


No 56
>tr|Q48986|Q48986_MYCCA Lys-tRNA synthetase (Fragment) OS=Mycoplasma capricolum PE=3 SV=1
Probab=93.62  E-value=0.046  Score=42.30  Aligned_cols=32  Identities=47%  Similarity=0.604  Sum_probs=27.9  Template_Neff=1.000

Q d12asa_         280 LRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       280 l~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++..+|.|.|-|||..||.|+|-....|..|-
T Consensus        18 lehampptagigigidrlvmlltncdsikdvl   49 (58)
T tr|Q48986|Q489   18 LEHAMPPTAGIGIGIDRLVMLLTNCDSIKDVL   49 (58)
Confidence            34568999999999999999999998888763


No 57
>tr|S2EWF1|S2EWF1_9ARCH Uncharacterized protein OS=Candidatus Nitrosoarchaeum limnia BG20 GN=BG20_I1147 PE=3 SV=1
Probab=93.57  E-value=0.048  Score=45.17  Aligned_cols=68  Identities=28%  Similarity=0.399  Sum_probs=52.7  Template_Neff=1.143

Q d12asa_         243 AFELSSMGIRVDADT-LKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       243 a~ElSSmGirVd~~~-L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .+||||-.-|+.+.+ |....+..+. .-....||--+-.--.|.-.|.|||.-||.|.|-.-+.|..|-
T Consensus         4 dlelssgstr~~k~~ele~rm~nkgm-kid~f~yhlg~fdygvpphagcgiglerl~maltg~enird~t   72 (83)
T tr|S2EWF1|S2EW    4 DLELSSGSTRMRKRHELEDRMNNKGM-KIDTFEYHLGVFDYGVPPHAGCGIGLERLMMALTGIENIRDTT   72 (83)
Confidence            368888888876643 4334433333 2356789999999999999999999999999999999998764


No 58
>tr|A0A0C2D2D0|A0A0C2D2D0_9BILA Uncharacterized protein OS=Ancylostoma duodenale GN=ANCDUO_13630 PE=3 SV=1
Probab=93.56  E-value=0.047  Score=43.48  Aligned_cols=36  Identities=42%  Similarity=0.552  Sum_probs=30.9  Template_Neff=1.306

Q d12asa_         275 WHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       275 ~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |..++.-| ||.|-|=|.|..||.|||-.|..|.||-
T Consensus         6 fctaleyg-lppt~gwgmgidrltm~l~d~nnik~v~   41 (65)
T tr|A0A0C2D2D0|    6 FCTALEYG-LPPTGGWGMGIDRLTMFLSDSNNIKDVI   41 (65)
Confidence            44455445 9999999999999999999999999984


No 59
>tr|A0A0G1WAG4|A0A0G1WAG4_9BACT Aspartyl-tRNA synthetase OS=Parcubacteria bacterium GW2011_GWC2_49_9 GN=UY52_C0014G0034 PE=3 SV=1
Probab=93.32  E-value=0.042  Score=60.07  Aligned_cols=257  Identities=19%  Similarity=0.287  Sum_probs=136.1  Template_Neff=4.900

Q d12asa_          17 SRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA--LPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMK   94 (327)
Q Consensus        17 ~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~--~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMn   94 (327)
                      -+-| .+.|.++|..|++... +|=      .-||  |.+..  .++...-.=+|..=|+..++- -||   +-.|- ..
T Consensus       152 R~fl-~~~gFlEVETPiL~~~-gga------~Arp--F~T~~n~l~g~~~~LrispELylKrLmV-gG~---ervfE-Ig  216 (463)
T tr|A0A0G1WAG4|  152 RRFL-DDRGFLEVETPMLQPI-GGA------GARP--FITHHNALDGMDLYLRIAPELYLKRLIV-GGF---ERVFE-IG  216 (463)
Confidence            3344 5679999999998877 542      2456  55554  323455455555555444331 233   22232 23


Q d12asa_          95 ALRPDEDRLSPLHSV-YVDQWDWERVMGDGERQFSTLKSTVEAIWAGIK-ATEAAVSEEFGLAPFLPDQIHFVHSQELLS  172 (327)
Q Consensus        95 AIRrDE~~ld~~HSi-yVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~-~te~~v~~~y~l~~~Lp~~I~FI~sqeL~~  172 (327)
                      ..=||| .+|..|.- -+-|+|++......+-.++...+.+..+...+. ... .+...|+          =|+..++.+
T Consensus       217 r~FRNE-gid~~hnP~EFT~lE~y~ay~d~~d~m~l~E~Li~~l~~~v~~~~~-~~~~pf~----------rit~~dai~  284 (463)
T tr|A0A0G1WAG4|  217 RNFRNE-GIDATHNPPEFTMLEFYMAYADYEDMMDLTEELIRHLAKEVFPEID-IFQKPFP----------RITMVEAIK  284 (463)
Confidence            333588 79999999 999999999999998888777777777777766 222 2222232          133344444


Q d12asa_         173 RYPDLDAKG-RERAIAKD-----L---GAVFLVG-IGGKLSDGHRHDVRAP---DYDDWSTPSELGHAGLNGDILVW--N  237 (327)
Q Consensus       173 ~YP~Ltpke-RE~~i~ke-----~---gAvFi~g-IG~~L~~G~~Hd~Rap---DYDDW~t~~~~~~~gLNGDilv~--n  237 (327)
                      .| ..+... |...+-..     .   ..+|++. ..+++.. .-|-.-+|   |-++.+     ..   ++++++.  +
T Consensus       285 ~~-g~D~~~~~~~~~~~~~l~~~~v~d~PtFv~~~~~~~~~~-~hp~~~sPla~~~~~~k-----~~---~~~~~~~rfe  354 (463)
T tr|A0A0G1WAG4|  285 KY-GIDKPDLREDKNDPKELAFCLIVDQPTFITDEEEGRWTA-SHPPEISPLAEDIELLK-----RN---RGKILTERFE  354 (463)
Confidence            44 222111 11100000     1   2333333 1111100 00000011   111111     00   2222221  1


Q d12asa_         238 PVLEDAFELSSMGIRV-DADTLKHQLALT----GDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       238 ~~l~~a~ElSSmGirV-d~~~L~~Ql~~~----~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      -++ .-.||+..-+|. |+....+-.+..    .+......+|-.++.-| +|.|-|=|||-.||+|+|.....|.||.+
T Consensus       355 lvi-~G~Ei~ng~selnDP~~Q~~~F~~~~~~d~Ea~~~d~dfl~AleyG-mPP~gG~GiGIDRLvMllt~~~sIRdVi~  432 (463)
T tr|A0A0G1WAG4|  355 LFI-NGREIANGFSELNDPIDQRERFEEKGKSDDEAMMMDEDFLEALEYG-MPPTGGLGIGIDRLVMLLTNSPSIRDVIL  432 (463)
Confidence            111 123444444443 333333333332    12223667777888888 99999999999999999999999999964


No 60
>tr|K1T7I3|K1T7I3_9ZZZZ Lysyl-tRNA synthetase (Fragment) OS=human gut metagenome GN=LEA_11318 PE=4 SV=1
Probab=93.30  E-value=0.056  Score=53.84  Aligned_cols=56  Identities=34%  Similarity=0.535  Sum_probs=45.4  Template_Neff=1.835

Q d12asa_         255 ADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+..|++.+..-+    -...+|-+++.=| +|.|-|=|||-.||.||+-....|.||-
T Consensus       118 ~eRF~~Qm~Ls~kGdDeam~IDqDFlrAL~yG-MPPTSG~GIG~DRL~MlmTgq~~IQEVl  177 (260)
T tr|K1T7I3|K1T7  118 EERFEEQMKLSEKGDDEAMFIDQDFLRALQYG-MPPTSGIGIGMDRLVMLMTGQQSIQEVL  177 (260)
Confidence            355677777665433    2567788888888 9999999999999999999999999984


No 61
>tr|A0A0G0L971|A0A0G0L971_9BACT Uncharacterized protein OS=candidate division TM6 bacterium GW2011_GWF2_38_10 GN=US69_C0008G0006 PE=4 SV=1
Probab=93.27  E-value=0.062  Score=42.02  Aligned_cols=28  Identities=39%  Similarity=0.637  Sum_probs=25.8  Template_Neff=1.000

Q d12asa_           3 IAKQRQISFVKSHFSRQLEERLGLIEVQ   30 (327)
Q Consensus         3 ~eTq~aI~~iK~~F~~~L~~~LnL~rVs   30 (327)
                      ..|+-+|.++|++|++.|+.+|||.||.
T Consensus        33 iktelgitfvkdtfqrllaqelnlmrvd   60 (61)
T tr|A0A0G0L971|   33 IKTELGITFVKDTFQRLLAQELNLMRVD   60 (61)
Confidence            4688899999999999999999999985


No 62
>tr|V2PR40|V2PR40_SALET Asparagine synthetase AsnA OS=Salmonella enterica subsp. enterica serovar Stanleyville str. CFSAN000624 GN=CFSAN000624_12912 PE=4 SV=1
Probab=93.26  E-value=0.063  Score=36.99  Aligned_cols=29  Identities=100%  Similarity=1.331  Sum_probs=28.2  Template_Neff=1.000

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEV   29 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rV   29 (327)
                      |||||||||||||||||||||||||||||
T Consensus         4 ayiakqrqisfvkshfsrqleerlgliev   32 (32)
T tr|V2PR40|V2PR    4 AYIAKQRQISFVKSHFSRQLEERLGLIEV   32 (32)
Confidence            79999999999999999999999999996


No 63
>tr|A0A0G0WIP2|A0A0G0WIP2_9BACT Aspartyl-tRNA synthetase OS=Microgenomates (Woesebacteria) bacterium GW2011_GWB1_41_10 GN=UU32_C0048G0005 PE=3 SV=1
Probab=92.84  E-value=0.076  Score=51.22  Aligned_cols=127  Identities=26%  Similarity=0.387  Sum_probs=89.8  Template_Neff=2.490

Q d12asa_         173 RYPDLDAKGRERAIA----KDLGAVFLVGIGGKLSDGHRHDVRAPDYDD--WSTPSELGHAGLNGDILVWNPVLEDAFEL  246 (327)
Q Consensus       173 ~YP~LtpkeRE~~i~----ke~gAvFi~gIG~~L~~G~~Hd~RapDYDD--W~t~~~~~~~gLNGDilv~n~~l~~a~El  246 (327)
                      .-|||+|.+ |+.||    ++||+-|+.---.+- +-+|. ---||.+|  ++         +.=|++|      +-+||
T Consensus        42 ~~~Dl~pe~-E~~i~~~a~~~~~sdfvfithyP~-~krpf-Ytmp~~~d~~~t---------~sfDlLF------rG~EI  103 (185)
T tr|A0A0G0WIP2|   42 EEPDLTPEE-EKLICKWAKKKYGSDFVFITHYPT-SKRPF-YTMPDPEDPEFT---------LSFDLLF------RGLEI  103 (185)
Confidence            348888865 55554    578887765322222 11111 11245555  77         6667766      57999


Q d12asa_         247 SSMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVR  320 (327)
Q Consensus       247 SSmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~  320 (327)
                      .+=|=|+. -+.|.+|.+..+.+...--.|-..-.-| +|.-=|=|+|.-||.|-||.-+.|.|+  |..|.++.
T Consensus       104 ttGsQRih~Y~~L~~~ik~kgl~~e~~k~Yl~~Fk~G-mPp~GG~~~GleRlt~~ll~~~nvrea--tlfpRD~~  175 (185)
T tr|A0A0G0WIP2|  104 TTGSQRIHNYDELVKNIKKKGLNPENFKFYLEAFKYG-MPPHGGLGFGLERLTMQLLGIDNVREA--TLFPRDMK  175 (185)
Confidence            99999985 5789999999998776666677777667 888889999999999999999999985  56665543


No 64
>tr|A0A0M1J520|A0A0M1J520_9GAMM Uncharacterized protein OS=Achromatium sp. WMS3 GN=TI05_18145 PE=4 SV=1
Probab=92.79  E-value=0.09  Score=41.10  Aligned_cols=35  Identities=43%  Similarity=0.613  Sum_probs=30.6  Template_Neff=1.000

Q d12asa_         277 QALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       277 ~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      .++..| ||.+-|-.||..||.|+|....||.||-+
T Consensus        23 ealeag-lpdcsgvaigldrllmwligarhidevla   57 (60)
T tr|A0A0M1J520|   23 EALEAG-LPDCSGVAIGLDRLLMWLIGARHIDEVLA   57 (60)
Confidence            344455 99999999999999999999999999965


No 65
>tr|A0A090PT73|A0A090PT73_9FLAO Lysyl-tRNA synthetase OS=Nonlabens ulvanivorans GN=JCM19298_1801 PE=3 SV=1
Probab=92.71  E-value=0.095  Score=45.81  Aligned_cols=38  Identities=45%  Similarity=0.638  Sum_probs=32.3  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .+|-+++.-| +|.|-|=|||-.||.|||-....|.||-
T Consensus         5 qdflraleyg-mpptsglgigmdrlimfltnnpsiqevl   42 (116)
T tr|A0A090PT73|    5 QDFLRALEYG-MPPTSGLGIGMDRLIMFLTNNPSIQEVL   42 (116)
Confidence            3455666656 8999999999999999999999999984


No 66
>tr|A0A0G7ZNA4|A0A0G7ZNA4_9MOLU Uncharacterized protein OS=Candidatus Hepatoplasma crinochetorum GN=HEPPS_03750 PE=4 SV=1
Probab=92.62  E-value=0.098  Score=52.99  Aligned_cols=290  Identities=14%  Similarity=0.225  Sum_probs=192.6  Template_Neff=1.320

Q d12asa_           3 IAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLS--GAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQ   80 (327)
Q Consensus         3 ~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLn--G~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~k   80 (327)
                      .|.|..|.+..-..-..+-.++.+.-|..|||....  +-|..+  ..-|.++|+..+. -.-.-+-.|..=|--.-+.+
T Consensus        15 leiq~linf~e~kii~l~y~~fd~hy~rkplfs~~k--~~d~ise~n~~rti~fdtind-yqiyn~yns~d~~f~kk~~~   91 (325)
T tr|A0A0G7ZNA4|   15 LEIQYLINFCEMKIIDLFYDRFDFHYIRKPLFSVDK--IIDQISEINYSRTISFDTIND-YQIYNIYNSFDFWFIKKISD   91 (325)
Confidence            355666777766667777889999999999997654  445554  3468999998764 34566788999999999999


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLP  159 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp  159 (327)
                      .++..+.|++|-.|.|-||- ++|.-.|+--.=+=+|..-+...--.+.....-++||+.+...=-++.---| ++-..|
T Consensus        92 lnlennhg~f~klnyivrd~-~id~nnsle~nilf~e~rsdn~~~~~~km~~lg~~iyqliydlffemhl~n~k~ki~ip  170 (325)
T tr|A0A0G7ZNA4|   92 LNLENNHGIFSKLNYIVRDS-KIDSNNSLERNILFFEMRSDNIENFFQKMQILGRQIYQLIYDLFFEMHLFNPKIKINIP  170 (325)
Confidence            99999999999999999998 8999999876666666665555555556666678899888876544432223 332333


Q d12asa_         160 --DQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWN  237 (327)
Q Consensus       160 --~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n  237 (327)
                        .+++|..--+.++.-.  +.+.-.+.+|-+||+|-+.   ..|.+....     |.  -+         .-..+-+|.
T Consensus       171 k~nk~~f~~~i~~~k~~q--~y~sfinqi~lk~gvvsvl---eel~~~e~e-----~f--ln---------~k~k~~fyy  229 (325)
T tr|A0A0G7ZNA4|  171 KINKVTFFNQIEFKKQNQ--NYQSFINQITLKYGVVSVL---EELNDNEDE-----NF--LN---------SKDKFSFYY  229 (325)
Confidence              4677776555555433  3455678899999987543   334332211     11  11         111222332


Q d12asa_         238 --PVLE---DAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       238 --~~l~---~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                        .-+.   +-+|+-.   |-..+.+....+ ..-..+++++|-..+.++---.||.--|...-|.+.-|.-+.+-|+|+
T Consensus       230 fsky~~enlkf~ei~k---rk~~~eil~~~~-k~i~nk~el~~ls~~i~dfky~tinikinldylilm~l~iasvfeiq~  305 (325)
T tr|A0A0G7ZNA4|  230 FSKYIRENLKFFEICK---RKEQSEILQEMQ-KNISNKIELDYLSEQIEDFKYQTINIKINLDYLILMVLNIASVFEIQS  305 (325)
Confidence              2222   2333322   222222222111 122346788888888888888899999999999999999999999998


Q d12asa_         313 GVWPAAVRE  321 (327)
Q Consensus       313 svW~~~~~~  321 (327)
                      ..=..++.+
T Consensus       306 n~nn~~v~e  314 (325)
T tr|A0A0G7ZNA4|  306 NSNNRKVNE  314 (325)
Confidence            766555443


No 67
>tr|U1QP96|U1QP96_9EURY Aspartyl/asparaginyl-tRNA synthetase OS=halophilic archaeon J07HX64 GN=J07HX64_00670 PE=3 SV=1
Probab=92.54  E-value=0.081  Score=47.51  Aligned_cols=68  Identities=31%  Similarity=0.350  Sum_probs=54.3  Template_Neff=4.000

Q d12asa_         243 AFELSSMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       243 a~ElSSmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      -.||.+=|-|++ .+.|.++++..+.+. ....|+-.+..=-.|.--|.|||--|+.|+|+.-.||.||-
T Consensus        32 ~~Ei~gGsqR~~~~~~L~~~~~~~gl~~-~~~~~Yld~~ryG~~PH~GfGlG~ER~~~~l~g~~nIRe~~  100 (111)
T tr|U1QP96|U1QP   32 GGEIIGGSQREHDYDILEERIKEKGLNP-EDYEWYLDLRRYGSPPHGGFGLGLERLLMWLLGLDNIREVI  100 (111)
Confidence            689999999986 467777888777644 45555555555557778899999999999999999999974


No 68
>tr|A0A0K0F1X5|A0A0K0F1X5_9BILA Uncharacterized protein OS=Strongyloides venezuelensis PE=3 SV=1
Probab=92.47  E-value=0.1  Score=46.71  Aligned_cols=41  Identities=39%  Similarity=0.531  Sum_probs=35.3  Template_Neff=2.067

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ..-.|.-++.-| ||.|-|=|.|..||.|||-.-+||.||--
T Consensus        50 ~DE~Fc~aLEyG-lpPTgGwGlGIDRl~mllTds~nikevll   90 (115)
T tr|A0A0K0F1X5|   50 LDETFCTALEYG-LPPTGGWGLGIDRLAMLLTDSQNIKEVLL   90 (115)
Confidence            445667777666 99999999999999999999999999953


No 69
>tr|W9LCI1|W9LCI1_FUSOX Aspartyl-tRNA synthetase OS=Fusarium oxysporum f. sp. lycopersici MN25 GN=FOWG_16590 PE=4 SV=1
Probab=92.45  E-value=0.086  Score=51.30  Aligned_cols=67  Identities=28%  Similarity=0.328  Sum_probs=53.3  Template_Neff=3.984

Q d12asa_         244 FELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       244 ~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .||.|=+-|| |.+.|.++.+..+-+...-..|-..-.-| -|..-|+|||-.|+.|++|.-.||.++-
T Consensus       100 ~EI~sGaqRih~~~~L~~~~e~~gi~~~~~~~Yid~f~yG-~pPHgG~GlGlER~~~~~~gl~nIR~~~  167 (185)
T tr|W9LCI1|W9LC  100 EEIVSGAQRIHDYEELRERAEEHGIDPESIYWYIDAFRYG-APPHGGFGLGLERLLMWLLGLGNIRDAS  167 (185)
Confidence            7899999998 56778888888777654444444444445 8999999999999999999999999863


No 70
>tr|W4YTD6|W4YTD6_STRPU Uncharacterized protein OS=Strongylocentrotus purpuratus PE=3 SV=1
Probab=92.39  E-value=0.12  Score=42.62  Aligned_cols=38  Identities=37%  Similarity=0.499  Sum_probs=32.3  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..|..++.-| ||.|=|-|+|..||.|||-....|.||-
T Consensus        14 enfctaleyg-lpptggfglgidrlamfltdsnnikevl   51 (80)
T tr|W4YTD6|W4YT   14 ENFCTALEYG-LPPTGGFGLGIDRLAMFLTDSNNIKEVL   51 (80)
Confidence            4455555555 9999999999999999999999999985


No 71
>tr|X0UCQ7|X0UCQ7_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S01H1_19828 PE=4 SV=1
Probab=92.25  E-value=0.1  Score=42.13  Aligned_cols=39  Identities=36%  Similarity=0.442  Sum_probs=30.6  Template_Neff=3.376

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..+|-.+-.--.|.--|.|||--|+.|++|...||.|+-
T Consensus        12 ~~~yl~~f~yG~pPHgG~GiGleR~~~~~lGl~nIRe~~   50 (61)
T tr|X0UCQ7|X0UC   12 FESYLDSFRYGMPPHGGFGIGLERIVMMLLGLKNIREAV   50 (61)
Confidence            334444444447778899999999999999999999974


No 72
>sp|O67258|SYK_AQUAE Lysine--tRNA ligase OS=Aquifex aeolicus (strain VF5) GN=lysS PE=3 SV=1
Probab=92.22  E-value=0.12  Score=57.32  Aligned_cols=263  Identities=25%  Similarity=0.373  Sum_probs=133.7  Template_Neff=2.276

Q d12asa_           7 RQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL---PDAQFEVVHSLAKWKRQTLGQHDF   83 (327)
Q Consensus         7 ~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~---~~~~~EIVhSLAKWKR~aL~ky~~   83 (327)
                      ++|+.++++|+.+     +.++|..|+.-.-.+|-|      -+|  |-.-..   .+...-|.--|. -||....  ||
T Consensus       275 k~I~~lr~fle~~-----gfiEVeTPilQpiasGA~------AkP--FiTyhN~Le~~lyLRIapELY-LKrLIVG--Gf  338 (597)
T sp|O67258|SYK_  275 KLIKELREFLESR-----GFIEVETPILQPIASGAN------AKP--FITYHNYLEQNLYLRIAPELY-LKRLIVG--GF  338 (597)
Confidence            5677777777654     789999999888887765      345  322111   022222222221 1222221  22


Q d12asa_          84 SAGEGLYTHMKALRPDEDRLSPLHS-------VYVDQWDWERVMGDGERQFSTLK-STVEAIWAGIKATEAAVSEEF---  152 (327)
Q Consensus        84 ~~geGiyTdMnAIRrDE~~ld~~HS-------iyVDQWDWEkvI~~~dRnl~~Lk-~tV~kIy~al~~te~~v~~~y---  152 (327)
                      +   -.|-=-...| .| -+|.+|-       .|.-=||+..++.=.++-..+|- ++|.+.-=.-..-|-.+...|   
T Consensus       339 ~---RVyElgknFR-NE-gvDtTHNPEFtMvEFY~AY~DY~dlm~~TE~lf~~~l~~~~g~lkit~~~~~ldf~~Pfk~~  413 (597)
T sp|O67258|SYK_  339 N---RVYELGKNFR-NE-GVDTTHNPEFTMVEFYAAYWDYKDLMKFTEELFSYLLKDTVGTLKITYQGQELDFKPPFRVY  413 (597)
Confidence            1   1111112223 34 4776664       68888999999987777665544 555443222222222222222   


Q d12asa_         153 ---G-LAPFLPDQIHFVH--SQELLSRYPDL--DAKGRER------A----IAK--DLGAVFLVGIGGKLSD-GHRHDVR  211 (327)
Q Consensus       153 ---~-l~~~Lp~~I~FI~--sqeL~~~YP~L--tpkeRE~------~----i~k--e~gAvFi~gIG~~L~~-G~~Hd~R  211 (327)
                         . ++.+-.++-.|+-  .++|.+.--.+  +..++..      +    ++.  --+.+|+|..-+.|+- .+.|- .
T Consensus       414 ~~fd~lkektgkDkdffl~~~e~~rkla~e~gip~a~~LTh~klidk~Fe~~vEe~L~qP~FvIdfPk~lsPLaKthR-~  492 (597)
T sp|O67258|SYK_  414 RYFDLLKEKTGKDKDFFLRDEEELRKLAKELGIPKAETLTHAKLIDKLFEHFVEEELIQPTFVIDFPKILSPLAKTHR-E  492 (597)
Confidence               1 1233333333332  33333322111  1111111      0    011  1356777777666642 12221 1


Q d12asa_         212 APDY-DDWS-----TPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDED----RLELEWHQALLR  281 (327)
Q Consensus       212 apDY-DDW~-----t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~----r~~~~~h~~ll~  281 (327)
                      .||- .-+.     .+..+.|.-||           +.+        .-++.+.+|++....-+    -..-+|-+++.-
T Consensus       493 dpdLVERFEL~i~~~EiaNaYtELN-----------dP~--------dQrerfleQ~kEk~mGdEEAm~mDEdFi~Aley  553 (597)
T sp|O67258|SYK_  493 DPDLVERFELFINKKEIANAYTELN-----------DPI--------DQRERFLEQLKEKEMGDEEAMDMDEDFITALEY  553 (597)
Confidence            1221 0111     00111122222           111        23567888887654433    256678888888


Q d12asa_         282 GEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       282 ~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      | +|.|.|-|||..||+|+|-.-..|.||-
T Consensus       554 G-mPPTaGeGIGiDRLvM~l~~~~SIreVI  582 (597)
T sp|O67258|SYK_  554 G-MPPTAGEGIGIDRLVMLLADVDSIREVI  582 (597)
Confidence            8 9999999999999999999999999985


No 73
>tr|A0A0E0KFX6|A0A0E0KFX6_ORYPU Lysine--tRNA ligase OS=Oryza punctata PE=3 SV=1
Probab=92.19  E-value=0.1  Score=59.13  Aligned_cols=56  Identities=30%  Similarity=0.451  Sum_probs=45.9  Template_Neff=3.600

Q d12asa_         255 ADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+.+|.+....-+    -....|..+|.-| ||.|-|=|+|..||+|||-...+|.||.
T Consensus       556 R~rFe~Qa~~k~~GDdEA~~~De~F~~ALEYG-LPPTgGwG~GIDRLvMlLTd~~~IkeV~  615 (652)
T tr|A0A0E0KFX6|  556 RERFEEQAKDKAAGDDEAMPIDEDFCTALEYG-LPPTGGWGMGIDRLVMLLTDSNNIKEVE  615 (652)
Confidence            456666765543332    3677899999998 9999999999999999999999999993


No 74
>tr|B7BCA8|B7BCA8_9PORP tRNA ligases class II (D, K and N) (Fragment) OS=Parabacteroides johnsonii DSM 18315 GN=PRABACTJOHN_02676 PE=3 SV=1
Probab=92.08  E-value=0.14  Score=48.69  Aligned_cols=56  Identities=38%  Similarity=0.540  Sum_probs=44.7  Template_Neff=1.000

Q d12asa_         255 ADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       255 ~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ++.+++||..+...+    ....+|-+++.=| +|.|-|-|||..||.|++-....|.||-
T Consensus        66 eerfkeqlrlsekgddeamfidqdflralqfg-mpptsgigigidrlvmlmtgqttiqevl  125 (208)
T tr|B7BCA8|B7BC   66 EERFKEQLRLSEKGDDEAMFIDQDFLRALQFG-MPPTSGIGIGIDRLVMLMTGQTTIQEVL  125 (208)
Confidence            456777877665443    2556777787777 8999999999999999999999999984


No 75
>tr|B7XRA1|B7XRA1_ENTBH Aspartate-ammonia ligase OS=Enterocytozoon bieneusi (strain H348) GN=EBI_24170 PE=4 SV=1
Probab=91.95  E-value=0.15  Score=43.11  Aligned_cols=80  Identities=40%  Similarity=0.641  Sum_probs=61.2  Template_Neff=1.000

Q d12asa_          40 GTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDR----LSPLHSVYVDQWD  115 (327)
Q Consensus        40 GlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~----ld~~HSiyVDQWD  115 (327)
                      |..|.|.|.|.   |+......--..||||||||||-.|-.-.-..-.|++..|+.+|.-||.    -|..||.-|.|+|
T Consensus         3 gvnddlngeep---vrfqtsdgewcsvvhslakwkrwmlwklrddgvsgiwcdmrgirkcedadvictsrmhsyqveqfd   79 (93)
T tr|B7XRA1|B7XR    3 GVNDDLNGEEP---VRFQTSDGEWCSVVHSLAKWKRWMLWKLRDDGVSGIWCDMRGIRKCEDADVICTSPMHSYQVEQFD   79 (93)
Confidence            44567777653   3344444455679999999999888665555667899999999987774    5889999999999


Q d12asa_         116 WERVMGD  122 (327)
Q Consensus       116 WEkvI~~  122 (327)
                      ||.|...
T Consensus        80 wekvips   86 (93)
T tr|B7XRA1|B7XR   80 WEKVNPS   86 (93)
Confidence            9999754


No 76
>tr|A0A0G1IIJ4|A0A0G1IIJ4_9BACT Asparagine-tRNA ligase OS=Parcubacteria (Giovannonibacteria) bacterium GW2011_GWA1_44_25 GN=UW53_C0025G0001 PE=3 SV=1
Probab=91.90  E-value=0.12  Score=55.50  Aligned_cols=80  Identities=25%  Similarity=0.184  Sum_probs=56.3  Template_Neff=4.000

Q d12asa_         227 AGLNGDI-LVWNPVLEDAFELSSMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       227 ~gLNGDi-lv~n~~l~~a~ElSSmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k  304 (327)
                      .-+|+|+ ++.-      =|+--+|=|+. .+.+.+-|+.-+....--..|-+|=..-.+-.|-|-|+|-=|+.|++|+.
T Consensus       291 ~A~NADliL~G~------GEvVGsGeR~~~~~e~~~sl~~h~v~~~dY~wYl~mr~~~p~~~TSGFGmGvERfi~WlL~~  364 (382)
T tr|A0A0G1IIJ4|  291 KALNADLILFGI------GEVVGSGERHTTAEETLESLERHNVPPEDYEWYIEMREFKPLKRTSGFGLGVERFILWLLQH  364 (382)
Confidence            5599999 5443      35666788864 56666667666664444444444433332336999999999999999999


Q d12asa_         305 PHIGQVQA  312 (327)
Q Consensus       305 ~HIgEVq~  312 (327)
                      ..|.+||.
T Consensus       365 ~dIrd~ql  372 (382)
T tr|A0A0G1IIJ4|  365 DDIRDVQL  372 (382)
Confidence            99999984


No 77
>tr|A0A0G1D8P3|A0A0G1D8P3_9BACT Asparagine-tRNA ligase (Fragment) OS=Parcubacteria bacterium GW2011_GWE2_43_12 GN=UV69_C0002G0045 PE=3 SV=1
Probab=91.54  E-value=0.13  Score=52.39  Aligned_cols=203  Identities=18%  Similarity=0.190  Sum_probs=127.8  Template_Neff=4.800

Q d12asa_          86 GEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKA----TEAAVSEEFGLAPFL---  158 (327)
Q Consensus        86 geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~----te~~v~~~y~l~~~L---  158 (327)
                      =.+.||=+-.-|..+ ..+.-|=.-.-.+|-|.-+...+-.++-+.+.+..|.+.+.+    --+.+...++..+.+   
T Consensus        22 ~~~Vf~igp~FRAE~-S~T~RHL~Ef~~lE~E~af~~~~dvm~~~E~~~~~i~~~vl~~~~~el~~~~~~~~~~~~~~~p  100 (256)
T tr|A0A0G1D8P3|   22 FGRVFTIGPTFRAEN-SNTRRHLAEFWMLEAEMAFIDYHDVMDLAEELIKYIFKYVLENCKEELEFFNKDYPFLPKIKKP  100 (256)
Confidence            378888888888665 566677667777777776666556665555555555555533    223333333311111   


Q d12asa_         159 PDQIHFVHSQELLSR-------YPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRA------PDYDDWSTPSELG  225 (327)
Q Consensus       159 p~~I~FI~sqeL~~~-------YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~Ra------pDYDDW~t~~~~~  225 (327)
                      -..|+|-.+-++++.       -.||+ .+-|+.+++.+|..+++       ..-|++.|+      |+...|+      
T Consensus       101 f~ritY~EAi~~L~e~~~~~~~g~Dl~-te~ER~L~e~~~~pv~v-------~~yP~~ikpFYm~~~p~d~k~~------  166 (256)
T tr|A0A0G1D8P3|  101 FPRITYTEAIEILKEAGDPVEWGDDLS-TEHERYLGELFGTPVFV-------TNYPKEIKPFYMRPNPEDPKTV------  166 (256)
Confidence            134555555555533       23443 24566777776644222       344555554      3334455      


Q d12asa_         226 HAGLNGDILVWNPVLEDAFELSSMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       226 ~~gLNGDilv~n~~l~~a~ElSSmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k  304 (327)
                         ..-|++-    ....-||-+-|=|++ .+.|.++++..+-+...-..|-.+-.=| -|..-|+|+|-=|++|+++.-
T Consensus       167 ---~s~DlL~----~~GvgEiigGsqRe~~~~~L~~r~~~~gl~~e~y~wYlDlrkyG-~~pH~GfGlG~ERlv~~~~gl  238 (256)
T tr|A0A0G1D8P3|  167 ---ASMDLLA----VPGIGEIIGGSQREHDYDELEERMKELGLDPEDYEWYLDLRKYG-SPPHGGFGLGFERLVMWLLGL  238 (256)
Confidence               6666661    133578999999998 7888999999888654444444444445 666789999999999999999


Q d12asa_         305 PHIGQVQ  311 (327)
Q Consensus       305 ~HIgEVq  311 (327)
                      .||.+|+
T Consensus       239 ~nIRd~~  245 (256)
T tr|A0A0G1D8P3|  239 DNIRDVI  245 (256)
Confidence            9999986


No 78
>tr|M8AA95|M8AA95_TRIUA Uncharacterized protein OS=Triticum urartu PE=3 SV=1
Probab=91.33  E-value=0.18  Score=53.30  Aligned_cols=41  Identities=44%  Similarity=0.604  Sum_probs=35.2  Template_Neff=2.477

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ...+|..++.-| +|.|-|=|+|..||.|||-.-.||.||-+
T Consensus       320 ~Dedf~~aLEyg-mpPtaGwGlGIDRlvmllTn~~~IreVl~  360 (368)
T tr|M8AA95|M8AA  320 LDEDFCTALEYG-MPPTSGWGLGIDRLVMLLTNSASIREVLA  360 (368)
Confidence            344566777666 99999999999999999999999999964


No 79
>tr|A0A078CAH0|A0A078CAH0_BRANA Lysine--tRNA ligase OS=Brassica napus GN=BnaA03g31900D PE=3 SV=1
Probab=91.25  E-value=0.18  Score=57.79  Aligned_cols=58  Identities=29%  Similarity=0.417  Sum_probs=46.9  Template_Neff=3.200

Q d12asa_         254 DADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       254 d~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      -++.+.+|.+....-+    -..-.|..+|.-| ||.|-|=|+|..||+|||-....|.||..
T Consensus       536 QR~RFeeQak~k~aGDDEA~~iDE~Fc~ALEYG-LPPTgGWG~GIDRL~MfLTds~nIKEVLl  597 (720)
T tr|A0A078CAH0|  536 QRERFEEQAKQKAAGDDEAQMIDENFCTALEYG-LPPTGGWGMGIDRLTMFLTDSNNIKEVLL  597 (720)
Confidence            3566777776544433    2566788888888 99999999999999999999999999974


No 80
>tr|A0A0G1M6Q9|A0A0G1M6Q9_9BACT Aspartyl-tRNA synthetase (Fragment) OS=Microgenomates (Woesebacteria) bacterium GW2011_GWE1_45_18 GN=UX03_C0006G0024 PE=3 SV=1
Probab=91.18  E-value=0.2  Score=51.70  Aligned_cols=228  Identities=21%  Similarity=0.240  Sum_probs=128.0  Template_Neff=2.381

Q d12asa_          61 DAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSV-YVDQWDWERVMGDGERQFSTLKSTVEAIWA  139 (327)
Q Consensus        61 ~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld~~HSi-yVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~  139 (327)
                      +...-.-||-..+|. .|.-.||..+--|   -+..| || +.-.-.+. -.-|+|.|.-.-.++--++........+.+
T Consensus        16 gkfyalpqspqqykq-llmvaGfErYfQI---arCfR-DE-D~R~DRqp~EftQlD~EMSFv~~~dil~l~E~m~~~lv~   89 (299)
T tr|A0A0G1M6Q9|   16 GKFYALPQSPQQYKQ-LLMVAGFERYFQI---ARCFR-DE-DARADRQPGEFTQLDMEMSFVTQEDVLDLTEGMFIELVE   89 (299)
Confidence            556666777777764 4555666544333   34444 77 56666676 788999998777666554443333322222


Q d12asa_         140 GIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGA----------------VFLVGIG-GKL  202 (327)
Q Consensus       140 al~~te~~v~~~y~l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gA----------------vFi~gIG-~~L  202 (327)
                      .+- .|+.+.     ...+|.    ++.+|-.+.|-.=-|.-|+++=.+.-=|                -|.+|-+ +++
T Consensus        90 ~~f-p~k~i~-----~~PfPr----lt~~e~mekyg~DKpDlR~~k~d~~elaFaWv~dfPlFe~~~~~~~~~~~~~~~~  159 (299)
T tr|A0A0G1M6Q9|   90 KLF-PEKKIT-----QSPFPR----LTYEESMEKYGTDKPDLRKNKNDPNELAFAWVVDFPLFEEQEAEDFFHGSGEKKW  159 (299)
Confidence            110 111111     112221    3444555555544444455443221111                1222223 222


Q d12asa_         203 SDGHRHDVRAPDYDD--WSTPSELGHAGLNGDIL--VWNPVLEDAFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQ  277 (327)
Q Consensus       203 ~~G~~Hd~RapDYDD--W~t~~~~~~~gLNGDil--v~n~~l~~a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~  277 (327)
                      . -.-|---||-.+|  |-...       -|+|.  -|+-+|+ -.||++=+||. +++..++-.++.+-+...+..|-+
T Consensus       160 ~-~~HhpFtapk~ed~~~l~~~-------~~~vr~~qyDlVlN-GyEiggGsiRih~pe~q~kif~~~g~s~~~~~~fgh  230 (299)
T tr|A0A0G1M6Q9|  160 A-PSHHPFTAPKEEDLELLMNE-------PGKVRAYQYDLVLN-GYEIGGGSIRIHDPEIQEKIFDLIGFSQEQKKRFGH  230 (299)
Confidence            1 1122233333333  33111       11221  3444554 57999999998 566667778888888888899999


Q d12asa_         278 ALLRGEM--PQTIGGGIGQSRLTMLLLQLPHIGQVQAG  313 (327)
Q Consensus       278 ~ll~~~l--P~TIgGGIgqSRl~M~lL~k~HIgEVq~s  313 (327)
                      +|.+=+.  |.-=|=+.|-.||.|.|+....|.||+|=
T Consensus       231 mleAf~yG~PPHgGia~G~DR~vm~l~gE~~irEviaF  268 (299)
T tr|A0A0G1M6Q9|  231 MLEAFEYGAPPHGGIAPGLDRLVMILCGEPNIREVIAF  268 (299)
Confidence            9987655  44444455899999999999999999973


No 81
>tr|Q9FEL1|Q9FEL1_TOBAC Lysyl-tRNA synthetase (Fragment) OS=Nicotiana tabacum GN=lysRS PE=2 SV=1
Probab=90.95  E-value=0.24  Score=46.63  Aligned_cols=58  Identities=34%  Similarity=0.449  Sum_probs=45.3  Template_Neff=1.461

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+|++.+..-+.    +.-.|.-++.-| ||.|=|=|.|..||.|||-.-+.|.||-
T Consensus        55 ~qrqrf~~q~k~rq~gddea~~lde~fc~aleyg-l~ptggwglgidrl~m~lt~s~nikevl  116 (169)
T tr|Q9FEL1|Q9FE   55 VQRQRFAEQLKDRQSGDDEAMALDETFCTALEYG-LPPTGGWGLGIDRLTMFLTDSQNIKEVL  116 (169)
Confidence            456677788876654432    445566666666 9999999999999999999999999984


No 82
>tr|G7DZZ8|G7DZZ8_MIXOS Uncharacterized protein OS=Mixia osmundae (strain CBS 9802 / IAM 14324 / JCM 22182 / KY 12970) GN=Mo02819 PE=3 SV=1
Probab=90.94  E-value=0.25  Score=53.60  Aligned_cols=112  Identities=34%  Similarity=0.493  Sum_probs=70.9  Template_Neff=1.000

Q d12asa_         191 GAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG-----HAGLNGDILVWN-----------------PVLED------  242 (327)
Q Consensus       191 gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~-----~~gLNGDilv~n-----------------~~l~~------  242 (327)
                      |.-++-.=|++- ||+    .--||||++||++..     ....|.|-++.+                 |.+..      
T Consensus       515 gvrllredgwkd-ngq----eigdyddfstptekrlgeivkakfntdyyildkfplavrpfytmpdgqdpklsnsydffm  589 (691)
T tr|G7DZZ8|G7DZ  515 GVRLLREDGWKD-NGQ----EIGDYDDFSTPTEKRLGEIVKAKFNTDYYILDKFPLAVRPFYTMPDGQDPKLSNSYDFFM  589 (691)
Confidence            444555445543 333    346999999999873     334566655332                 22222      


Q d12asa_         243 -AFELSSMGIRVDADTLKH-QLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIG  308 (327)
Q Consensus       243 -a~ElSSmGirVd~~~L~~-Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIg  308 (327)
                       .-||-|-.-||....|.+ ....-+-..+.-.+|-.-..=| -|.-.|||||.-|+.|++|.--.|.
T Consensus       590 rgeeilsgaqrvhtadllesrmrevgidpkemeeyvdgfrwg-tpphagggiglervlmlflklgnir  656 (691)
T tr|G7DZZ8|G7DZ  590 RGEEILSGAQRVHTADLLESRMREVGIDPKEMEEYVDGFRWG-TPPHAGGGIGLERVLMLFLKLGNIR  656 (691)
Confidence             235666667777665543 4444455445555666666666 5778899999999999999886664


No 83
>tr|G5DY35|G5DY35_9PIPI Putative lysyl-trna synthetase (Fragment) OS=Hymenochirus curtipes PE=2 SV=1
Probab=90.86  E-value=0.21  Score=52.14  Aligned_cols=206  Identities=22%  Similarity=0.311  Sum_probs=118.2  Template_Neff=3.426

Q d12asa_          90 YTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEF-------G---L--APF  157 (327)
Q Consensus        90 yTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y-------~---l--~~~  157 (327)
                      +.+-+..|..+  +|..|..-+--.+|=.+-....-.++..+..+..+...+  ....|.  |       +   +  .+.
T Consensus        42 ~~~gr~frnEg--~d~~hnpeF~~~e~y~~YaDy~dlm~~~e~l~~~~v~~i--~~~~i~--~~~~~~~~~~~~i~f~~p  115 (302)
T tr|G5DY35|G5DY   42 MANGRQFRNEG--IDLTHNPEFTTCEFYQAYADYNDLMDLTEKLLSSMVKHI--GDKKIT--YHPEGPENDKHEIDFTPP  115 (302)
Confidence            44556677554  888898888888888887777777777777777777666  223332  3       1   1  111


Q d12asa_         158 LPDQIHFVHSQE---LLSRYPDLDA-----KGRERAIAKDLGAV---------FLVGI-GGKLSDGHRHDVRAPDYDDWS  219 (327)
Q Consensus       158 Lp~~I~FI~sqe---L~~~YP~Ltp-----keRE~~i~ke~gAv---------Fi~gI-G~~L~~G~~Hd~RapDYDDW~  219 (327)
                      + ..|+|++.=+   +.+..+..++     ..-...+|+++|.-         ++..| |+.+.-...+-..-=||.-=.
T Consensus       116 f-~ri~~~~~l~~~~~~~~~~~~~~~~~~~~~~L~~~~~~~gi~~~~~~s~~~l~d~i~~~~~ep~~i~PTFiid~P~~~  194 (302)
T tr|G5DY35|G5DY  116 F-RRITMIEALEKYAGIKLPAKCSIDTDESRAELEDICKRFGIECPPPESTARLLDKIFKEFLEPKCINPTFIIDYPQIM  194 (302)
Confidence            1 2344444332   2234444442     33335556665532         11111 111112222222223333222


Q d12asa_         220 TPSELGHAGLNGDILVWNPVLEDAFELSSMGIRV------------DADTLKHQLALTGDEDR----LELEWHQALLRGE  283 (327)
Q Consensus       220 t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirV------------d~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~  283 (327)
                      +|...-..        =||.+-.++||--+|+-+            -++.+.+|++..++.+.    ..-.|--++.=| 
T Consensus       195 sPlAK~hr--------~~~~l~erf~l~~~g~EI~na~sELNDPl~Qr~rf~~Q~~~~~~gd~EA~~~De~F~~ALEYG-  265 (302)
T tr|G5DY35|G5DY  195 SPLAKWHR--------SNPGLVERFELFVFGKEIANAYSELNDPLEQRERFEEQAKQKKAGDDEAMPIDEDFCEALEYG-  265 (302)
Confidence            22211111        134455566666666654            24677888877665543    566788888777 


Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +|.|-|=|||.-||+|+|..-..|.||-
T Consensus       266 mPP~ag~g~gIdrL~m~Ltd~~~i~evl  293 (302)
T tr|G5DY35|G5DY  266 MPPTAGWGIGIDRLTMLLTDSPNIKEVI  293 (302)
Confidence            9999999999999999999999999884


No 84
>tr|A0A079FAQ6|A0A079FAQ6_ECOLX tRNA synthetases class II family protein OS=Escherichia coli 2-460-02_S3_C3 GN=AC49_5683 PE=3 SV=1
Probab=90.46  E-value=0.31  Score=39.85  Aligned_cols=39  Identities=36%  Similarity=0.526  Sum_probs=31.6  Template_Neff=1.000

Q d12asa_         272 ELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       272 ~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..+|-.++.-| ||.|.|=|||..|+.|++-....|..|-
T Consensus        27 dedyvtaleyg-lpptaglgigidrmimlftnshtirdvi   65 (74)
T tr|A0A079FAQ6|   27 DEDYVTALEYG-LPPTAGLGIGIDRMIMLFTNSHTIRDVI   65 (74)
Confidence            34555566555 9999999999999999999888888764


No 85
>tr|R5CXA5|R5CXA5_9FIRM Aspartate--ammonia ligase OS=Firmicutes bacterium CAG:555 GN=BN705_00405 PE=4 SV=1
Probab=90.45  E-value=0.31  Score=38.72  Aligned_cols=28  Identities=29%  Similarity=0.388  Sum_probs=25.7  Template_Neff=1.000

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIE   28 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~r   28 (327)
                      |+.+||+||..+|..|+.+|..-|||..
T Consensus        14 svydtqkaisllkrlfedrlgallnl~~   41 (63)
T tr|R5CXA5|R5CX   14 SVYDTQKAISLLKRLFEDRLGALLNLXX   41 (63)
Confidence            5779999999999999999999999864


No 86
>tr|A0A0G0IMH4|A0A0G0IMH4_9BACT Elongation factor P-(R)-beta-lysine ligase OS=Microgenomates (Roizmanbacteria) bacterium GW2011_GWC2_37_13 GN=US40_C0008G0013 P
Probab=90.37  E-value=0.26  Score=53.23  Aligned_cols=276  Identities=20%  Similarity=0.284  Sum_probs=147.3  Template_Neff=3.200

Q d12asa_           6 QRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQ-------TL   78 (327)
Q Consensus         6 q~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~-------aL   78 (327)
                      .+.|+.|..+|.+     -|..+|..|+++.... --=.|+    +....+..   ....-.++=.|+|.+       +|
T Consensus        88 ~kVi~~IR~FFk~-----~gf~EVetP~L~p~p~-~EpyLe----vFeT~l~~---~~~~~~~~~~~~k~yL~tSPEy~L  154 (404)
T tr|A0A0G0IMH4|   88 EKVIDAIREFFKK-----QGFHEVETPILVPAPS-PEPYLE----VFETELRD---RKEKTALPGKKEKAYLITSPEYAL  154 (404)
Confidence            3456666666654     5889999999887541 111111    11111111   111112223344443       44


Q d12asa_          79 GQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPF  157 (327)
Q Consensus        79 ~ky~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~  157 (327)
                      +|.=-..-.-+|.=-.+-|..| ..+.+|+--+-=+.|=+|=...-.-++-..+....|...++.-+..+  .|. -.-.
T Consensus       155 KKLLaaG~Gncf~i~K~FRN~E-~~s~~HNPEFTMLEWYRv~adY~~im~D~e~Ll~~I~~~l~~~~~~~--~Yqgk~~~  231 (404)
T tr|A0A0G0IMH4|  155 KKLLAAGFGNCFEITKSFRNEE-SFSSLHNPEFTMLEWYRVGADYFDIMDDVEKLFLFINKKLKPKADKL--KYQGKKID  231 (404)
Confidence            4432222234677778889888 89999999999999998877776666666666666666665443222  132 1112


Q d12asa_         158 LPDQIHFVHSQELLSRYPDLDAKGRE-----RAIAKDLG-------------AVFLVGIGGKLS-DGHRHDVRAPDYDDW  218 (327)
Q Consensus       158 Lp~~I~FI~sqeL~~~YP~LtpkeRE-----~~i~ke~g-------------AvFi~gIG~~L~-~G~~Hd~RapDYDDW  218 (327)
                      |. .-.=||-+|.-++|-++++.+-.     .+.|+++|             -+|+-.|-..|. +|+|-=.     -|+
T Consensus       232 l~-pweriSv~eAF~Kyagi~~d~~~~~~~l~~~a~~kGY~v~~~t~edlf~qIflneIEP~Lg~~~~Ptii-----YDY  305 (404)
T tr|A0A0G0IMH4|  232 LS-PWERISVAEAFQKYAGIDLDELFDEEKLKKKAKKKGYKVANTTWEDLFYQIFLNEIEPHLGKQGKPTII-----YDY  305 (404)
Confidence            22 33346677777777777755542     33444432             123333333332 1221100     000


Q d12asa_         219 STPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDA-----DTLKHQLALTGD---EDRLELEWHQALLRGEMPQTIGG  290 (327)
Q Consensus       219 ~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~-----~~L~~Ql~~~~~---~~r~~~~~h~~ll~~~lP~TIgG  290 (327)
                      ..+-. ....++     =|+..-..+|+=-.|+-+--     ..-++|.+....   ......+|-.+|..| ||.|-|-
T Consensus       306 P~~~A-ALAK~~-----~D~r~AERFE~YiaGlELgN~fsELtD~kEQ~kRf~~eKt~y~ID~dFIeALk~G-lP~~sGI  378 (404)
T tr|A0A0G0IMH4|  306 PKPQA-ALAKKK-----SDPRFAERFEVYIAGLELGNCFSELTDAKEQEKRFKAEKTVYPIDEDFIEALKSG-LPECSGI  378 (404)
Confidence            00000 000011     12333344444444544311     123444444443   334556667777777 9999999


Q d12asa_         291 GIGQSRLTMLLLQLPHIGQV  310 (327)
Q Consensus       291 GIgqSRl~M~lL~k~HIgEV  310 (327)
                      +||--||.|++..-..|.++
T Consensus       379 AlGvDRLvMLf~d~~sI~d~  398 (404)
T tr|A0A0G0IMH4|  379 AVGVDRLVMLFADVKNIDDT  398 (404)
Confidence            99999999999999988876


No 87
>tr|A0A095BRF7|A0A095BRF7_SCHHA Lysine--tRNA ligase OS=Schistosoma haematobium GN=MS3_11247 PE=3 SV=1
Probab=90.22  E-value=0.34  Score=46.66  Aligned_cols=58  Identities=31%  Similarity=0.414  Sum_probs=41.3  Template_Neff=1.198

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+|.+...+.+.    ..-.|-.++--| ||.|.|=|+|..||.|||-....|.||-
T Consensus        39 vqrerf~eqakdkaagd~eamvtdeqfltalgyg-lpptagwglgidrltmfltntnnikevi  100 (206)
T tr|A0A095BRF7|   39 VQRERFMEQAKDKAAGDAEAMVTDEQFLTALGYG-LPPTAGWGLGIDRLTMFLTNTNNIKEVI  100 (206)
Confidence            456667777665433321    222344444444 9999999999999999999999999985


No 88
>tr|X1CXE8|X1CXE8_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S01H4_56259 PE=4 SV=1
Probab=90.14  E-value=0.35  Score=35.06  Aligned_cols=31  Identities=32%  Similarity=0.625  Sum_probs=28.2  Template_Neff=1.000

Q d12asa_         282 GEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       282 ~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +.||.+.|...|..||.|++-....|.||-.
T Consensus         3 ndlpeaagnalgidrltmlftdsrqidevvt   33 (39)
T tr|X1CXE8|X1CX    3 NDLPEAAGNALGIDRLTMLFTDSRQIDEVVT   33 (39)
Confidence            4599999999999999999999999999853


No 89
>tr|A0A0A0L4F7|A0A0A0L4F7_CUCSA Uncharacterized protein OS=Cucumis sativus GN=Csa_4G619050 PE=4 SV=1
Probab=89.55  E-value=0.44  Score=39.93  Aligned_cols=40  Identities=43%  Similarity=0.566  Sum_probs=32.3  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +..+|-.++.-| +|.+-|=|||..||.|+|-.-+.|..|-
T Consensus        43 ldddfltaleyg-mppasglgigidrlvmlltnsasirdvi   82 (84)
T tr|A0A0A0L4F7|   43 LDDDFLTALEYG-MPPASGLGIGIDRLVMLLTNSASIRDVI   82 (84)
Confidence            444555566555 8999999999999999999999888774


No 90
>tr|A0A059WSR4|A0A059WSR4_9BACT tRNA synthetases class II (D, K and N) (Fragment) OS=uncultured bacterium PE=3 SV=1
Probab=88.88  E-value=0.47  Score=45.98  Aligned_cols=46  Identities=39%  Similarity=0.568  Sum_probs=39.9  Template_Neff=2.564

Q d12asa_         275 WHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRE  321 (327)
Q Consensus       275 ~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~~~~  321 (327)
                      |-.++..| ||.+.|=-.|-.||.|+-+.-.||.+|-+-.||.++.+
T Consensus       131 fl~al~~g-LPdC~GVAlG~DRL~Mla~g~~~i~~v~af~~~~~~~~  176 (176)
T tr|A0A059WSR4|  131 FLAALEAG-LPDCTGVALGFDRLVMLALGANTIADVTAFPWPHATME  176 (176)
Confidence            34455555 99999999999999999999999999999999987653


No 91
>tr|A0A0J9VCE3|A0A0J9VCE3_PLAVI Lysine-tRNA ligase OS=Plasmodium vivax Brazil I GN=PVBG_00933 PE=3 SV=1
Probab=88.85  E-value=0.49  Score=53.53  Aligned_cols=45  Identities=36%  Similarity=0.584  Sum_probs=38.2  Template_Neff=2.239

Q d12asa_         266 GDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       266 ~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +.+.....+|--+|.-| ||.|=|=|||..||||+|-.-.-|..|-
T Consensus       625 n~n~eiDyDyvTaLahG-lPPtgGlGiGIDRLcMLlTnsssIknvl  669 (677)
T tr|A0A0J9VCE3|  625 SDNHEIDYDYVTALAHG-LPPTGGLGIGIDRLCMLLTNSSSIKNIL  669 (677)
Confidence            34445788899999998 9999999999999999999888887664


No 92
>tr|A0A0G2UQE3|A0A0G2UQE3_THAPS Lysine--tRNA ligase OS=Thalassiosira pseudonana PE=2 SV=1
Probab=88.84  E-value=0.44  Score=54.68  Aligned_cols=59  Identities=36%  Similarity=0.464  Sum_probs=48.5  Template_Neff=3.470

Q d12asa_         253 VDADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      |-++.+.+|++.....+    ...-.|-++|..| +|.|=|=|||..||.|||=..+.|.||-+
T Consensus       599 ~QRe~Fe~Q~~~ka~GD~Ea~~iDE~Fl~ALE~g-mPPtGG~GiGIDRLVMlLT~s~sIkeVi~  661 (669)
T tr|A0A0G2UQE3|  599 VQRERFEKQARDKAAGDEEAMDIDEDFLQALEYG-MPPTGGWGIGIDRLVMLLTNSANIREVIL  661 (669)
Confidence            45677888887655333    3567888898888 99999999999999999999999999864


No 93
>tr|G1XRE9|G1XRE9_ARTOA Lysine--tRNA ligase OS=Arthrobotrys oligospora (strain ATCC 24927 / CBS 115.81 / DSM 1491) GN=AOL_s00193g148 PE=3 SV=1
Probab=88.65  E-value=0.44  Score=54.46  Aligned_cols=42  Identities=36%  Similarity=0.520  Sum_probs=38.2  Template_Neff=4.000

Q d12asa_         270 RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       270 r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ...-.|-+++.=| ||.|-|=|+|-.||||+|=.+..|++|..
T Consensus       546 ~~De~Yl~ALEwG-LPPTGGWGcGIDRLvMLftG~~RI~dVL~  587 (626)
T tr|G1XRE9|G1XR  546 EVDESYLEALEWG-LPPTGGWGCGIDRLVMLFTGAKRIGDVLP  587 (626)
Confidence            4667888999888 99999999999999999999999999964


No 94
>tr|A0A0G4LDB2|A0A0G4LDB2_9PEZI Uncharacterized protein (Fragment) OS=Verticillium longisporum GN=BN1708_003268 PE=3 SV=1
Probab=88.64  E-value=0.56  Score=41.64  Aligned_cols=40  Identities=33%  Similarity=0.427  Sum_probs=32.8  Template_Neff=1.463

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAG  313 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~s  313 (327)
                      -.|.-.+.-| ||.|.|=|.|-.|+.|||-....|.||-+-
T Consensus        39 en~c~~leyg-lpptagwg~g~~r~~mfltd~~~i~evl~f   78 (110)
T tr|A0A0G4LDB2|   39 ENFCMSLEYG-LPPTAGWGMGIDRMVMFLTDHYTIREVLAF   78 (110)
Confidence            3444444444 999999999999999999999999999753


No 95
>tr|Q05FQ1|Q05FQ1_CARRP Putative phenylalanyl-tRNA synthetase alpha subunit OS=Carsonella ruddii (strain PV) GN=CRP_089 PE=4 SV=1
Probab=88.21  E-value=0.62  Score=40.04  Aligned_cols=25  Identities=40%  Similarity=0.632  Sum_probs=22.9  Template_Neff=1.724

Q d12asa_         286 QTIGGGIGQSRLTMLLLQLPHIGQV  310 (327)
Q Consensus       286 ~TIgGGIgqSRl~M~lL~k~HIgEV  310 (327)
                      ..+.||||-.||+|..-.+.||.+|
T Consensus        61 ~g~aggigidrl~mi~k~~~~ik~i   85 (87)
T tr|Q05FQ1|Q05F   61 KGIAGGIGIDRLIMIKKKKKHIKYI   85 (87)
Confidence            4688999999999999999999886


No 96
>sp|P43151|GPA_LEIDO Putative guanine nucleotide-binding protein subunit alpha OS=Leishmania donovani PE=2 SV=1
Probab=87.87  E-value=0.73  Score=48.61  Aligned_cols=58  Identities=29%  Similarity=0.387  Sum_probs=44.9  Template_Neff=1.000

Q d12asa_         253 VDADTLKHQLALTGDEDRL----ELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r~----~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+||..+...+..    ...|-.+ ++..||.|=|=|+|..||.|||-....|.||-
T Consensus       258 vqreefvkqlrnrdkgddesmeidegfvaa-lehalpptggwglgidrlvmfltsqsnikevl  319 (464)
T sp|P43151|GPA_  258 VQREEFVKQLRNRDKGDDESMEIDEGFVAA-LEHALPPTGGWGLGIDRLVMFLTSQSNIKEVL  319 (464)
Confidence            5678888999877654432    2334433 45669999999999999999999999999985


No 97
>tr|B9Y6V6|B9Y6V6_9FIRM Aspartate--tRNA ligase OS=Holdemania filiformis DSM 12042 GN=aspS PE=3 SV=1
Probab=87.86  E-value=0.52  Score=53.61  Aligned_cols=76  Identities=26%  Similarity=0.356  Sum_probs=58.9  Template_Neff=4.800

Q d12asa_         235 VWNPVLEDAFELSSMGIRVDADTLKHQ-LALTGDED---RLEL-EWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       235 v~n~~l~~a~ElSSmGirVd~~~L~~Q-l~~~~~~~---r~~~-~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      -||-+++ -.||+|=+||.....+.++ ++..+-..   .... -|-.++.-| -|.--|-|+|--||+|+|+...+|.|
T Consensus       468 ~YDLVlN-G~EiGgGSiRih~~~~Q~~if~~lg~~~e~~~~~Fg~lL~Al~yG-~PPHgGiAlGlDRlvmll~g~~sIRD  545 (580)
T tr|B9Y6V6|B9Y6  468 AYDLVLN-GEEIGGGSIRIHDPELQEKVFEILGISPEEAEEKFGFYLDAFKYG-APPHGGFALGLDRLVMLLTGEDNIRD  545 (580)
Confidence            4555666 7899999999998888887 55555543   2222 334455668 89999999999999999999999999


Q d12asa_         310 VQA  312 (327)
Q Consensus       310 Vq~  312 (327)
                      |.|
T Consensus       546 VIA  548 (580)
T tr|B9Y6V6|B9Y6  546 VIA  548 (580)
Confidence            986


No 98
>tr|A6DS64|A6DS64_9BACT Lysyl-tRNA synthetase-related protein OS=Lentisphaera araneosa HTCC2155 GN=LNTAR_23684 PE=3 SV=1
Probab=87.61  E-value=0.53  Score=49.17  Aligned_cols=266  Identities=20%  Similarity=0.232  Sum_probs=140.0  Template_Neff=5.200

Q d12asa_           5 KQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS   84 (327)
Q Consensus         5 Tq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~   84 (327)
                      -.+.++.|+++|.++     |.++|..|+.+....     ....-.|  |.++   +  .=..-|-.    +++|+.=-.
T Consensus        10 Ra~~l~~iR~FF~~r-----~~lEVeTP~l~~~~~-----~d~~i~~--~~~~---~--~yL~TSPE----~~MKrLLaa   68 (292)
T tr|A6DS64|A6DS   10 RARILAAIRRFFAER-----GVLEVETPALSQAPV-----TDPHLDP--FQTK---G--LYLHTSPE----FAMKRLLAA   68 (292)
Confidence            345667777777654     889999999987651     1111111  2210   0  11111111    233333222


Q d12asa_          85 AGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHF  164 (327)
Q Consensus        85 ~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l~~~Lp~~I~F  164 (327)
                      ...-||.=-.+-|.+|  .+..|+--+-=+.|=++-.....-++...+.|+.+...+..-+          ..+.....-
T Consensus        69 G~~~Ifqi~kvFR~~E--~g~~HnpEFTMLEWYr~g~~~~~lm~e~~~L~~~~~~~~~~~~----------~~~~~~~~~  136 (292)
T tr|A6DS64|A6DS   69 GSGRIYQICKVFRNGE--RGRRHNPEFTMLEWYRPGFDYEQLMDETEDLLQQVLECAGRKR----------CDPFAPPER  136 (292)
Confidence            3356777778888876  9999999998888888776665555555555554444433221          134445566


Q d12asa_         165 VHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAF  244 (327)
Q Consensus       165 I~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~  244 (327)
                      +|.+|+-.+|=+++|.+-.  ..--..-+|.-.|-..|..+.|  ..--||.--.    .....++-|    ++..-..|
T Consensus       137 ~s~~eaF~~~~gid~~~~~--~~~~~~~l~~~~VEP~L~~~~p--~fl~dyPa~~----AaLAr~~~~----~~~vaeRF  204 (292)
T tr|A6DS64|A6DS  137 LSYQEAFQRYAGIDPLSAS--RDDLFDLLFSHKVEPHLGQDRP--TFLYDYPASQ----AALARISPE----DPRVAERF  204 (292)
Confidence            7777777777777766532  1111222333333333321110  0001111100    000001100    12233344


Q d12asa_         245 ELSSMGIRV--------DADT----LKHQLALTGD----EDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIG  308 (327)
Q Consensus       245 ElSSmGirV--------d~~~----L~~Ql~~~~~----~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIg  308 (327)
                      ||=..|+-.        |+..    +.+..+.+..    .-.....|-+++..| +|.+-|=.||--||+|+++.+.+|.
T Consensus       205 Ely~~G~ELANgy~ELtD~~eqr~RF~~~~~~R~~~g~~~~p~De~fl~al~~g-lP~csGvAlGiDRLlMl~~g~~~I~  283 (292)
T tr|A6DS64|A6DS  205 ELYICGIELANGFHELTDAAEQRRRFEADNAERARLGLEQYPIDEDFLAALEAG-MPDCSGVALGVDRLLMLALGAESID  283 (292)
Confidence            444444432        3333    2222222221    112334445555444 9999999999999999999999999


Q d12asa_         309 QVQAGVWP  316 (327)
Q Consensus       309 EVq~svW~  316 (327)
                      +|.+..|.
T Consensus       284 ~V~~F~~~  291 (292)
T tr|A6DS64|A6DS  284 DVIAFPVE  291 (292)
Confidence            99998775


No 99
>tr|A6MK47|A6MK47_CALJA Lysyl-tRNA synthetase-like protein (Fragment) OS=Callithrix jacchus PE=2 SV=1
Probab=87.29  E-value=0.85  Score=43.10  Aligned_cols=40  Identities=35%  Similarity=0.554  Sum_probs=34.2  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....|..++.-| ||.|.|=|.|..|+.|||-....|.||-
T Consensus       112 idenfctaleyg-lpptagwgmgidrvamfltdsnnikevl  151 (180)
T tr|A6MK47|A6MK  112 IDENFCTALEYG-LPPTAGWGMGIDRVTMFLTDSHNIKEVL  151 (180)
Confidence            445566666666 9999999999999999999999999985


No 100
>tr|C5LM59|C5LM59_PERM5 Lysyl-tRNA synthetase, putative OS=Perkinsus marinus (strain ATCC 50983 / TXsc) GN=Pmar_PMAR008972 PE=3 SV=1
Probab=87.17  E-value=0.88  Score=47.20  Aligned_cols=51  Identities=33%  Similarity=0.474  Sum_probs=39.3  Template_Neff=1.000

Q d12asa_         260 HQLALTGDED--RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       260 ~Ql~~~~~~~--r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .|.+..+..+  .....|-.++.-| ||.|=|=|+|-.||||||-..-.|.||-
T Consensus       316 aqakaagddeacdvdngfvtaleyg-lpptggwglgvdrlcmfltdninikevi  368 (388)
T tr|C5LM59|C5LM  316 AQAKAAGDDEACDVDNGFVTALEYG-LPPTGGWGLGVDRLCMFLTDNINIKEVI  368 (388)
Confidence            4445555444  2455666666666 9999999999999999999999999985


No 101
>tr|F2TYB2|F2TYB2_SALR5 Lysine--tRNA ligase OS=Salpingoeca rosetta (strain ATCC 50818 / BSB-021) GN=PTSG_01071 PE=3 SV=1
Probab=87.17  E-value=0.86  Score=51.12  Aligned_cols=59  Identities=27%  Similarity=0.476  Sum_probs=45.5  Template_Neff=1.257

Q d12asa_         252 RVDADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       252 rVd~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .|.++.+..|.+.....+    ...-.|-+++..| ||.|.|=|||..|+.|||-.-..|.||-
T Consensus       489 ~vqrerfeqq~~dka~gddeaq~ide~fvkalehg-lpptagwg~g~dr~~mfltdsnnikevl  551 (796)
T tr|F2TYB2|F2TY  489 QIQRERFEQQMQDKAMGDDEAQDIDEAFIKALEHG-LPPTAGWGMGIDRMTMFLTDSNNIKEVL  551 (796)
Confidence            345666666665433222    2556788888887 9999999999999999999999999985


No 102
>tr|A0A0M7R9W6|A0A0M7R9W6_ECOLX Lysyl-tRNA synthetase OS=Escherichia coli GN=BN1843_28000 PE=4 SV=1
Probab=87.12  E-value=0.89  Score=34.87  Aligned_cols=29  Identities=45%  Similarity=0.669  Sum_probs=22.9  Template_Neff=1.000

Q d12asa_         275 WHQALLRGEMPQTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       275 ~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k  304 (327)
                      |..++.-| ||.|.|=|.|..||.|+|-..
T Consensus         7 fctaleyg-lpptaglgmgidrlamlltds   35 (50)
T tr|A0A0M7R9W6|    7 FCTALEYG-LPPTAGLGMGIDRLAMLLTDS   35 (50)
Confidence            44444444 999999999999999998653


No 103
>tr|A0A0C3FBM0|A0A0C3FBM0_9HOMO Uncharacterized protein OS=Piloderma croceum F 1598 GN=PILCRDRAFT_91269 PE=3 SV=1
Probab=87.02  E-value=0.91  Score=50.86  Aligned_cols=71  Identities=37%  Similarity=0.515  Sum_probs=51.1  Template_Neff=1.000

Q d12asa_         244 FELSSMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAA  318 (327)
Q Consensus       244 ~ElSSmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~  318 (327)
                      -|+-|-|-|-- +-.|.+||+.+--+.....+|-.+..=| +|-| |||||.-||.|++|+--...  .+|.+|.+
T Consensus       739 eeilsggqrshsayelekqleeacinseemrdyvnafrwg-mpht-gggiglerlvmlflelgnlr--naslfphd  810 (895)
T tr|A0A0C3FBM0|  739 EEILSGGQRSHSAYELEKQLEEACINSEEMRDYVNAFRWG-MPHT-GGGIGLERLVMLFLELGNLR--NASLFPHD  810 (895)
Confidence            36677777754 4578999998766655556677666666 6765 99999999999999875543  35666543


No 104
>tr|J8USK7|J8USK7_PSEPU Uncharacterized protein OS=Pseudomonas putida S11 GN=PPS11_43328 PE=4 SV=1
Probab=86.81  E-value=0.96  Score=38.40  Aligned_cols=18  Identities=50%  Similarity=0.787  Sum_probs=15.6  Template_Neff=1.000

Q d12asa_         290 GGIGQSRLTMLLLQLPHI  307 (327)
Q Consensus       290 GGIgqSRl~M~lL~k~HI  307 (327)
                      .||||||++..||+...|
T Consensus         7 sgigqsriaavllradai   24 (88)
T tr|J8USK7|J8US    7 SGIGQSRIAAVLLRADAI   24 (88)
Confidence            489999999999997655


No 105
>tr|W7K2V9|W7K2V9_PLAFA Asparagine-tRNA ligase OS=Plasmodium falciparum UGT5.1 GN=C923_01075 PE=3 SV=1
Probab=86.47  E-value=1  Score=51.11  Aligned_cols=58  Identities=29%  Similarity=0.383  Sum_probs=51.6  Template_Neff=1.000

Q d12asa_         254 DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       254 d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +-..|++|.+..+.+-++-.||-+.-..|.+|- .|-|||-.||.|||-.-..|..|-.
T Consensus       991 nlknlkkqmenkklnmklydpyiqlrkhgniph-agfgigmdrlimflssinnirdvvt 1048 (1058)
T tr|W7K2V9|W7K2  991 NLKNLKKQMENKKLNMKLYDPYIQLRKHGNIPH-AGFGIGMDRLIMFLSSINNIRDVVT 1048 (1058)
Confidence            456799999999998899999999999999996 5999999999999999888888753


No 106
>tr|A0A0N5CP41|A0A0N5CP41_THECL Uncharacterized protein OS=Thelazia callipaeda PE=4 SV=1
Probab=86.39  E-value=1.1  Score=34.23  Aligned_cols=36  Identities=39%  Similarity=0.623  Sum_probs=28.4  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHI  307 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HI  307 (327)
                      ....|..++.-| ||.|.|=|.|..||.|.|-....|
T Consensus        12 idenfctaleyg-lpptagwgmgidrltmiltdsnni   47 (48)
T tr|A0A0N5CP41|   12 IDENFCTALEYG-LPPTAGWGMGIDRLTMILTDSNNI   47 (48)
Confidence            445666666666 999999999999999998765544


No 107
>tr|A0A0A9R9K1|A0A0A9R9K1_ARUDO OVA5 OS=Arundo donax PE=3 SV=1
Probab=86.28  E-value=1  Score=37.19  Aligned_cols=39  Identities=38%  Similarity=0.536  Sum_probs=32.0  Template_Neff=1.366

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      -+|-.++.-| +|++.|=|.|..||.|+|-.-+.|..|-|
T Consensus        25 edf~t~leyg-mp~a~gmglgi~rlvmll~~sasirdvia   63 (71)
T tr|A0A0A9R9K1|   25 EDFVTSLEYG-MPPASGMGLGINRLVMLLKNSASIRDVIA   63 (71)
Confidence            3444455445 89999999999999999999999998865


No 108
>tr|A0A098ED06|A0A098ED06_9ZZZZ Uncharacterized protein OS=groundwater metagenome GN=MSIBF_A3110004 PE=4 SV=1
Probab=85.94  E-value=1.2  Score=35.69  Aligned_cols=32  Identities=44%  Similarity=0.614  Sum_probs=23.0  Template_Neff=1.000

Q d12asa_         229 LNGDILVWN-PVLEDAFELSSMG-IRVDADTLKH  260 (327)
Q Consensus       229 LNGDilv~n-~~l~~a~ElSSmG-irVd~~~L~~  260 (327)
                      -|..|+--| +.+.+++|+|||| |.|||+-|+.
T Consensus        24 knteifntnreilsksiemssmgmidvdpellkt   57 (62)
T tr|A0A098ED06|   24 KNTEIFNTNREILSKSIEMSSMGMIDVDPELLKT   57 (62)
Confidence            344444333 4578999999999 7899988764


No 109
>tr|W7TMK7|W7TMK7_9STRA Lysyl-trna synthetase (Fragment) OS=Nannochloropsis gaditana GN=Naga_101432g2 PE=3 SV=1
Probab=85.77  E-value=1.2  Score=43.06  Aligned_cols=38  Identities=37%  Similarity=0.448  Sum_probs=31.7  Template_Neff=1.345

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .+|.-++.. -||.|=|=|-|..||.|||-.|..|.||-
T Consensus        23 eefcva~e~-~lpptggwgcgidr~amfl~~k~nikevl   60 (190)
T tr|W7TMK7|W7TM   23 EEFCVALEH-ALPPTGGWGCGIDRLAMFLADKHNIKEVL   60 (190)
Confidence            345555444 49999999999999999999999999984


No 110
>tr|T0YBN7|T0YBN7_9ZZZZ tRNA synthetase class II (D K and N) (Fragment) OS=mine drainage metagenome GN=B2A_13458 PE=4 SV=1
Probab=85.54  E-value=1.3  Score=41.89  Aligned_cols=134  Identities=25%  Similarity=0.300  Sum_probs=68.5  Template_Neff=1.000

Q d12asa_         156 PFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLS-DGHRHDVRAPDYDDWSTPSELGHAGLNGDIL  234 (327)
Q Consensus       156 ~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~-~G~~Hd~RapDYDDW~t~~~~~~~gLNGDil  234 (327)
                      +.+...-.-.++.||+.+|--    .-|.++.|.|.--|.-     |- --+-.|...||...--         ||-|++
T Consensus        18 prfgpawktytthelearygp----nwemaaskahdqpfwa-----lchrrefydredpdepghf---------lnydli   79 (175)
T tr|T0YBN7|T0YB   18 PRFGPAWKTYTTHELEARYGP----NWEMAASKAHDQPFWA-----LCHRREFYDREDPDEPGHF---------LNYDLI   79 (175)
Confidence            455555556677888888732    2345555555554431     11 1122333344433333         888876


Q d12asa_         235 VWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       235 v~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +=|.. +.|+.=...--+.|  .+...+...-.. -..+.-.-.|..| |-.+.|||+|.-||.-||-+..|+|.||.
T Consensus        80 ypngy-gealsggerewqpd--rirdrisrdpip-aaalaayleiahg-lvssaggglgierllrfltrsphvgdvqm  152 (175)
T tr|T0YBN7|T0YB   80 YPNGY-GEALSGGEREWQPD--RIRDRISRDPIP-AAALAAYLEIAHG-LVSSAGGGLGIERLLRFLTRSPHVGDVQM  152 (175)
Confidence            54432 33332222111111  111111000000 0111111223333 67789999999999999999999999995


No 111
>tr|T2JYP4|T2JYP4_CROWT Lysyl-tRNA synthetase (Class II) OS=Crocosphaera watsonii WH 0402 GN=CWATWH0402_5103 PE=3 SV=1
Probab=85.36  E-value=1.3  Score=39.36  Aligned_cols=37  Identities=49%  Similarity=0.574  Sum_probs=30.5  Template_Neff=1.000

Q d12asa_         275 WHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       275 ~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      |-.++.-| +|.|=|=|||..||.|+|-..+.|..|-+
T Consensus        10 fltaleyg-mpptgglgigidrlimlltdsasirdvia   46 (116)
T tr|T2JYP4|T2JY   10 FLTALEYG-MPPTGGLGIGIDRLIMLLTDSASIRDVIA   46 (116)
Confidence            33344444 89999999999999999999999998865


No 112
>tr|H9JJU4|H9JJU4_BOMMO Uncharacterized protein OS=Bombyx mori PE=3 SV=1
Probab=85.35  E-value=1.3  Score=36.89  Aligned_cols=30  Identities=40%  Similarity=0.802  Sum_probs=24.1  Template_Neff=1.000

Q d12asa_         286 QTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPA  317 (327)
Q Consensus       286 ~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~  317 (327)
                      -.+|||||.-|+.|+.|.-..|..  .|.+|.
T Consensus        42 hgvgggiglervvmlylgldnirk--tsmfpr   71 (78)
T tr|H9JJU4|H9JJ   42 HGVGGGIGLERVVMLYLGLDNIRK--TSMFPR   71 (78)
Confidence            358999999999999999888865  344544


No 113
>tr|S9QX59|S9QX59_9RHOB Asparagine synthetase A OS=Rubellimicrobium thermophilum DSM 16684 GN=ruthe_02357 PE=4 SV=1
Probab=84.79  E-value=1.5  Score=37.66  Aligned_cols=41  Identities=32%  Similarity=0.669  Sum_probs=35.8  Template_Neff=1.000

Q d12asa_         278 ALLRGEMPQTIGGGIGQSRLTMLLLQL--PHIGQVQAGVWPAA  318 (327)
Q Consensus       278 ~ll~~~lP~TIgGGIgqSRl~M~lL~k--~HIgEVq~svW~~~  318 (327)
                      ++--|..|.+.|||=|..+..-|||+.  .|..|--.|-||+.
T Consensus        46 avglgrvpsafgggegegqftsflleaqarhlaepggsgwpnr   88 (92)
T tr|S9QX59|S9QX   46 AVGLGRVPSAFGGGEGEGQFTSFLLEAQARHLAEPGGSGWPNR   88 (92)
Confidence            444578999999999999999999975  69999999999974


No 114
>tr|A0A0D2X2W4|A0A0D2X2W4_CAPO3 Lysine--tRNA ligase OS=Capsaspora owczarzaki (strain ATCC 30864) GN=CAOG_004076 PE=3 SV=1
Probab=84.58  E-value=1.5  Score=47.98  Aligned_cols=55  Identities=31%  Similarity=0.479  Sum_probs=43.4  Template_Neff=1.000

Q d12asa_         256 DTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       256 ~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..+.+|.+.....+.    ....|..++.-| ||.|+|=|.|..||+|+|-..+.|.||-
T Consensus       564 arfaeqakardagddeaqpmdsafcdaleyg-lpptvgwgmgldrlvmlltdsanirevi  622 (652)
T tr|A0A0D2X2W4|  564 ARFAEQAKARDAGDDEAQPMDSAFCDALEYG-LPPTVGWGMGLDRLVMLLTDSANIREVI  622 (652)
Confidence            345666666554442    456677777777 9999999999999999999999999985


No 115
>tr|A0A067QMI8|A0A067QMI8_9HOMO Uncharacterized protein OS=Jaapia argillacea MUCL 33604 GN=JAAARDRAFT_29893 PE=4 SV=1
Probab=84.41  E-value=1.6  Score=35.60  Aligned_cols=42  Identities=24%  Similarity=0.416  Sum_probs=34.7  Template_Neff=1.000

Q d12asa_          92 HMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKST  133 (327)
Q Consensus        92 dMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~t  133 (327)
                      .+..||-.|...|-...+-|||||..+||-..-|++...-+.
T Consensus        12 tfaqirleetpmdltstikvdqwdlsrvilanlrtlgwvsea   53 (68)
T tr|A0A067QMI8|   12 TFAQIRLEETPMDLTSTIKVDQWDLSRVILANLRTLGWVSEA   53 (68)
Confidence            456788888667888889999999999999999998776554


No 116
>tr|M1CK15|M1CK15_SOLTU Uncharacterized protein OS=Solanum tuberosum GN=PGSC0003DMG400026909 PE=4 SV=1
Probab=84.09  E-value=1.7  Score=35.09  Aligned_cols=38  Identities=39%  Similarity=0.548  Sum_probs=31.0  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      +...|..++.-| ||.|=|=|.|..||.|+|-..+.|.-
T Consensus        16 ldenfctaleyg-lpptggwglgidrlamlltdsqnikv   53 (64)
T tr|M1CK15|M1CK   16 LDENFCTALEYG-LPPTGGWGLGIDRLAMLLTDSQNIKV   53 (64)
Confidence            445666676666 99999999999999999998887753


No 117
>tr|A0A0G1BLE2|A0A0G1BLE2_9BACT Lysine-tRNA ligase (Fragment) OS=Microgenomates (Gottesmanbacteria) bacterium GW2011_GWC2_42_8 GN=UV46_C0051G0012 PE=4 SV=1
Probab=84.06  E-value=1.7  Score=33.99  Aligned_cols=43  Identities=28%  Similarity=0.410  Sum_probs=35.3  Template_Neff=1.000

Q d12asa_         268 EDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       268 ~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .-+-..+|--++.+| ||++-|-..|..|.+|++-....|..|-
T Consensus         7 pvkwdrefimalksg-lplssgvalgidrtamlfcdtknirdvl   49 (54)
T tr|A0A0G1BLE2|    7 PVKWDREFIMALKSG-LPLSSGVALGIDRTAMLFCDTKNIRDVL   49 (54)
Confidence            334455666677777 9999999999999999999999998874


No 118
>tr|A0A0E8TK72|A0A0E8TK72_STREE Sulfite exporter TauE/SafE OS=Streptococcus pneumoniae GN=ERS020526_05916 PE=4 SV=1
Probab=83.09  E-value=2  Score=35.26  Aligned_cols=29  Identities=45%  Similarity=0.683  Sum_probs=18.4  Template_Neff=1.000

Q d12asa_         283 EMPQTIGGGIGQSR---LTMLLLQLPHIGQVQ  311 (327)
Q Consensus       283 ~lP~TIgGGIgqSR---l~M~lL~k~HIgEVq  311 (327)
                      -||-+||||||-|-   |--.||-.--||..-
T Consensus         4 ilpiaigggigyssegyldyilllqvligtml   35 (70)
T tr|A0A0E8TK72|    4 ILPIAIGGGIGYSSEGYLDYILLLQVLIGTML   35 (70)
Confidence            37999999999985   333344334444443


No 119
>tr|R5QIY9|R5QIY9_9FIRM Lysine--tRNA ligase OS=Ruminococcus torques CAG:61 GN=BN734_01976 PE=3 SV=1
Probab=82.92  E-value=2  Score=45.93  Aligned_cols=55  Identities=35%  Similarity=0.515  Sum_probs=41.5  Template_Neff=1.572

Q d12asa_         256 DTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       256 ~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +.+..|++.+..-+    ....+|-.++.-| +|.|=|=|+|..||||||-+-+.|..|-
T Consensus       183 erf~~q~~~~~~g~eean~~d~df~~ale~g-mpp~ggig~gidrl~m~~t~s~sirdvl  241 (396)
T tr|R5QIY9|R5QI  183 ERFEAQMKLREQGDEEANMIDMDFVNALEYG-MPPTGGIGYGIDRLVMLLTNSASIRDVL  241 (396)
Confidence            44566666554333    2556677777767 8999999999999999999999988773


No 120
>tr|A0A0R3TSL9|A0A0R3TSL9_HYMNN Uncharacterized protein OS=Hymenolepis nana PE=4 SV=1
Probab=82.82  E-value=2.1  Score=40.62  Aligned_cols=42  Identities=36%  Similarity=0.545  Sum_probs=31.5  Template_Neff=1.000

Q d12asa_         285 PQTIGGGIGQSRLTMLLLQLPHIG-QVQAGVWPAAVRESVPSL  326 (327)
Q Consensus       285 P~TIgGGIgqSRl~M~lL~k~HIg-EVq~svW~~~~~~~~~~~  326 (327)
                      -.|-.||.|.|-|..----++.|| +-|+|.||++.+..|..+
T Consensus        49 vvtksgglgtseltynqhlrqkigpqeqvsywpederrrcpai   91 (176)
T tr|A0A0R3TSL9|   49 VVTKSGGLGTSELTYNQHLRQKIGPQEQVSYWPEDERRRCPAI   91 (176)
Confidence            356789999999876433334454 569999999999999764


No 121
>tr|D4JFP2|D4JFP2_9FIRM Lysine--tRNA ligase OS=Faecalitalea cylindroides T2-87 GN=lysS PE=3 SV=1
Probab=82.64  E-value=1.8  Score=49.98  Aligned_cols=201  Identities=21%  Similarity=0.220  Sum_probs=105.9  Template_Neff=3.100

Q d12asa_          93 MKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQIHFVHSQELL  171 (327)
Q Consensus        93 MnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~FI~sqeL~  171 (327)
                      |+.+=|.| -+|..||--+-.+-+=+--...+--++.-.+.++.+-.....|...   .|. ..-.|...-..|+--+-.
T Consensus       247 IgR~FRNE-GmD~~HnPEFT~~E~YqAY~Dy~~mm~lte~li~~~A~~v~Gt~~i---~~~G~eiDl~~pw~ritm~dav  322 (649)
T tr|D4JFP2|D4JF  247 IGRVFRNE-GMDSTHNPEFTMLELYQAYGDYNDMMDLTEELIQEVADEVLGTTQI---TYDGTEIDLDGPWKRITMYDAL  322 (649)
Confidence            44444566 5888888655555444444444444444444444444433332111   232 222333333344444444


Q d12asa_         172 SRYPD-----LDAKGRERAIAKDLGAVF---LVGIGGKLSDG----------------------HRHDVRAPDYDDWSTP  221 (327)
Q Consensus       172 ~~YP~-----LtpkeRE~~i~ke~gAvF---i~gIG~~L~~G----------------------~~Hd~RapDYDDW~t~  221 (327)
                      +.|-.     -++.+..+++|+++|.-+   -.+.|+.+.-=                      +|---|.|+.......
T Consensus       323 ke~~G~d~~~~~~~eea~~~A~~~~ie~~~~~~~~Gkii~e~Fee~vE~~LiqPTFV~d~PveiSPLaK~~~~dP~~ter  402 (649)
T tr|D4JFP2|D4JF  323 SEALGEDITPITTDEELRAIAKEHGIEVEEKGWGKGKIIEELFEEFVEDKLIQPTFVTDYPVEISPLAKRHRSDPGLTER  402 (649)
Confidence            44332     233455788999998887   35555544210                      1111222222222100


Q d12asa_         222 SEL---GHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQ  294 (327)
Q Consensus       222 ~~~---~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgq  294 (327)
                      -+-   |.---|+=--.+||+-             -++.+.+|.+.+..-+    ....+|-.++.-| +|.|=|=|||.
T Consensus       403 fElfI~G~E~anAfSELnDPid-------------QreRF~~Q~~~~a~GDdEA~~~DeDFl~ALE~G-mPPtGGlG~GI  468 (649)
T tr|D4JFP2|D4JF  403 FELFINGWELANAFSELNDPID-------------QRERFEAQARLAAAGDDEAMQLDEDFLRALEYG-MPPTGGLGIGI  468 (649)
Confidence            000   0000111111223322             2456777777665444    3667888888888 99999999999


Q d12asa_         295 SRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       295 SRl~M~lL~k~HIgEVq  311 (327)
                      .||+|+|-.-.-|.+|-
T Consensus       469 DRlvMllT~~~sIRdvl  485 (649)
T tr|D4JFP2|D4JF  469 DRLVMLLTDSASIRDVI  485 (649)
Confidence            99999999999998874


No 122
>tr|A0A017T1J2|A0A017T1J2_9DELT Translation elongation factor P Lys34:lysine transferase OS=Chondromyces apiculatus DSM 436 GN=CAP_6182 PE=3 SV=1
Probab=82.26  E-value=2  Score=47.09  Aligned_cols=278  Identities=22%  Similarity=0.304  Sum_probs=146.1  Template_Neff=2.400

Q d12asa_           5 KQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS   84 (327)
Q Consensus         5 Tq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~   84 (327)
                      -+.+...++..|-.     =+..+|..|+.|.. .|+.=.|+-..--=.|-|.+          -=..-||+.-.  || 
T Consensus       116 Ra~~l~avR~~f~~-----~gFlEveTP~lv~~-~~~d~HLd~~~~~~~~LiTS----------~EyQmKRll~g--G~-  176 (431)
T tr|A0A017T1J2|  116 RAAILRAVRGYFAE-----RGFLEVETPLLVPS-PGLDLHLDPFQAGPGWLITS----------PEYQMKRLLSG--GF-  176 (431)
Confidence            34455566666654     37889999999865 46655554322111222222          11234555543  33 


Q d12asa_          85 AGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLP-DQI  162 (327)
Q Consensus        85 ~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp-~~I  162 (327)
                        +-+||=-...|++|  +...|-.-+.-+.|-+-...-|--+.--...|...-.++..+-..+.-.=. +  .|. ..-
T Consensus       177 --ervftl~~~fR~ge--~G~~HnpEFTMlEW~RA~a~~d~i~~DtEalV~~~~~al~~~~~~~~~~~r~i--~Ls~~~w  250 (431)
T tr|A0A017T1J2|  177 --ERVFTLCRCFRRGE--LGARHNPEFTMLEWYRAFAGLDDIMNDTEALVARAAEALGPGARCISYRGRQI--RLSGPPW  250 (431)
Confidence              67899999999998  889999999999998888877766666677777777777665554321111 1  111 111


Q d12asa_         163 HFVHSQELLSRYPDLDAKGRE-----RAIAKDLGAVFLVGIGGKL-----------SDG-HRHDVRAP--DYDDWSTPSE  223 (327)
Q Consensus       163 ~FI~sqeL~~~YP~LtpkeRE-----~~i~ke~gAvFi~gIG~~L-----------~~G-~~Hd~Rap--DYDDW~t~~~  223 (327)
                      .=+|-.|-..+|-+..-..-+     .+.|..-+    ++||-..           -|. +||=++.-  =--||-.+-.
T Consensus       251 ~Rltv~eA~aR~~~V~~~~~~s~~~~~~~a~~~~----~~i~~~~~~d~~~~~~ll~d~~~P~l~~~~PVfL~dwPA~ma  326 (431)
T tr|A0A017T1J2|  251 ERLTVREAFARYAGVAVDDAFSAASLAQAARRAD----LGIELAFRSDEERFFTLLVDQVEPALGAERPVFLVDWPAPMA  326 (431)
Confidence            112223333333221111000     00000000    1222220           010 23333311  1235543322


Q d12asa_         224 LGHAGLNGDILVWNPVLEDAFELSSMGIRV--------DAD----TLKHQLALTGDED----RLELEWHQALLRGEMPQT  287 (327)
Q Consensus       224 ~~~~gLNGDilv~n~~l~~a~ElSSmGirV--------d~~----~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~T  287 (327)
                      ..-.-.-|     .|.+-..+|+=-.||-.        |+.    -+..-+......-    -+...|-+++.+| +|.+
T Consensus       327 s~A~l~P~-----~p~l~eR~E~fvagiEl~ngF~eL~D~~eQ~a~~~~~~a~R~~~G~p~~~lD~r~l~AL~eG-~Pp~  400 (431)
T tr|A0A017T1J2|  327 SLARLKPG-----SPALAERFELFVAGIELCNGFGELTDPAEQRARFAADQAARAARGLPVYPLDERFLAALEEG-MPPG  400 (431)
Confidence            22111222     13333344443333332        221    1112222222211    2556788888888 9999


Q d12asa_         288 IGGGIGQSRLTMLLLQLPHIGQVQAGVWPA  317 (327)
Q Consensus       288 IgGGIgqSRl~M~lL~k~HIgEVq~svW~~  317 (327)
                      .|--.|-.||.|+|+...||++|-+=.|+.
T Consensus       401 aG~ALGvDRLvmL~~Ga~~I~~vl~~s~~~  430 (431)
T tr|A0A017T1J2|  401 AGNALGVDRLVMLLTGAEHIRDVLAFSWDE  430 (431)
Confidence            999999999999999999999998877753


No 123
>tr|H2KSF2|H2KSF2_CLOSI Lysine--tRNA ligase OS=Clonorchis sinensis GN=CLF_109131 PE=3 SV=1
Probab=82.19  E-value=2.2  Score=47.43  Aligned_cols=58  Identities=28%  Similarity=0.383  Sum_probs=44.9  Template_Neff=1.314

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+|.+...+.+-    ..-+|-.++--| ||.|.|=|.|..||.|||-....|.||-
T Consensus       462 vqr~rf~eqak~kaagd~eam~tde~fltalgyg-lpptagwg~gidrl~mfltntnn~kevi  523 (647)
T tr|H2KSF2|H2KS  462 VQRQRFVEQAKAKAAGDGEAMPTDEPFLTALGYG-LPPTAGWGMGIDRLAMFLTNTNNLKEVI  523 (647)
Confidence            667777888765544332    445566666555 9999999999999999999999999985


No 124
>tr|D8LM77|D8LM77_ECTSI Lysine--tRNA ligase OS=Ectocarpus siliculosus GN=Esi_0392_0007 PE=3 SV=1
Probab=81.79  E-value=2.5  Score=47.64  Aligned_cols=139  Identities=24%  Similarity=0.267  Sum_probs=78.8  Template_Neff=1.000

Q d12asa_         170 LLSRYPDLDAKG---RERAIAKDLGAVF---------L-VGIGGKLSDGHRHDVRAPDYDDWSTPSEL---GHAGLNGDI  233 (327)
Q Consensus       170 L~~~YP~Ltpke---RE~~i~ke~gAvF---------i-~gIG~~L~~G~~Hd~RapDYDDW~t~~~~---~~~gLNGDi  233 (327)
                      |.-..|.|+..|   ..++..-|||.-+         + --+|.-|.|.-.|-..--|-..-.+|-..   ...||.-..
T Consensus       390 lkvklpalddpeidgklqalltehglecapplttarlldklvgdfledkcvhptfitdhpeimsplakyhrsrpglterf  469 (849)
T tr|D8LM77|D8LM  390 LKVKLPALDDPEIDGKLQALLTEHGLECAPPLTTARLLDKLVGDFLEDKCVHPTFITDHPEIMSPLAKYHRSRPGLTERF  469 (849)
Confidence            344566666544   3445555666322         1 23577777777887776665554433322   233333221


Q d12asa_         234 LVWNPVLEDAFELSSMG---IRVDADTLKHQLALTGDED----RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPH  306 (327)
Q Consensus       234 lv~n~~l~~a~ElSSmG---irVd~~~L~~Ql~~~~~~~----r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~H  306 (327)
                      -+  .++++.++-++--   --|.+..+..|.+.+-..+    -...+|..++.-| ||.|.|=|+|-.||.|||-.|..
T Consensus       470 el--fvcgrelcnaytelnnpmvqrqrfldqskasvagddeaqvhdeefctameyg-lpptagwgvgvdrltmflsdknn  546 (849)
T tr|D8LM77|D8LM  470 EL--FVCGRELCNAYTELNNPMVQRQRFLDQSKASVAGDDEAQVHDEEFCTAMEYG-LPPTAGWGVGVDRLTMFLSDKNN  546 (849)
Confidence            11  1122222222211   1245556666665543322    1334566665555 99999999999999999999999


Q d12asa_         307 IGQVQ  311 (327)
Q Consensus       307 IgEVq  311 (327)
                      |.||-
T Consensus       547 ikevl  551 (849)
T tr|D8LM77|D8LM  547 IKEVL  551 (849)
Confidence            99985


No 125
>tr|M7SLK5|M7SLK5_EUTLA Putative aspartyl-trna cytoplasmic protein OS=Eutypa lata (strain UCR-EL1) GN=UCREL1_7942 PE=3 SV=1
Probab=81.57  E-value=2.1  Score=51.20  Aligned_cols=254  Identities=25%  Similarity=0.328  Sum_probs=148.5  Template_Neff=3.200

Q d12asa_          24 LGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRL  103 (327)
Q Consensus        24 LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~l  103 (327)
                      .+.++|-.|=+-...|     =+|   +=-|.|.-- +..+=..||--=-|.|++.- +|.    =+=..-++=|-|  =
T Consensus       271 ~gFiEIHTPKLQg~AT-----ESG---AsVFkV~YF-GR~AFLAQSPQLaKQM~IaA-Df~----RVyEIGpVFRAE--N  334 (914)
T tr|M7SLK5|M7SL  271 QGFIEIHTPKLQGGAT-----ESG---ASVFKVNYF-GRPAFLAQSPQLAKQMAIAA-DFE----RVYEIGPVFRAE--N  334 (914)
Confidence            4677777775543332     123   333777666 77888888866666666642 221    111112223333  3


Q d12asa_         104 SPLH--SVYVDQWDWERVMGDGE-RQFSTLKSTVEAIWAGIKA----TEAAVSEEFG-LAPFLPDQIHFVHSQE---LL-  171 (327)
Q Consensus       104 d~~H--SiyVDQWDWEkvI~~~d-Rnl~~Lk~tV~kIy~al~~----te~~v~~~y~-l~~~Lp~~I~FI~sqe---L~-  171 (327)
                      +|+|  =--+-.+|-|..|.+.- -.|+.+-.+.+.||+.+++    --+.|...|| -.-..+++-.-|+..|   |+ 
T Consensus       335 SNThRHLTEYTGLDlEMaI~~hYHE~l~viD~~lK~IFk~iy~r~r~Eie~vk~~fPheDlvwldeTpii~F~egIkmL~  414 (914)
T tr|M7SLK5|M7SL  335 SNTHRHLTEYTGLDLEMAIEEHYHEALDVIDGMLKNIFKGIYERYRREIEIVKHQFPHEDLVWLDETPIIPFAEGIKMLN  414 (914)
Confidence            4554  44567899999998543 3455555555555555544    4445566676 2212222222222222   11 


Q d12asa_         172 ------------SRYPDLDA--KGRERAIAKD-LGAVFLVGIGGKLSDGHRHDVRA------PDYDDWSTPSELGHAGLN  230 (327)
Q Consensus       172 ------------~~YP~Ltp--keRE~~i~ke-~gAvFi~gIG~~L~~G~~Hd~Ra------pDYDDW~t~~~~~~~gLN  230 (327)
                                  ..+-||+-  .-|.=+++|| ||.=+.|      -|.-|-..|.      |+.+.|+         -.
T Consensus       415 dsGW~~edG~~~~e~EDL~Tr~EiRLGeLVKEKY~TDYYI------LDKFP~sARPFYTmpDp~d~~~T---------NS  479 (914)
T tr|M7SLK5|M7SL  415 DSGWRDEDGSPLSEDEDLSTRDEIRLGELVKEKYGTDYYI------LDKFPASARPFYTMPDPEDPRFT---------NS  479 (914)
Confidence                        11222321  2244445544 4443322      1333444443      4445555         34


Q d12asa_         231 GDILVWNPVLEDAFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       231 GDilv~n~~l~~a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      =||+|=      -=||-|=|=|+ |+..|.+.++..+-+...-.+|-+.-.=| -|.-.|||||-=||.|++|+--.|. 
T Consensus       480 FDIFvR------GQEI~SGGQRIHd~~~Le~rM~~~gi~p~~M~EYlegFrwG-aPPHaG~GIGLERlvmL~L~LGnIR-  551 (914)
T tr|M7SLK5|M7SL  480 FDIFVR------GQEILSGGQRIHDAKMLEERMRKAGIDPDTMEEYLEGFRWG-APPHAGAGIGLERLVMLFLKLGNIR-  551 (914)
Confidence            466653      45899999998 67888899999999888888999988888 6778899999999999999986554 


Q d12asa_         310 VQAGVWPA  317 (327)
Q Consensus       310 Vq~svW~~  317 (327)
                       -||.+|.
T Consensus       552 -~ASLFpR  558 (914)
T tr|M7SLK5|M7SL  552 -LASLFPR  558 (914)
Confidence             3555553


No 126
>tr|A0A0F6TGI3|A0A0F6TGI3_9BETA A38.5 OS=Rat cytomegalovirus ALL-03 GN=a38.5 PE=4 SV=1
Probab=81.04  E-value=2.7  Score=38.80  Aligned_cols=37  Identities=32%  Similarity=0.536  Sum_probs=30.9  Template_Neff=1.000

Q d12asa_         204 DGHRHDVRAPDYDDWSTPSELGHAGL-NGDILVWNPVL  240 (327)
Q Consensus       204 ~G~~Hd~RapDYDDW~t~~~~~~~gL-NGDilv~n~~l  240 (327)
                      .-+.||.|.-|...+..+...--.|| ||||++|+...
T Consensus        95 krqehdvrerdgqgfrreggevdtglgngdiyfwdgss  132 (144)
T tr|A0A0F6TGI3|   95 KRQEHDVRERDGQGFRREGGEVDTGLGNGDIYFWDGSS  132 (144)
Confidence            34789999999999987777767777 99999999754


No 127
>tr|W4F9U4|W4F9U4_9STRA Uncharacterized protein OS=Aphanomyces astaci GN=H257_19384 PE=4 SV=1
Probab=80.73  E-value=2.9  Score=40.33  Aligned_cols=58  Identities=28%  Similarity=0.383  Sum_probs=42.5  Template_Neff=1.000

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+|.+.+...+.    ....|..++.-| ||.|=|=|-|-.|+.|||-.|-.|.||-
T Consensus        45 vqrerfteqakqaadgddeaqphdeafctameyg-lpptggwgcgvdriamfltnkfnikevl  106 (194)
T tr|W4F9U4|W4F9   45 VQRERFTEQAKQAADGDDEAQPHDEAFCTAMEYG-LPPTGGWGCGVDRIAMFLTNKFNIKEVL  106 (194)
Confidence            456667777766544333    223344454444 9999999999999999999999999985


No 128
>tr|A0A0N8EBA7|A0A0N8EBA7_9CRUS Putative Aspartate--tRNA ligase, cytoplasmic OS=Daphnia magna PE=4 SV=1
Probab=80.53  E-value=3  Score=42.38  Aligned_cols=70  Identities=27%  Similarity=0.439  Sum_probs=43.9  Template_Neff=1.000

Q d12asa_         245 ELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPA  317 (327)
Q Consensus       245 ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~  317 (327)
                      |+-|-.-|+ |++-|.+....-+-.-..-..|-....-| -|.-.|||||--|++|++|.-..|..  .|.+|.
T Consensus       215 eimsgaqrihdpeflserarhhnidlskikayidsfryg-cpphagggigmervamlflgldnirk--tsmfpr  285 (292)
T tr|A0A0N8EBA7|  215 EIMSGAQRIHDPEFLSERARHHNIDLSKIKAYIDSFRYG-CPPHAGGGIGMERVAMLFLGLDNIRK--TSMFPR  285 (292)
Confidence            444444444 55555554443333323334455555555 68889999999999999999988875  344544


No 129
>tr|G5ANI0|G5ANI0_HETGA Lysine--tRNA ligase OS=Heterocephalus glaber GN=GW7_04322 PE=3 SV=1
Probab=80.31  E-value=2.8  Score=45.95  Aligned_cols=40  Identities=35%  Similarity=0.475  Sum_probs=35.2  Template_Neff=1.976

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....|..++.-| ||.|.|=|.|..|+.|||..-..|.||-
T Consensus       379 ~d~~fc~aleyg-lppt~gwg~gid~~~mfl~dsn~i~evl  418 (447)
T tr|G5ANI0|G5AN  379 VDETFCNALEYG-LPPTGGWGCGIDRLAMFLTDSNTIREVL  418 (447)
Confidence            455677777777 9999999999999999999999999985


No 130
>tr|A0A0G4LE24|A0A0G4LE24_9PEZI Uncharacterized protein (Fragment) OS=Verticillium longisporum GN=BN1723_012086 PE=3 SV=1
Probab=79.78  E-value=3.3  Score=42.80  Aligned_cols=42  Identities=29%  Similarity=0.362  Sum_probs=34.1  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVW  315 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW  315 (327)
                      ..|.-.|.- -||.|.|=|.|..|++|||-..-.|.||-+-.+
T Consensus       266 enfcmsley-glpptagwgmgidrmvmfltdhytirevlafpf  307 (337)
T tr|A0A0G4LE24|  266 ENFCMSLEY-GLPPTAGWGMGIDRMVMFLTDHYTIREVLAFPF  307 (337)
Confidence            344444444 499999999999999999999999999987554


No 131
>tr|V4KXH3|V4KXH3_EUTSA Uncharacterized protein OS=Eutrema salsugineum GN=EUTSA_v10015942mg PE=3 SV=1
Probab=79.57  E-value=2.8  Score=48.27  Aligned_cols=258  Identities=20%  Similarity=0.275  Sum_probs=155.3  Template_Neff=2.900

Q d12asa_          24 LGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRL  103 (327)
Q Consensus        24 LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~l  103 (327)
                      -+.++|-.|=++...|      .|.  .--|.+.-- +..+=..||--=.|.+|+- -+|..   .++ .-..=|-|  -
T Consensus       334 ~~FiEIHTPKLiggaS------EGG--a~VFkl~Yf-gq~AcLAQSPQL~KQMAIc-gDf~R---VFE-IGPVFRAE--d  397 (628)
T tr|V4KXH3|V4KX  334 EGFIEIHTPKLIGGSS------EGG--AAVFKLDYF-GQPACLAQSPQLYKQMAIC-GDFER---VFE-IGPVFRAE--D  397 (628)
Confidence            3678888887776653      121  223666654 5577788998777877753 12211   111 11112223  3


Q d12asa_         104 SPL--------HSVYVDQWDWERVMGDGERQ-FSTLKSTVEAIWAGIKATE----AAVSEEFGL--APFLP--DQIHFVH  166 (327)
Q Consensus       104 d~~--------HSiyVDQWDWEkvI~~~dRn-l~~Lk~tV~kIy~al~~te----~~v~~~y~l--~~~Lp--~~I~FI~  166 (327)
                      +|+        |=.-+-.+|.|..|.+.-.. ++.+-+...-||+.|.+--    +.+.+.||.  -..|+  ..|+|=.
T Consensus       398 SnTHR~~~~~~HLcEFvGLD~EM~ik~hY~evvD~vd~lFv~IF~~L~~~~~kEle~v~~QyP~eplk~l~~Tlrltf~E  477 (628)
T tr|V4KXH3|V4KX  398 SNTHRHPMKHGHLCEFVGLDLEMAIKEHYHEVVDVVDRLFVFIFDGLNKNCKKELEAVKRQYPFEPFKYLEKTLRLTFEE  477 (628)
Confidence            444        44667799999999876542 5566666667777665432    345667773  11333  2344433


Q d12asa_         167 SQELLSR-------YPDLD-AKGRE-RAI-AKDLGAVFLVGIGGKLSDGHRHDVRAPDYDD--WSTPSELGHAGLNGDIL  234 (327)
Q Consensus       167 sqeL~~~-------YP~Lt-pkeRE-~~i-~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDD--W~t~~~~~~~gLNGDil  234 (327)
                      .-+|++.       +-||+ ..||. -++ -.+|+.-|.|-.-.+|+ -.|- ---||+||  |+         --=|.+
T Consensus       478 gi~mLkeagv~i~~~~DLsTE~Er~LG~LVkeKY~TDFYil~ryPla-vRPF-YTMP~~~d~~yS---------NSfD~F  546 (628)
T tr|V4KXH3|V4KX  478 GIQMLREAGVEIDDLGDLSTEMEKKLGKLVKEKYNTDFYILHRYPLA-VRPF-YTMPCPDDPRYS---------NSYDFF  546 (628)
Confidence            3333331       33443 33332 122 33577777776555552 1111 11366666  44         111222


Q d12asa_         235 VWNPVLEDAFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAG  313 (327)
Q Consensus       235 v~n~~l~~a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~s  313 (327)
                            =+--||.|=+=|| |++-|.++.+..+-.-..-..|-.+..-|. |.--|+|||.=|+.|++|.-..|..+  |
T Consensus       547 ------mRGEEIiSGAQRIHdp~lL~era~~~gId~~ti~~YidsFryGa-pPHgG~GiGLERVvMLflgL~NIRkt--S  617 (628)
T tr|V4KXH3|V4KX  547 ------MRGEEILSGAQRIHDPELLEERAKECGIDVKTISTYIDSFRYGA-PPHGGFGIGLERVVMLFLGLNNIRKT--S  617 (628)
Confidence                  2345899999998 788888899999988788889999999995 55789999999999999999888764  4


Q d12asa_         314 VWPA  317 (327)
Q Consensus       314 vW~~  317 (327)
                      ..|.
T Consensus       618 LFPR  621 (628)
T tr|V4KXH3|V4KX  618 LFPR  621 (628)
Confidence            4443


No 132
>tr|A0A0F7UAI7|A0A0F7UAI7_NEOCL Lysyl-tRNA synthetase, related OS=Neospora caninum (strain Liverpool) GN=BN1204_012465 PE=3 SV=1
Probab=79.53  E-value=3.3  Score=48.51  Aligned_cols=41  Identities=41%  Similarity=0.559  Sum_probs=36.7  Template_Neff=1.414

Q d12asa_         270 RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       270 r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      --+.+|-.+++-| ||.|.|=|||..||+|+|-.-.-|..|-
T Consensus      1016 vpnqef~~almlg-lpptaglgigldrlvmlltd~~~irdvv 1056 (1065)
T tr|A0A0F7UAI7| 1016 VPNQEFVSALMLG-LPPTAGLGIGLDRLVMLLTDSAFIRDVV 1056 (1065)
Confidence            3678899999988 9999999999999999999998888774


No 133
>tr|A0A089QLH2|A0A089QLH2_9PROC Uncharacterized protein OS=Prochlorococcus sp. MIT 0604 GN=EW14_1407 PE=4 SV=1
Probab=79.38  E-value=3.4  Score=33.40  Aligned_cols=22  Identities=36%  Similarity=0.509  Sum_probs=18.6  Template_Neff=1.000

Q d12asa_         166 HSQELLSRYPDLDAKGRERAIA  187 (327)
Q Consensus       166 ~sqeL~~~YP~LtpkeRE~~i~  187 (327)
                      +++||+++||+|.|=+|...-+
T Consensus         6 tseeleemypdlepwqrdelev   27 (63)
T tr|A0A089QLH2|    6 TSEELEEMYPDLEPWQRDELEV   27 (63)
Confidence            6899999999999999865433


No 134
>tr|A0A060SS49|A0A060SS49_PYCCI Lysine--tRNA ligase OS=Pycnoporus cinnabarinus GN=BN946_scf184497.g4 PE=3 SV=1
Probab=78.99  E-value=3.6  Score=48.32  Aligned_cols=40  Identities=38%  Similarity=0.557  Sum_probs=34.8  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....|-.++..| ||.|=|=|||..||.|||-.-..|.||-
T Consensus       552 idetfvdalehg-lpptggwgigidrlvmfltdstnikevl  591 (1431)
T tr|A0A060SS49|  552 IDETFVDALEHG-LPPTGGWGIGIDRLVMFLTDSTNIKEVL  591 (1431)
Confidence            455677777777 9999999999999999999999999984


No 135
>tr|A0A0G4MVL4|A0A0G4MVL4_9PEZI Uncharacterized protein (Fragment) OS=Verticillium longisporum GN=BN1708_016604 PE=4 SV=1
Probab=78.94  E-value=3.4  Score=35.97  Aligned_cols=41  Identities=34%  Similarity=0.493  Sum_probs=34.6  Template_Neff=1.614

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +...|-+++..| +|.|-|=|.|-.||+|+|-.-..|..|-+
T Consensus        42 ld~~~~~a~e~g-~pp~~g~g~g~drlvml~~~~~~i~d~~~   82 (88)
T tr|A0A0G4MVL4|   42 LDKSYIKAMEYG-MPPTSGWGCGIDRLVMLFSDSNRISDCIT   82 (88)
Confidence            455677888888 99999999999999999998888877643


No 136
>tr|D1NED1|D1NED1_HAEIF Lysine--tRNA ligase OS=Haemophilus influenzae HK1212 GN=lysU PE=3 SV=1
Probab=78.17  E-value=4  Score=39.15  Aligned_cols=40  Identities=38%  Similarity=0.604  Sum_probs=32.5  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ...+|--++..| ||.|.|-|.|..||.|+.-....|..|-
T Consensus       136 kdedfvvalehg-lpptageglgidrlamlyanapsirdvi  175 (183)
T tr|D1NED1|D1NE  136 KDEDFVVALEHG-LPPTAGEGLGIDRLAMLYANAPSIRDVI  175 (183)
Confidence            344566666666 9999999999999999998888887764


No 137
>tr|A0A078F524|A0A078F524_BRANA BnaA05g27370D protein OS=Brassica napus GN=BnaA05g27370D PE=3 SV=1
Probab=78.09  E-value=4  Score=35.08  Aligned_cols=38  Identities=37%  Similarity=0.464  Sum_probs=31.3  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..|..++.-| |..|=|=|+|..||.|+|-..+.|.||-
T Consensus        27 etfcnaleyg-laptggwglgidrlamlltdsqnikevi   64 (89)
T tr|A0A078F524|   27 ETFCNALEYG-LAPTGGWGLGIDRLAMLLTDSQNIKEVI   64 (89)
Confidence            3456666656 6678888999999999999999999985


No 138
>tr|A0A086M299|A0A086M299_TOXGO tRNA ligases class II (D, K and N) domain-containing protein (Fragment) OS=Toxoplasma gondii RUB GN=TGRUB_220350B PE=3 SV=1
Probab=77.50  E-value=4.3  Score=46.02  Aligned_cols=40  Identities=43%  Similarity=0.590  Sum_probs=35.7  Template_Neff=1.096

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .+.+|-.+++-| ||.|.|=|||..||+|+|-.-+.|..|-
T Consensus       770 pnqefvsalmlg-lpptaglgigldrlvmlltdsafirdvv  809 (818)
T tr|A0A086M299|  770 PNQEFVSALMLG-LPPTAGLGIGLDRLVMLLTDSAFIRDVV  809 (818)
Confidence            566788888888 9999999999999999999999998874


No 139
>tr|E9ICF4|E9ICF4_SOLIN Putative uncharacterized protein (Fragment) OS=Solenopsis invicta GN=SINV_00256 PE=3 SV=1
Probab=77.34  E-value=4.4  Score=45.04  Aligned_cols=66  Identities=33%  Similarity=0.456  Sum_probs=48.0  Template_Neff=1.000

Q d12asa_         245 ELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       245 ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |+-|-.-|+ |++-|.+..+.-+-....-..|-.+..-| -|...|||||--|+.|+.|.-..|..|.
T Consensus       443 eiisgaqrihdpdflterakhhgidiekiksyidafryg-cpphagggigmervvmlylgldnirkvs  509 (686)
T tr|E9ICF4|E9IC  443 EIISGAQRIHDPDFLTERAKHHGIDIEKIKSYIDAFRYG-CPPHAGGGIGMERVVMLYLGLDNIRKVS  509 (686)
Confidence            444444443 66777766666666555555666666666 6888999999999999999999998764


No 140
>tr|W7W9Z0|W7W9Z0_9BURK Lysine--tRNA ligase OS=Methylibium sp. T29 GN=acoD PE=3 SV=1
Probab=76.75  E-value=4.7  Score=46.69  Aligned_cols=40  Identities=43%  Similarity=0.676  Sum_probs=34.7  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ...+|-+++..| ||.|=|.|||..||.|+|-....|..|-
T Consensus       460 fdadyiralehg-lpptggcgigidrlmmlltdspnirdvi  499 (1129)
T tr|W7W9Z0|W7W9  460 FDADYIRALEHG-LPPTGGCGIGIDRLMMLLTDSPNIRDVI  499 (1129)
Confidence            456788888777 9999999999999999999988888774


No 141
>tr|W4FAF2|W4FAF2_9STRA Uncharacterized protein OS=Aphanomyces astaci GN=H257_19384 PE=3 SV=1
Probab=76.67  E-value=4.8  Score=40.50  Aligned_cols=58  Identities=28%  Similarity=0.383  Sum_probs=42.7  Template_Neff=1.000

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.++.+.+|.+.+...+.    ....|..++.-| ||.|=|=|-|-.|+.|||-.|-.|.||-
T Consensus       112 vqrerfteqakqaadgddeaqphdeafctameyg-lpptggwgcgvdriamfltnkfnikevl  173 (261)
T tr|W4FAF2|W4FA  112 VQRERFTEQAKQAADGDDEAQPHDEAFCTAMEYG-LPPTGGWGCGVDRIAMFLTNKFNIKEVL  173 (261)
Confidence            456667777766544433    223344454444 9999999999999999999999999985


No 142
>tr|L1LCI1|L1LCI1_BABEQ Lysine--tRNA ligase OS=Babesia equi strain WA GN=BEWA_014220 PE=3 SV=1
Probab=76.38  E-value=4.9  Score=46.95  Aligned_cols=38  Identities=37%  Similarity=0.538  Sum_probs=32.6  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ..|..++.-| ||.|.|=|.|..||.||+-.|..|.||.
T Consensus      1206 esfctaleyg-lpptagwgmgidrltmfmadknnikevi 1243 (1262)
T tr|L1LCI1|L1LC 1206 ESFCTALEYG-LPPTAGWGMGIDRLTMFMADKNNIKEVI 1243 (1262)
Confidence            3456666555 9999999999999999999999999985


No 143
>tr|A0A0D7BRC0|A0A0D7BRC0_9HOMO Uncharacterized protein OS=Cylindrobasidium torrendii FP15055 ss-10 GN=CYLTODRAFT_19643 PE=4 SV=1
Probab=76.38  E-value=4.9  Score=33.14  Aligned_cols=24  Identities=25%  Similarity=0.584  Sum_probs=22.2  Template_Neff=1.000

Q d12asa_         282 GEMPQTIGGGIGQSRLTMLLLQLP  305 (327)
Q Consensus       282 ~~lP~TIgGGIgqSRl~M~lL~k~  305 (327)
                      .+||+.+|.|+..|--.||+++|.
T Consensus         8 sklplvvgagvalsaaamfmmqkd   31 (69)
T tr|A0A0D7BRC0|    8 SKLPLVVGAGVALSAAAMFMMQKD   31 (69)
Confidence            469999999999999999999985


No 144
>tr|A0A067BUS2|A0A067BUS2_SAPPC Uncharacterized protein (Fragment) OS=Saprolegnia parasitica (strain CBS 223.65) GN=SPRG_21703 PE=4 SV=1
Probab=76.32  E-value=5  Score=35.96  Aligned_cols=44  Identities=30%  Similarity=0.498  Sum_probs=35.4  Template_Neff=1.000

Q d12asa_          54 VKVKALPDAQFEVVHS--LAKWKRQTLGQHDFSAGEGLYTHMKALRPDE  100 (327)
Q Consensus        54 F~v~~~~~~~~EIVhS--LAKWKR~aL~ky~~~~geGiyTdMnAIRrDE  100 (327)
                      |.|...+=...|||..  |-||||....+.|.   .|+-|.|-.+|.-.
T Consensus        21 fvvstidfppieivsqeylikwkrsvmsrwgh---rglatsmvsvrsta   66 (113)
T tr|A0A067BUS2|   21 FVVSTIDFPPIEIVSQEYLIKWKRSVMSRWGH---RGLATSMVSVRSTA   66 (113)
Confidence            6666655556788754  88999999999987   89999999998554


No 145
>tr|A0A0F8YFG0|A0A0F8YFG0_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=LCGC14_2825740 PE=4 SV=1
Probab=76.07  E-value=5.1  Score=42.37  Aligned_cols=49  Identities=39%  Similarity=0.569  Sum_probs=34.3  Template_Neff=1.000

Q d12asa_         261 QLALTGDEDR--LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQV  310 (327)
Q Consensus       261 Ql~~~~~~~r--~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEV  310 (327)
                      ||...+..+.  ...+|-.++.-| +|.|-|-|||..||+|.|-....|..|
T Consensus        53 qlrkkgqlethqvdedfleameyg-mpptsgfgigidrlvmtltnspsiqdv  103 (403)
T tr|A0A0F8YFG0|   53 QLRKKGQLETHQVDEDFLEAMEYG-MPPTSGFGIGIDRLVMTLTNSPSIQDV  103 (403)
Confidence            4444444332  333454555445 899999999999999999988887766


No 146
>tr|A0A0G2YPK2|A0A0G2YPK2_BACIU Aspartyl/asparaginyl-tRNA synthetase OS=Bacillus subtilis GN=ABA10_03500 PE=3 SV=1
Probab=75.74  E-value=4.2  Score=44.13  Aligned_cols=260  Identities=17%  Similarity=0.195  Sum_probs=151.1  Template_Neff=3.400

Q d12asa_          22 ERLGLIEVQAPILSRVGDG------------------TQDNLSGAEKAVQVKVKALPDAQ-FEVVHSLAKWKRQTLGQHD   82 (327)
Q Consensus        22 ~~LnL~rVsaPLfv~~~sG------------------lNDdLnG~ErpV~F~v~~~~~~~-~EIVhSLAKWKR~aL~ky~   82 (327)
                      ...+.+++.+|++-..-+-                  .|||+.    |-..+|.-. |.. ..+..|.-=.|..++.-++
T Consensus        34 d~~GF~EilP~iigpvTDPg~rga~~~~~~~~~~~~~~~~~~~----~~~~~vd~Y-G~~~Y~l~~S~IlyKQ~a~~~~~  108 (347)
T tr|A0A0G2YPK2|   34 DEHGFVELLPPIIGPVTDPGARGAKIVLGDRILPEIRENQDDD----PSQADVDYY-GHEYYKLMTSAILYKQAALLAFD  108 (347)
Confidence            3457888888887554321                  122222    212333333 455 8899998888998888774


Q d12asa_          83 FSAGEGLYTHMKALRP---DEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAV--SEEFG--LA  155 (327)
Q Consensus        83 ~~~geGiyTdMnAIRr---DE~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v--~~~y~--l~  155 (327)
                      -     ||+=--.+|=   ++ ....-|=+-+-|+|-|..-...+--++.-..+|..+.+.++......  .+.++  +.
T Consensus       109 K-----Iy~iaPNVRLEp~e~-a~TgRHL~EF~QiDvE~~~asredvm~laE~ll~~~v~~v~~~~~~~~~L~~lgr~~~  182 (347)
T tr|A0A0G2YPK2|  109 K-----IFCIAPNVRLEPEET-ASTGRHLVEFFQIDVEWAGASREDVMDLAEDLLIHVVKYVLSEHPEERRLEVLGRDLP  182 (347)
Confidence            3     6665555553   33 46678889999999998888777776666666666665555554443  22232  21


Q d12asa_         156 P--FLP-DQIHFVHSQELLSRYP-DLDAK-----GRERAIAKDLG-AVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG  225 (327)
Q Consensus       156 ~--~Lp-~~I~FI~sqeL~~~YP-~Ltpk-----eRE~~i~ke~g-AvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~  225 (327)
                      .  .|- ....=|+.+|-.+++- ..+|+     +-|..+.+++. .+||+.+  +..+-.--|...|  ..--      
T Consensus       183 ~~~~l~~~pf~rit~~eAve~~~~~~~~~~ei~w~~E~~LS~~~~~PFwI~dy--P~gsRgFYdrE~p--pGvl------  252 (347)
T tr|A0A0G2YPK2|  183 AKELLAKGPFPRITHAEAVERLGRSQSPDAEIDWEGEAILSAQADRPFWITDY--PKGSRGFYDREDP--PGVL------  252 (347)
Confidence            1  011 1111133333333333 11111     23344444444 3455444  2222233444444  1111      


Q d12asa_         226 HAGLNGDILVWNPVLEDAFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       226 ~~gLNGDilv~n~~l~~a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k  304 (327)
                         .|=|++ |-.-.   =|++|=|+|- +.+.+..|.+..   ...-..|+.+...|-+|+|.|-|||-=||.=|+..-
T Consensus       253 ---~dfDLi-~PeG~---GEl~SGgeRe~e~e~Iv~R~re~---p~ky~wyl~~akeg~~~~SAGFGIGvERLtRyi~Gl  322 (347)
T tr|A0A0G2YPK2|  253 ---RDFDLI-YPEGF---GELISGGEREYEYERIIERMREK---PAKYKWYLELAKEGGIPQSAGFGIGVERLTRYICGL  322 (347)
Confidence               333433 33322   3889999996 457788888877   223445667777778999999999999999999999


Q d12asa_         305 PHIGQVQA  312 (327)
Q Consensus       305 ~HIgEVq~  312 (327)
                      .||-++.+
T Consensus       323 ~~V~~~~~  330 (347)
T tr|A0A0G2YPK2|  323 DAVWDARP  330 (347)
Confidence            88766544


No 147
>tr|A0A0D2IX76|A0A0D2IX76_9CHLO Uncharacterized protein (Fragment) OS=Monoraphidium neglectum GN=MNEG_15402 PE=4 SV=1
Probab=75.24  E-value=5.6  Score=34.94  Aligned_cols=39  Identities=46%  Similarity=0.610  Sum_probs=30.6  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ++|-..+.-| +|.+=|=|||-.||.|+|-....|..|-+
T Consensus        56 ydfitaleyg-mppcgglgigvdrlvmlltdspsirdvia   94 (99)
T tr|A0A0D2IX76|   56 YDFITALEYG-MPPCGGLGIGVDRLVMLLTDSPSIRDVIA   94 (99)
Confidence            3444444444 88888999999999999999988888754


No 148
>tr|W0U5P1|W0U5P1_9FIRM Aspartyl/glutamyl-tRNA(Asn/Gln) amidotransferase,B subunit OS=Ruminococcus bicirculans GN=gatB PE=3 SV=1
Probab=74.74  E-value=5  Score=48.23  Aligned_cols=261  Identities=18%  Similarity=0.221  Sum_probs=146.3  Template_Neff=2.645

Q d12asa_          23 RLGLIEVQAPILSRVGDGTQDNLSGAEKAV-QVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDED  101 (327)
Q Consensus        23 ~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV-~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~  101 (327)
                      +=+.++|..|.++..+.         |..- -|.+..- +..+-.-||-.-+|.++..-|+-      |-..-..=|||+
T Consensus       675 ~e~F~EI~TP~i~A~~~---------EG~~n~F~v~~~-g~~~~L~QsPQ~fKQ~~V~~fek------~FqIA~~fR~E~  738 (957)
T tr|W0U5P1|W0U5  675 EEEFTEIETPILVAGTD---------EGARNEFIVPYR-GKFYTLPQAPQQFKQMMVAGFEK------YFQIARCFRDED  738 (957)
Confidence            44788999999998762         3333 4777776 78899999999999996654432      222223345664


Q d12asa_         102 RLSPLHSVYVDQWDWERVMGDG-----ERQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQIHFVHSQELLSRYP  175 (327)
Q Consensus       102 ~ld~~HSiyVDQWDWEkvI~~~-----dRnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~FI~sqeL~~~YP  175 (327)
                      .-..-|=--+-|+|+|.---.+     +-+...|+.+|.+||..-+.+-+.+...-| +.  .=+-|++...++.+..  
T Consensus       739 ~~~~Rh~~Ef~~LD~EM~Y~dSM~dim~m~T~m~~~v~~~i~~~Y~~~m~~~~~D~P~lr--~v~cIkv~~~~q~~r~--  814 (957)
T tr|W0U5P1|W0U5  739 SRGDRHQPEFTQLDFEMAYIDSMQDIIDMNTAMFNEVVEKIYGTYKEAMDKYGCDRPDLR--IVKCIKVSAAEQIKRG--  814 (957)
Confidence            5556677778899998754422     345567899999999955554444433333 32  1233555544442211  


Q d12asa_         176 DLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTP---SE-LGHAGLN---GDILVWNP-VLEDAFELS  247 (327)
Q Consensus       176 ~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~---~~-~~~~gLN---GDilv~n~-~l~~a~ElS  247 (327)
                        +--+-+.   ++...--.+.++.--+       .--|.|-|.-.   .+ ..+..+|   |.++.|-. .+=.-+||.
T Consensus       815 --~KGq~~k---~eL~P~~~V~~~~fak-------ek~d~~~~~~TH~P~S~P~~Y~m~~~~~~~iA~qyDLi~~G~EI~  882 (957)
T tr|W0U5P1|W0U5  815 --SKGQIEK---KELRPAWVVDFPMFAK-------EKTDEGRWTFTHNPFSMPPFYDMNKHEGPIIAQQYDLILNGYEIG  882 (957)
Confidence              1111111   1111111111111100       00111223200   00 0111122   44444432 234578999


Q d12asa_         248 SMGIRVD-ADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAA  318 (327)
Q Consensus       248 SmGirVd-~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~  318 (327)
                      |=|.|.- .+.+..-.+..|-+.---..|+++-.=| +|.-=|=|||--||.|.||.|..|.||-  .+|..
T Consensus       883 ~G~~R~H~~e~~~a~~k~~G~~~ed~~~~~~A~~YG-~PPHGGi~~GldRLmM~L~kkssIRe~~--~FPk~  951 (957)
T tr|W0U5P1|W0U5  883 GGSQRAHKYEIQEATYKNMGYDKEDMKTMYKAFKYG-APPHGGIAWGLDRLMMILEKKSSIREVM--AFPKT  951 (957)
Confidence            9999984 4445555555454333233555665556 8888888999999999999999999874  44543


No 149
>tr|W4KAS8|W4KAS8_9HOMO Uncharacterized protein OS=Heterobasidion irregulare TC 32-1 GN=HETIRDRAFT_163241 PE=4 SV=1
Probab=74.47  E-value=6  Score=31.58  Aligned_cols=41  Identities=27%  Similarity=0.357  Sum_probs=32.5  Template_Neff=1.000

Q d12asa_          75 RQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDW  116 (327)
Q Consensus        75 R~aL~ky~~~~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDW  116 (327)
                      |+.|-+.+|+.-+|++|.|..|-+.. +-..----|..|+-|
T Consensus        17 ryilinnnfplfeglftsmypiykcq-eakknpvwyysqikw   57 (57)
T tr|W4KAS8|W4KA   17 RYILINNNFPLFEGLFTSMYPIYKCQ-EAKKNPVWYYSQIKW   57 (57)
Confidence            77888899999999999999999887 444444556667666


No 150
>tr|A0A0A2MJI3|A0A0A2MJI3_9FLAO Uncharacterized protein OS=Flavobacterium subsaxonicum WB 4.1-42 = DSM 21790 GN=Q766_16485 PE=4 SV=1
Probab=73.59  E-value=6.6  Score=31.87  Aligned_cols=31  Identities=16%  Similarity=0.315  Sum_probs=27.8  Template_Neff=1.000

Q d12asa_           9 ISFVKSHFSRQLEERLGLIEVQAPILSRVGD   39 (327)
Q Consensus         9 I~~iK~~F~~~L~~~LnL~rVsaPLfv~~~s   39 (327)
                      -+.||+.|++++.-+|.|++|+.-+|+.+..
T Consensus        15 akliknifqkqfvieldltkiskrmflnpsk   45 (62)
T tr|A0A0A2MJI3|   15 AKLIKNIFQKQFVIELDLTKISKRMFLNPSK   45 (62)
Confidence            3678999999999999999999999988763


No 151
>tr|A0A037UZF0|A0A037UZF0_9RHIZ Uncharacterized protein (Fragment) OS=Rhodomicrobium udaipurense JA643 GN=T281_07980 PE=4 SV=1
Probab=73.41  E-value=6.7  Score=37.39  Aligned_cols=37  Identities=38%  Similarity=0.608  Sum_probs=32.5  Template_Neff=1.000

Q d12asa_          61 DAQFEVVHSLAK-----WKRQTLGQHDFSAGEGLYTHMKALR   97 (327)
Q Consensus        61 ~~~~EIVhSLAK-----WKR~aL~ky~~~~geGiyTdMnAIR   97 (327)
                      .-.+..+|||+.     |.|..|.++|-..|-|.+|-|-+--
T Consensus        46 avqakmlhslariygvewnrrtlsefgaclgagvvtrmaaaf   87 (169)
T tr|A0A037UZF0|   46 AVQAKMLHSLARIYGVEWNRRTLSEFGACLGAGVVTRMAAAF   87 (169)
Confidence            346778999995     9999999999999999999998754


No 152
>tr|H9IUD4|H9IUD4_BOMMO Uncharacterized protein OS=Bombyx mori PE=4 SV=1
Probab=72.74  E-value=7.1  Score=45.37  Aligned_cols=109  Identities=29%  Similarity=0.340  Sum_probs=74.0  Template_Neff=1.000

Q d12asa_          95 ALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFS----------TLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQIH  163 (327)
Q Consensus        95 AIRrDE~~ld~~HSiyVDQWDWEkvI~~~dRnl~----------~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~  163 (327)
                      |+..|| +..|+||+.||...+..-+.+..|-+.          -|+.||.|.-+++..--+.= +.|. -...|.+-  
T Consensus       909 avntde-dwanlhsvvvdrmsydaeveknkrlmktieelrykkqdlkntvtkmqkamekytkkd-kefeakrkeledc--  984 (1109)
T tr|H9IUD4|H9IU  909 AVNTDE-DWANLHSVVVDRMSYDAEVEKNKRLMKTIEELRYKKQDLKNTVTKMQKAMEKYTKKD-KEFEAKRKELEDC--  984 (1109)
Confidence            567788 799999999999999888877777543          35677777766653221111 1122 11111110  


Q d12asa_         164 FVHSQELLSRYPDLDA---------KGRERAIAKDLGAVFLVGIGGKLSDGHR  207 (327)
Q Consensus       164 FI~sqeL~~~YP~Ltp---------keRE~~i~ke~gAvFi~gIG~~L~~G~~  207 (327)
                      --..+||..+|..|+.         ++||...-+-..|-+-..|..+|++.+.
T Consensus       985 kaeleelkqrykeldeecetcaeylkqreeqckrlkeakialeivdklsnqkv 1037 (1109)
T tr|H9IUD4|H9IU  985 KAELEELKQRYKELDEECETCAEYLKQREEQCKRLKEAKIALEIVDKLSNQKV 1037 (1109)
Confidence            1235789999999975         6788887777888999999999988764


No 153
>tr|D0MYK2|D0MYK2_PHYIT Aspartyl-tRNA synthetase, putative OS=Phytophthora infestans (strain T30-4) GN=PITG_03806 PE=3 SV=1
Probab=72.43  E-value=7.3  Score=39.73  Aligned_cols=55  Identities=29%  Similarity=0.440  Sum_probs=40.5  Template_Neff=1.000

Q d12asa_         254 DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       254 d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      |+.-|.+.....+...-.-..|--.+.-+.||-. |||||.-|+.|+.|.--.|..
T Consensus       211 dpkllmermdelgvpqvsmcnyiyslrldslphg-gggiglervvmlylglgnisg  265 (284)
T tr|D0MYK2|D0MY  211 DPKLLMERMDELGVPQVSMCNYIYSLRLDSLPHG-GGGIGLERVVMLYLGLGNISK  265 (284)
Confidence            5666666666666655555566666777778864 899999999999998876653


No 154
>tr|V9W5Z7|V9W5Z7_9BACL Putative membrane protein OS=Paenibacillus larvae subsp. larvae DSM 25430 GN=ERIC2_c15220 PE=4 SV=1
Probab=72.40  E-value=7.4  Score=34.27  Aligned_cols=19  Identities=53%  Similarity=0.833  Sum_probs=15.4  Template_Neff=1.000

Q d12asa_         286 QTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       286 ~TIgGGIgqSRl~M~lL~k  304 (327)
                      +.+|=|.|||||+.+||=-
T Consensus        48 levglglgqsrlailllfv   66 (99)
T tr|V9W5Z7|V9W5   48 LEVGLGLGQSRLAILLLFV   66 (99)
Confidence            3578899999999988743


No 155
>tr|W4XMC8|W4XMC8_STRPU Uncharacterized protein OS=Strongylocentrotus purpuratus PE=4 SV=1
Probab=71.43  E-value=8.1  Score=38.01  Aligned_cols=35  Identities=37%  Similarity=0.766  Sum_probs=28.9  Template_Neff=1.000

Q d12asa_         190 LGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELG  225 (327)
Q Consensus       190 ~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~  225 (327)
                      ...||-+.-|... -|..|..|+-.||||+.|+-..
T Consensus        61 aapvftvktgrsv-pgnrhnkrstnyddweeptyhd   95 (210)
T tr|W4XMC8|W4XM   61 AAPVFTVKTGRSV-PGNRHNKRSTNYDDWEEPTYHD   95 (210)
Confidence            3468999999887 5889999999999999776543


No 156
>tr|G6CVH5|G6CVH5_DANPL Aspartyl-tRNA synthetase OS=Danaus plexippus GN=KGM_10156 PE=3 SV=1
Probab=71.17  E-value=8.2  Score=40.04  Aligned_cols=65  Identities=31%  Similarity=0.382  Sum_probs=41.6  Template_Neff=1.000

Q d12asa_         244 FELSSMGIRVDADTLKHQLAL-TGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       244 ~ElSSmGirVd~~~L~~Ql~~-~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      -||-|-.-|+-+-+|...-.. -+.....-..|-.+..-| -|.-.|||||.-|++|+.|.-..|..
T Consensus       168 eeilsgaqrihepallteraqhhgidiskiaayieafrlg-cpphagggiglervvmlylgldnirk  233 (325)
T tr|G6CVH5|G6CV  168 EEILSGAQRIHEPALLTERAQHHGIDISKIAAYIEAFRLG-CPPHAGGGIGLERVVMLYLGLDNIRK  233 (325)
Confidence            355555566666665433222 222222334455555555 67889999999999999999888865


No 157
>tr|D7G3Y0|D7G3Y0_ECTSI Lysine--tRNA ligase OS=Ectocarpus siliculosus GN=LysRS PE=3 SV=1
Probab=70.96  E-value=8.4  Score=43.26  Aligned_cols=41  Identities=46%  Similarity=0.625  Sum_probs=35.2  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ...+|-.++..| +|.|.|=|||..||.|+|-.-..|..|.+
T Consensus       648 vdddfldalerg-mpptaglgigidrlimllagatsirdvia  688 (710)
T tr|D7G3Y0|D7G3  648 VDDDFLDALERG-MPPTAGLGIGIDRLIMLLAGATSIRDVIA  688 (710)
Confidence            455677777776 99999999999999999999999988865


No 158
>tr|B9SM60|B9SM60_RICCO Putative uncharacterized protein OS=Ricinus communis GN=RCOM_0878920 PE=4 SV=1
Probab=70.40  E-value=8.8  Score=31.67  Aligned_cols=34  Identities=35%  Similarity=0.449  Sum_probs=27.3  Template_Neff=1.000

Q d12asa_         274 EWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHI  307 (327)
Q Consensus       274 ~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HI  307 (327)
                      |.|+|+.....-.-=|||=|-|||.+.||---.|
T Consensus         2 plhrmlnkkqvvvygggggggsrlllilllwldi   35 (67)
T tr|B9SM60|B9SM    2 PLHRMLNKKQVVVYGGGGGGGSRLLLILLLWLDI   35 (67)
Confidence            6789988888877778889999999988754443


No 159
>tr|A0A0F7U2K7|A0A0F7U2K7_9EURO Uncharacterized protein OS=Penicillium brasilianum GN=PMG11_10446 PE=4 SV=1
Probab=70.26  E-value=8.9  Score=37.79  Aligned_cols=25  Identities=52%  Similarity=0.796  Sum_probs=22.3  Template_Neff=1.000

Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIG  308 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIg  308 (327)
                      -|.-.|.|+|.-||.|+||+-.|-.
T Consensus        46 apphagagvglerllmlllqlrhpe   70 (212)
T tr|A0A0F7U2K7|   46 APPHAGAGVGLERLLMLLLQLRHPE   70 (212)
Confidence            5778899999999999999998863


No 160
>tr|A0A0G5IVB7|A0A0G5IVB7_PSEAI Uncharacterized protein OS=Pseudomonas aeruginosa GN=PAERUG_E15_London_28_01_14_06260 PE=4 SV=1
Probab=70.18  E-value=8.6  Score=29.11  Aligned_cols=29  Identities=31%  Similarity=0.523  Sum_probs=26.1  Template_Neff=1.496

Q d12asa_         293 GQSRLTMLLLQLPHIGQVQAGVWPAAVRE  321 (327)
Q Consensus       293 gqSRl~M~lL~k~HIgEVq~svW~~~~~~  321 (327)
                      |.-|..|.|++-..--||--|+||+..++
T Consensus        11 g~~rftmilmedcdp~~vvksiwpegrve   39 (40)
T tr|A0A0G5IVB7|   11 GAKRFTMILMEDCDPVAVVKSIWPEGRVE   39 (40)
Confidence            67899999999999999999999987654


No 161
>tr|B9TDH8|B9TDH8_RICCO Putative uncharacterized protein OS=Ricinus communis GN=RCOM_1851920 PE=4 SV=1
Probab=69.42  E-value=9.6  Score=33.80  Aligned_cols=18  Identities=56%  Similarity=1.023  Sum_probs=13.9  Template_Neff=1.000

Q d12asa_         281 RGEMPQTIGGGIGQSRLTM  299 (327)
Q Consensus       281 ~~~lP~TIgGGIgqSRl~M  299 (327)
                      +|.-|.| |||||.||-|-
T Consensus        78 sgdqpet-gggighsracs   95 (102)
T tr|B9TDH8|B9TD   78 SGDQPET-GGGIGHSRACS   95 (102)
Confidence            4555654 99999999984


No 162
>tr|A0A060T8I3|A0A060T8I3_BLAAD ARAD1D03784p OS=Blastobotrys adeninivorans GN=GNLVRS02_ARAD1D03784g PE=4 SV=1
Probab=69.41  E-value=9.6  Score=38.61  Aligned_cols=50  Identities=34%  Similarity=0.463  Sum_probs=36.8  Template_Neff=1.000

Q d12asa_         246 LSSMGIRVDADTLKHQLALTGDED-------RLELEWHQALLRGEMPQTIGGGIGQSRL  297 (327)
Q Consensus       246 lSSmGirVd~~~L~~Ql~~~~~~~-------r~~~~~h~~ll~~~lP~TIgGGIgqSRl  297 (327)
                      -+..||-+|-+-|..|-+....+.       .....||+|+++++||+|-  -|-|-||
T Consensus        31 aadsgiaidnellrqqhkqqtqnaaavqaavqssntfhqmllnnelpqtq--qilqerl   87 (262)
T tr|A0A060T8I3|   31 AADSGIAIDNELLRQQHKQQTQNAAAVQAAVQSSNTFHQMLLNNELPQTQ--QILQERL   87 (262)
Confidence            456799999998888766543322       3457899999999999984  4666665


No 163
>tr|F0VCT3|F0VCT3_NEOCL Lysyl-tRNA synthetase, related OS=Neospora caninum (strain Liverpool) GN=NCLIV_012460 PE=3 SV=1
Probab=69.02  E-value=9.9  Score=46.86  Aligned_cols=61  Identities=36%  Similarity=0.521  Sum_probs=44.9  Template_Neff=1.000

Q d12asa_         250 GIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       250 GirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.......-+.|.+......--+.+|-.+++-| ||.|.|=|||..||.|+|-...-|..|-
T Consensus      2478 gltpetqdtknqaeakadrtvpnqefiaalmlg-lpptaglgigldrlvmlltdaplirdvv 2538 (2547)
T tr|F0VCT3|F0VC 2478 GLTPETQDTKNQAEAKADRTVPNQEFIAALMLG-LPPTAGLGIGLDRLVMLLTDAPLIRDVV 2538 (2547)
Confidence            333344445556554444334677888888888 9999999999999999998888777764


No 164
>tr|R8NW90|R8NW90_BACCE Uncharacterized protein OS=Bacillus cereus VD136 GN=IIW_02773 PE=4 SV=1
Probab=68.54  E-value=10  Score=30.00  Aligned_cols=34  Identities=35%  Similarity=0.427  Sum_probs=26.8  Template_Neff=1.000

Q d12asa_         274 EWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHI  307 (327)
Q Consensus       274 ~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HI  307 (327)
                      ....+-+.|+....-|||||..|-..+||-.+.-
T Consensus         8 gvmtmklagkvavvtgggigigrstalllaeqga   41 (53)
T tr|R8NW90|R8NW    8 GVMTMKLAGKVAVVTGGGIGIGRSTALLLAEQGA   41 (53)
Confidence            3344567788888999999999999999876543


No 165
>tr|X1K747|X1K747_9ZZZZ Uncharacterized protein OS=marine sediment metagenome GN=S03H2_64332 PE=4 SV=1
Probab=67.60  E-value=8.9  Score=36.60  Aligned_cols=75  Identities=29%  Similarity=0.424  Sum_probs=49.5  Template_Neff=3.464

Q d12asa_         236 WNPVLEDAFELSSMGIRVDADTLKHQ-LAL---TGDEDRLELEWH-QALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQV  310 (327)
Q Consensus       236 ~n~~l~~a~ElSSmGirVd~~~L~~Q-l~~---~~~~~r~~~~~h-~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEV  310 (327)
                      |+-+++ -.||++=.||+-...|.+- .+.   ...+...+..|- .++.-| -|.-=|=++|-.||+|+|..-..|.+|
T Consensus        20 YDlVLN-G~EiGgGSIRIh~~elQ~kvF~~lG~~~eea~ekFGFLleAlkyG-aPPHGGiAfGlDRL~mLl~g~~sIRdV   97 (133)
T tr|X1K747|X1K7   20 YDLVLN-GTEIGGGSIRIHNRELQEKVFEILGITDEEAEEKFGFLLEALKYG-APPHGGIAFGLDRLVMLLTGAESIRDV   97 (133)
Confidence            444443 4699999999987776543 222   333323333332 223333 576777788999999999999999999


Q d12asa_         311 QA  312 (327)
Q Consensus       311 q~  312 (327)
                      -|
T Consensus        98 IA   99 (133)
T tr|X1K747|X1K7   98 IA   99 (133)
Confidence            76


No 166
>tr|A0A093XZD4|A0A093XZD4_9PEZI Lysine--tRNA ligase OS=Pseudogymnoascus pannorum VKM F-3808 GN=O988_07417 PE=3 SV=1
Probab=67.40  E-value=9.5  Score=44.64  Aligned_cols=55  Identities=29%  Similarity=0.406  Sum_probs=44.1  Template_Neff=2.913

Q d12asa_         256 DTLKHQLALTGDED--RLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       256 ~~L~~Ql~~~~~~~--r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      .-+.+|+.....+.  ...-.|-+++..| ||.|=|=|-|--||||+|=....|++|-
T Consensus       511 rKF~eQ~~~k~~~~~~~iDesYv~AL~~G-LPPTGGWGCGveRLVMl~~Ga~RIsDvL  567 (645)
T tr|A0A093XZD4|  511 RKFVQQARWKDEENEAVIDESYVEALEWG-LPPTGGWGCGIDRLVMLFSGAKRISDVL  567 (645)
Confidence            34566666555443  3556788899888 9999999999999999999999999874


No 167
>tr|A0A093V0G4|A0A093V0G4_TALMA Aspartate--tRNA ligase, cytoplasmic OS=Talaromyces marneffei PM1 GN=GQ26_0350630 PE=3 SV=1
Probab=67.29  E-value=11  Score=42.06  Aligned_cols=213  Identities=20%  Similarity=0.298  Sum_probs=100.4  Template_Neff=1.000

Q d12asa_          54 VKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKAL------RPDEDRLSPL--HSVYVDQWDWERVMGDGER  125 (327)
Q Consensus        54 F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAI------RrDE~~ld~~--HSiyVDQWDWEkvI~~~dR  125 (327)
                      |.|+-. +..+-+.||-.--|.+|+.           -||..+      -|-|  -+|+  |=--+..+|||+.....-.
T Consensus       292 fevkyf-nkkgylaqspqlmkqmaia-----------gdmesvfevgpvfrae--nsnthrhlteftgldfektfrhhyh  357 (659)
T tr|A0A093V0G4|  292 FEVKYF-NKKGYLAQSPQLMKQMAIA-----------GDMESVFEVGPVFRAE--NSNTHRHLTEFTGLDFEKTFRHHYH  357 (659)
Confidence            677665 5677777887777776654           233221      1222  3444  4445678899987765433


Q d12asa_         126 -QFSTLKSTVEAIWAGIKAT----EAAVSEEFG-L-APFLPDQ-----IHFVHSQELLSRYPDLDAKGRERAIAKDLGAV  193 (327)
Q Consensus       126 -nl~~Lk~tV~kIy~al~~t----e~~v~~~y~-l-~~~Lp~~-----I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAv  193 (327)
                       -++|-.+.+-=|...|++-    -..+.+.|| - .-+||++     +.++.--.|++.-- .+-.|.| ++-..+-..
T Consensus       358 evldfaeellvfiltelkerykdeiaviqksypkagdfklpkdgkalrlnymdgvallkeag-vdtseqe-afendftta  435 (659)
T tr|A0A093V0G4|  358 EVLDFAEELLVFILTELKERYKDEIAVIQKSYPKAGDFKLPKDGKALRLNYMDGVALLKEAG-VDTSEQE-AFENDFTTA  435 (659)
Confidence             3555555554454444432    223455566 3 4466654     33333333333221 1111111 111111111


Q d12asa_         194 FLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDA-DTLKHQLALTGDEDR--  270 (327)
Q Consensus       194 Fi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~-~~L~~Ql~~~~~~~r--  270 (327)
                      .--.+|..+...-.....+-.||-+.          .            --||-|-.-|++. ..|.+-.-..+...+  
T Consensus       436 mekklgqiirekprrstfshsydffm----------r------------geeimsgaqrindvkeleesmikkgvdpkte  493 (659)
T tr|A0A093V0G4|  436 MEKKLGQIIREKPRDPTFSHSYDFFM----------R------------GEEIMSGAQRINDVKELEESMIKKGVDPKTE  493 (659)
Confidence            12223333333223333444454443          1            1122222223221 122222222222211  


Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLP  305 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~  305 (327)
                      --.+|-.+...| -|.-.|||+|..|++||+|.-.
T Consensus       494 gfedylnafrqg-cpphaggglglnrivmfflglp  527 (659)
T tr|A0A093V0G4|  494 GFEDYLNAFRQG-CPPHAGGGLGLNRIVMFFLGLP  527 (659)
Confidence            123455556666 6888999999999999999754


No 168
>tr|K3YKL8|K3YKL8_SETIT Uncharacterized protein OS=Setaria italica GN=Si014787m.g PE=4 SV=1
Probab=66.16  E-value=12  Score=31.29  Aligned_cols=16  Identities=44%  Similarity=0.515  Sum_probs=14.1  Template_Neff=1.000

Q d12asa_         294 QSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       294 qSRl~M~lL~k~HIgE  309 (327)
                      |.|+||.||+.+.|.|
T Consensus         6 qarvcmcllqeqgide   21 (72)
T tr|K3YKL8|K3YK    6 QARVCMCLLQEQGIDE   21 (72)
Confidence            7899999999887776


No 169
>tr|F6I4L4|F6I4L4_VITVI Lysine--tRNA ligase OS=Vitis vinifera GN=VIT_14s0060g02350 PE=3 SV=1
Probab=65.93  E-value=13  Score=44.72  Aligned_cols=58  Identities=34%  Similarity=0.448  Sum_probs=44.0  Template_Neff=1.000

Q d12asa_         253 VDADTLKHQLALTGDEDR----LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       253 Vd~~~L~~Ql~~~~~~~r----~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      |.+..+.+||+.+...+.    +...|.-++.-| ||.|=|=|+|..||.|.|-.-+.|.||-
T Consensus      1492 vqrqrfaeqlkdrqsgddeamaldetfcmaleyg-lpptggwglgidrltmmltdsqnikevl 1553 (1578)
T tr|F6I4L4|F6I4 1492 VQRQRFAEQLKDRQSGDDEAMALDETFCMALEYG-LPPTGGWGLGIDRLTMMLTDSQNIKEVL 1553 (1578)
Confidence            556677788877654442    334455555444 9999999999999999999999999984


No 170
>tr|C5LAS7|C5LAS7_PERM5 Aspartyl-tRNA synthetase, putative OS=Perkinsus marinus (strain ATCC 50983 / TXsc) GN=Pmar_PMAR008285 PE=3 SV=1
Probab=64.79  E-value=14  Score=37.52  Aligned_cols=80  Identities=33%  Similarity=0.474  Sum_probs=58.9  Template_Neff=1.000

Q d12asa_         229 LNGDILVWNPVLEDAFELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHI  307 (327)
Q Consensus       229 LNGDilv~n~~l~~a~ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HI  307 (327)
                      -.-||++-.      -|+.|-.-|+ |++.|.+.....+.+-..-.+|-.+..-|..|- .|.|+|.-|++|++|.-..|
T Consensus       170 nsydifirg------eevtsgaqrihdadmlaekaqrlgvelstiqpyvdafkygsyph-agagvglervvmlflglpni  242 (257)
T tr|C5LAS7|C5LA  170 NSYDIFIRG------EEVTSGAQRIHDADMLAEKAQRLGVELSTIQPYVDAFKYGSYPH-AGAGVGLERVVMLFLGLPNI  242 (257)
Confidence            444665543      3455555555 566666666677777778889999999999986 59999999999999999888


Q d12asa_         308 GQVQAGVWPA  317 (327)
Q Consensus       308 gEVq~svW~~  317 (327)
                      .=|  |.+|.
T Consensus       243 rlv--smfpr  250 (257)
T tr|C5LAS7|C5LA  243 RLV--SMFPR  250 (257)
Confidence            765  44443


No 171
>tr|A3TSB1|A3TSB1_OCEBH Putative agmatinase protein (Fragment) OS=Oceanicola batsensis (strain ATCC BAA-863 / DSM 15984 / HTCC2597) GN=OB2597_04630 PE=4 SV=1
Probab=63.47  E-value=15  Score=36.38  Aligned_cols=71  Identities=25%  Similarity=0.416  Sum_probs=55.3  Template_Neff=1.000

Q d12asa_         232 DILVWNPVLED----AFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPH  306 (327)
Q Consensus       232 Dilv~n~~l~~----a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~H  306 (327)
                      .|--||.++.-    ++-++..|    .--+......+.|.+..+..|.+-+-.|-||+++||--..||-...-|-+.|
T Consensus       102 rigpwnrvlgvaplqalrfadlg----dvpmrsrydltschedieafydrvvgagllplsvggdhsisrpviaalarrh  176 (213)
T tr|A3TSB1|A3TS  102 RIGPWNRVLGVAPLQALRFADLG----DVPMRSRYDLTSCHEDIEAFYDRVVGAGLLPLSVGGDHSISRPVIAALARRH  176 (213)
Confidence            34468877763    44555444    2345566778888888999999999999999999999999998888887776


No 172
>tr|A0A068VDE3|A0A068VDE3_COFCA Lysine--tRNA ligase OS=Coffea canephora GN=GSCOC_T00012421001 PE=3 SV=1
Probab=63.37  E-value=14  Score=43.05  Aligned_cols=41  Identities=46%  Similarity=0.584  Sum_probs=33.7  Template_Neff=2.100

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +..+|-.+|.=| +|.|-|=|||..||.|+|-.-+-|..|-|
T Consensus       616 lDeDFltALEYG-MPPt~GmGiGIDRLVMLLTnSaSIRDVIa  656 (664)
T tr|A0A068VDE3|  616 LDEDFLTALEYG-MPPTGGMGIGIDRLVMLLTDSASIRDVIA  656 (664)
Confidence            344555555555 99999999999999999999999998865


No 173
>tr|A0A0C3J1A7|A0A0C3J1A7_9RHOO Uncharacterized protein OS=Thauera sp. SWB20 GN=PO78_3329 PE=4 SV=1
Probab=63.15  E-value=16  Score=29.61  Aligned_cols=32  Identities=31%  Similarity=0.385  Sum_probs=28.1  Template_Neff=1.000

Q d12asa_          21 EERLGLIEVQAPILSRVGDGTQDNLSGAEKAV   52 (327)
Q Consensus        21 ~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV   52 (327)
                      +.+--|.||+-||+.-..+||---|||+|++-
T Consensus        23 sgkttlrrvsiplisghqsglyrvlngvekaq   54 (58)
T tr|A0A0C3J1A7|   23 SGKTTLRRVSIPLISGHQSGLYRVLNGVEKAQ   54 (58)
Confidence            44567899999999999999999999999873


No 174
>tr|A0A0F9MLR7|A0A0F9MLR7_9ZZZZ Uncharacterized protein OS=marine sediment metagenome GN=LCGC14_1139150 PE=4 SV=1
Probab=62.34  E-value=17  Score=27.71  Aligned_cols=28  Identities=21%  Similarity=0.345  Sum_probs=23.6  Template_Neff=1.000

Q d12asa_         120 MGDGERQFSTLKSTVEAIWAGIKATEAAVS  149 (327)
Q Consensus       120 I~~~dRnl~~Lk~tV~kIy~al~~te~~v~  149 (327)
                      |.+.||.  ||++|-.||.+.|+++++...
T Consensus         3 ikkkdre--flketyqkilkelkeadknyt   30 (42)
T tr|A0A0F9MLR7|    3 IKKKDRE--FLKETYQKILKELKEADKNYT   30 (42)
Confidence            5677887  999999999999999887543


No 175
>tr|U2R9P7|U2R9P7_9FIRM Uncharacterized protein OS=Oscillibacter sp. KLE 1745 GN=HMPREF1546_01368 PE=4 SV=1
Probab=61.92  E-value=17  Score=30.79  Aligned_cols=25  Identities=32%  Similarity=0.598  Sum_probs=21.0  Template_Neff=1.253

Q d12asa_         277 QALLRGEMPQTIGGGIGQSRLTMLL  301 (327)
Q Consensus       277 ~~ll~~~lP~TIgGGIgqSRl~M~l  301 (327)
                      ++...|-|-+-.|=|||.||.|.|-
T Consensus         7 kmf~rgilllll~vgi~isrfc~ft   31 (71)
T tr|U2R9P7|U2R9    7 KMFIRGILLLLLGVGIGISRFCIFT   31 (71)
Confidence            4566788888899999999999874


No 176
>tr|A0A0F8B681|A0A0F8B681_CERFI Aspartate--tRNA ligase cytoplasmic OS=Ceratocystis fimbriata f. sp. platani GN=drs-1 PE=3 SV=1
Probab=61.18  E-value=17  Score=43.10  Aligned_cols=247  Identities=25%  Similarity=0.349  Sum_probs=127.8  Template_Neff=1.313

Q d12asa_          24 LGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRL  103 (327)
Q Consensus        24 LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~l  103 (327)
                      .+.++|..|-+-...+     =+|.|   -|.|.-- +..+=..||-.--|.+++.- +|+.-.-|=   -..| -|  -
T Consensus       443 ~gfveintpklqpaat-----esgae---vfkvnyf-grraflaqspql~kqm~iaa-dfkrvyeig---pvfr-ae--n  506 (1054)
T tr|A0A0F8B681|  443 QGFVEINTPKLQPAAT-----ESGAE---VFKVNYF-GRRAFLAQSPQLVKQMTIAA-DFKRVYEIG---PVFR-AE--N  506 (1054)
Confidence            4556666554433332     12322   3666554 56666777765555555532 332111111   1122 22  4


Q d12asa_         104 SPLHSVYVD--QWDWERVMGDGERQ-FSTLKSTVEAIWAGIKATEA--AVSEEFG-LAP-FL--PDQIHFVHSQELLSR-  173 (327)
Q Consensus       104 d~~HSiyVD--QWDWEkvI~~~dRn-l~~Lk~tV~kIy~al~~te~--~v~~~y~-l~~-~L--p~~I~FI~sqeL~~~-  173 (327)
                      +|+|--.-.  .+|.|.-.....-. |....++.+.|+.+....-+  .+.+++| -.- .|  -.-|.|.+--+++.. 
T Consensus       507 snthrhlteyigld~emelq~dhyell~ivd~~lkni~~avq~mpela~vr~rwps~dlv~l~ktpii~f~~gi~mlrad  586 (1054)
T tr|A0A0F8B681|  507 SNTHRHLTEYIGLDLEMELQHDHYELLAIVDEMLKNIFAAVQSMPELAEVRKRWPSTDLVWLEKTPIISFETGIAMLRAD  586 (1054)
Confidence            556554333  34555443332222 23344556667777665532  3455565 211 11  123555554444332 


Q d12asa_         174 -----YPDLDAKG--RERAIAK-DLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNG-DILVWNPVLEDAF  244 (327)
Q Consensus       174 -----YP~Ltpke--RE~~i~k-e~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNG-Dilv~n~~l~~a~  244 (327)
                           --||+-..  |.-+.+| +||.-+.|      -|.-|...|.. |-.--   . ...-.|| ||++      +--
T Consensus       587 g~~ievadlstrdeirlgelv~~~ygtdyyi------ldkfpssarpf-yaqr~---g-ds~ftngfdifl------rgq  649 (1054)
T tr|A0A0F8B681|  587 GADIEVADLSTRDEIRLGELVKAEYGTDYYI------LDKFPSSARPF-YAQRL---G-DSNFTNGFDIFL------RGQ  649 (1054)
Confidence                 33455433  3334443 45554332      13223333321 11100   0 0001233 3433      456


Q d12asa_         245 ELSSMGIRV-DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       245 ElSSmGirV-d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k  304 (327)
                      ||||-|-|+ |+..|.+-+..++-.+--..+|-.....| -|.-=|.|.|.-|+.|+||+-
T Consensus       650 eissggqrihdaa~lrrs~r~sgia~agmeey~~~fdag-apphggag~glervlmlll~l  709 (1054)
T tr|A0A0F8B681|  650 EISSGGQRIHDAAKLRRSMRESGIAEAGMEEYLSGFDAG-APPHGGAGLGLERVLMLLLEL  709 (1054)
Confidence            999999998 66778888888888777777888888888 466679999999999999987


No 177
>tr|F3GQQ4|F3GQQ4_PSESJ Nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase (Fragment) OS=Pseudomonas syringae pv. pisi str. 1704B GN=cobT PE=4 SV=1
Probab=61.04  E-value=18  Score=28.14  Aligned_cols=22  Identities=32%  Similarity=0.527  Sum_probs=18.5  Template_Neff=1.000

Q d12asa_          23 RLGLIEVQAPILSRVGDGTQDN   44 (327)
Q Consensus        23 ~LnL~rVsaPLfv~~~sGlNDd   44 (327)
                      ...|.+..+||+|-+++|||-.
T Consensus         9 acsllecaapllvgpgtglnae   30 (47)
T tr|F3GQQ4|F3GQ    9 ACSLLECAAPLLVGPGTGLNAE   30 (47)
Confidence            3467889999999999999853


No 178
>tr|A0A0L9TWR3|A0A0L9TWR3_PHAAN Uncharacterized protein OS=Phaseolus angularis GN=LR48_Vigan02g116300 PE=4 SV=1
Probab=61.01  E-value=18  Score=35.87  Aligned_cols=32  Identities=25%  Similarity=0.455  Sum_probs=29.2  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLL  302 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL  302 (327)
                      .+.+|-+.+.||.|-++-||||-.--||-||-
T Consensus       139 sqsdyvqrvvngrlvyasgggitkeslcsfls  170 (211)
T tr|A0A0L9TWR3|  139 SQSDYVQRVVNGRLVYASGGGITKESLCSFLS  170 (211)
Confidence            56789999999999999999999999999863


No 179
>tr|I3KV14|I3KV14_ORENI Uncharacterized protein OS=Oreochromis niloticus GN=ST3GAL1 (9 of 9) PE=4 SV=1
Probab=59.33  E-value=20  Score=35.07  Aligned_cols=43  Identities=23%  Similarity=0.461  Sum_probs=37.9  Template_Neff=1.470

Q d12asa_         100 EDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIK  142 (327)
Q Consensus       100 E~~ld~~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~  142 (327)
                      |..|+++.-+--|..||=+-|....||..|.+.||.|+|+.+-
T Consensus        23 epfls~t~~~s~d~~~ww~~~q~~r~nf~~~~~t~~k~f~i~~   65 (166)
T tr|I3KV14|I3KV   23 EPFLSPTNNISEDEFNWWRRIQYNRQNFNFFNKTVDKVFQIFP   65 (166)
Confidence            3357888888899999999999999999999999999998764


No 180
>tr|A0A0P6D0T2|A0A0P6D0T2_9CRUS Putative Lysine--tRNA ligase OS=Daphnia magna PE=4 SV=1
Probab=59.32  E-value=20  Score=39.79  Aligned_cols=40  Identities=30%  Similarity=0.497  Sum_probs=33.3  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....|...+.-| ||.|=|=|.|..|++|||-....|.||-
T Consensus       507 idetfctsleyg-lpptggwgmgidrmamfltnsnnikevl  546 (568)
T tr|A0A0P6D0T2|  507 IDETFCTSLEYG-LPPTGGWGMGIDRMAMFLTNSNNIKEVL  546 (568)
Confidence            344566666555 9999999999999999999999999985


No 181
>tr|C0NGY7|C0NGY7_AJECG Aspartyl-tRNA synthetase OS=Ajellomyces capsulatus (strain G186AR / H82 / ATCC MYA-2454 / RMSCC 2432) GN=HCBG_02609 PE=3 SV=1
Probab=59.07  E-value=20  Score=42.78  Aligned_cols=220  Identities=18%  Similarity=0.257  Sum_probs=108.7  Template_Neff=1.452

Q d12asa_          54 VKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKA-------LRPDEDRLSPLHSVYVDQWDWERVMGDGERQ  126 (327)
Q Consensus        54 F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnA-------IRrDE~~ld~~HSiyVDQWDWEkvI~~~dRn  126 (327)
                      |.|+-- ...+-..||-.-+|.+.+.           -||..       .|..+ .-..-|=--+-.+|+|+.....-..
T Consensus       734 fev~yf-kr~aylaqspqlykqmcia-----------gdme~vfei~pvfrae~-snthrhltef~gldfektf~shyhe  800 (1009)
T tr|C0NGY7|C0NG  734 FEVTYF-KRKAYLAQSPQLYKQMCIA-----------GDMESVFEIGPVFRAEE-SNTHRHLTEFVGLDFEKTFQSHYHE  800 (1009)
Confidence            666654 4566677777666666553           23432       24333 2333344556689999987665433


Q d12asa_         127 -FSTLKSTVEAIWAGIKA----TEAAVSEEFG-L-APFLPDQ-----IHFVHSQ-ELLSRYPDLDAKGRE-RAIA----K  188 (327)
Q Consensus       127 -l~~Lk~tV~kIy~al~~----te~~v~~~y~-l-~~~Lp~~-----I~FI~sq-eL~~~YP~LtpkeRE-~~i~----k  188 (327)
                       |+|-....-=|...||.    --..+.+.|| - .-+||++     +.+..-- -|.+.--|+|.+||- +.++    |
T Consensus       801 vlefae~llvfils~lk~r~k~~i~iiq~sypkagdf~lpkdgkalrl~ym~gvallkeagvd~seqerfend~tt~mek  880 (1009)
T tr|C0NGY7|C0NG  801 VLEFAEDLLVFILSELKIRFKKEIEVIQRSYPKAGDFRLPKDGKALRLKYMEGVALLKEAGVDVTEQERFENDLSTAMEK  880 (1009)
Confidence             45555555555555543    2234566787 3 4567764     3333332 344555677877773 2221    2


Q d12asa_         189 --------DLGAVFLVG------IGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVD  254 (327)
Q Consensus       189 --------e~gAvFi~g------IG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd  254 (327)
                              +|..-|.+-      +-.--..--|.|.|-...-|+-         .-            --||-|-.-|++
T Consensus       881 qlg~iirekydtdfyvldkfpmavrpfytk~cp~dp~fsnsydff---------mr------------geeimsgaqrin  939 (1009)
T tr|C0NGY7|C0NG  881 QLGRIIREKYDTDFYVLDKFPMAVRPFYTKPCPDDPTFSNSYDFF---------MR------------GEEIMSGAQRIN  939 (1009)
Confidence                    222222210      0000001122333322222222         11            123333333333


Q d12asa_         255 A-DTLKHQLALTGDEDR--LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIG  308 (327)
Q Consensus       255 ~-~~L~~Ql~~~~~~~r--~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIg  308 (327)
                      . ..|..-....+....  --.+|-.+...| -|.-.|||.|-.|++||+|.-..|.
T Consensus       940 ~~kele~sm~akgldp~~egfedyi~afrqg-cpphaggglglnrivmfflglpnvr  995 (1009)
T tr|C0NGY7|C0NG  940 EAKELEAAMSAKGLDPNAEGFEDYIGAFRQG-CPPHAGGGLGLNRIVMFFLGLPNIR  995 (1009)
Confidence            2 122222222222111  223455666666 6888999999999999999876654


No 182
>tr|W7TE52|W7TE52_9STRA Aspartyl-trna cytoplasmic OS=Nannochloropsis gaditana GN=Naga_100431g3 PE=3 SV=1
Probab=57.51  E-value=20  Score=42.00  Aligned_cols=256  Identities=20%  Similarity=0.306  Sum_probs=149.1  Template_Neff=2.738

Q d12asa_          25 GLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLS  104 (327)
Q Consensus        25 nL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMnAIRrDE~~ld  104 (327)
                      +...|..|=++...+     =.| -.-  |.+.-- +..+=..||-.-.|.+|+.- ++.   |.+ ..--.=|-|..-+
T Consensus       312 gF~eIhTPklI~~~s-----egg-a~v--F~l~YF-g~~a~LAQSpQlyKQMa~~~-dl~---rVF-EIGPVFRAE~S~t  377 (601)
T tr|W7TE52|W7TE  312 GFHEIHTPKLITATS-----EGG-AEV--FHLDYF-GQDAWLAQSPQLYKQMALSM-DMQ---RVF-EIGPVFRAEVSKS  377 (601)
Confidence            788888887776553     122 112  555444 77888999999999997642 221   111 1111112221344


Q d12asa_         105 PLHSVYVDQWDWERVMGDGER-QFSTLKSTVEAIWAGIKAT----EAAVSEEFG-LAPFLPDQ---IHFVHSQELLSRYP  175 (327)
Q Consensus       105 ~~HSiyVDQWDWEkvI~~~dR-nl~~Lk~tV~kIy~al~~t----e~~v~~~y~-l~~~Lp~~---I~FI~sqeL~~~YP  175 (327)
                      .-|=.-+.-+|.|.-|.+.-. .++.+......||..|.+-    -..|.+.|| -.+.+|.+   |+|--.-++++.--
T Consensus       378 ~RHl~EFTglD~em~i~~hY~Evle~~~~~f~~iF~~L~~r~a~~l~~i~~qyp~~~~~~~~~~~~~~~~Ea~~~Lre~g  457 (601)
T tr|W7TE52|W7TE  378 RRHMTEFTGLDIEMAIQDHYHEVLELIESMFKFIFRGLQERYAPLLEIVQHQYPSAPPFEHGKPPRITFLEAKRILREEG  457 (601)
Confidence            456667788899988876543 3566667777777776543    234567787 34445544   44544444444333


Q d12asa_         176 D---------LDAKG--RERAIAK-DLGAVFLVGIGGKLSDGHRHDVRA------PDYDDWSTPSELGHAGLNGDILVWN  237 (327)
Q Consensus       176 ~---------Ltpke--RE~~i~k-e~gAvFi~gIG~~L~~G~~Hd~Ra------pDYDDW~t~~~~~~~gLNGDilv~n  237 (327)
                      +         +|..+  +.-++.| +||+-|.+-      |..|-+.|.      ||..-.+         -.=|+++  
T Consensus       458 ~~~e~~~~~d~~~~~E~~LGr~vkek~~tD~f~i------D~YP~s~RpFytm~~Pddp~~S---------Ns~D~f~--  520 (601)
T tr|W7TE52|W7TE  458 GMLESDDGKDFTDQEEAALGRHVREKYGTDFFTI------DQYPASMRPFYTMPCPDDPGFS---------NSYDTFV--  520 (601)
Confidence            2         22222  1222333 356555442      555666553      4443333         2223433  


Q d12asa_         238 PVLEDAFELSSMGIRVDA-DTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWP  316 (327)
Q Consensus       238 ~~l~~a~ElSSmGirVd~-~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~  316 (327)
                          .-=||.|=.=|+.. +-|.+-+...+-.-..-.+|-.+...| .|.--|+|||--|+.|++|.-..|.++  |.+|
T Consensus       521 ----rG~EI~SGaQRihd~dll~e~~~a~g~~~e~l~~Y~~af~~G-~ppHgG~GiGleRvv~lylgL~nvRk~--slFP  593 (601)
T tr|W7TE52|W7TE  521 ----RGREICSGAQRIHDYDLLREAMRAGGLDPDKWQPYLAAFKSG-MPPHGGCGLGLERVVQLFLGLDNVRKV--SLFP  593 (601)
Confidence                23467776667644 444454444443324556788888888 677789999999999999999988875  5555


Q d12asa_         317 AA  318 (327)
Q Consensus       317 ~~  318 (327)
                      .+
T Consensus       594 RD  595 (601)
T tr|W7TE52|W7TE  594 RD  595 (601)
Confidence            54


No 183
>tr|R5GV28|R5GV28_9BACT Uncharacterized protein OS=Prevotella sp. CAG:755 GN=BN773_00558 PE=4 SV=1
Probab=54.74  E-value=27  Score=29.57  Aligned_cols=27  Identities=30%  Similarity=0.535  Sum_probs=24.3  Template_Neff=1.000

Q d12asa_         297 LTMLLLQLPHIGQVQAGVWPAAVRESV  323 (327)
Q Consensus       297 l~M~lL~k~HIgEVq~svW~~~~~~~~  323 (327)
                      .+|+|+.|..|.|.-.|+|...++++-
T Consensus        41 famllmhkkkieeagtsvwassvleyq   67 (73)
T tr|R5GV28|R5GV   41 FAMLLMHKKKIEEAGTSVWASSVLEYQ   67 (73)
Confidence            689999999999999999999888763


No 184
>tr|L1JTG7|L1JTG7_GUITH Uncharacterized protein OS=Guillardia theta CCMP2712 GN=GUITHDRAFT_134183 PE=4 SV=1
Probab=52.68  E-value=29  Score=38.00  Aligned_cols=43  Identities=30%  Similarity=0.523  Sum_probs=33.9  Template_Neff=1.629

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAA  318 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~  318 (327)
                      -.|-.-+..|-+|.   |||||||+.-+--+..-|.||-+.+|..+
T Consensus        90 ~~~~~~~~~gii~~---~~~gqs~~~w~~r~h~~i~~~~~~~w~t~  132 (383)
T tr|L1JTG7|L1JT   90 GNWPRCVEGGIIPF---CGAGQSRCAWYVRSHHKIREVFAGLWQTE  132 (383)
Confidence            34445566666664   78999999988888888999999999754


No 185
>tr|E4TU80|E4TU80_MARTH Uncharacterized protein OS=Marivirga tractuosa (strain ATCC 23168 / DSM 4126 / NBRC 15989 / NCIMB 1408 / VKM B-1430 / H-43) GN=Ftrac_1011 PE=4 SV=1
Probab=52.50  E-value=32  Score=32.66  Aligned_cols=44  Identities=30%  Similarity=0.472  Sum_probs=31.9  Template_Neff=1.000

Q d12asa_          33 ILSRVGDGTQDNLSGAEKAV-QVKVKALPDAQFEVVHSL----AKWKRQT   77 (327)
Q Consensus        33 Lfv~~~sGlNDdLnG~ErpV-~F~v~~~~~~~~EIVhSL----AKWKR~a   77 (327)
                      -|...++|||||+..+-+-| +|.+.-. +..+||+.|-    .||--..
T Consensus        63 aftpngdglndnfgavakgvedfkiiiy-nrfgevifsstdietkwdgki  111 (145)
T tr|E4TU80|E4TU   63 AFTPNGDGLNDNFGAVAKGVEDFKIIIY-NRFGEVIFSSTDIETKWDGKI  111 (145)
Confidence            36678999999998655544 5777765 7888888875    4776543


No 186
>tr|A0A0G1C013|A0A0G1C013_9BACT Aspartyl-tRNA synthetase OS=Microgenomates (Beckwithbacteria) bacterium GW2011_GWA2_43_10 GN=UV54_C0047G0008 PE=3 SV=1
Probab=52.49  E-value=30  Score=37.94  Aligned_cols=120  Identities=27%  Similarity=0.363  Sum_probs=79.9  Template_Neff=1.513

Q d12asa_         187 AKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSEL--GHAGLNGDILVWNPVLEDAFELSSMGIRV-DADTLKHQLA  263 (327)
Q Consensus       187 ~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~--~~~gLNGDilv~n~~l~~a~ElSSmGirV-d~~~L~~Ql~  263 (327)
                      +|+-.--|.+|-|..---...|=--+|..||-.--..+  -..||.-|++     | .-+|..+-.||+ |+..+.+-..
T Consensus        78 tkqskedfyygsgqakfa~~~~~~t~~~p~dvplldknplk~~~l~~dlv-----~-~g~e~~~g~~rih~~~i~~~~f~  151 (394)
T tr|A0A0G1C013|   78 TKQSKEDFYYGSGQAKFAPSHHMFTAPHPDDVPLLDKNPLKVRGLQHDLV-----L-NGYEVGGGSIRIHDPKIQEKVFE  151 (394)
Confidence            34555556677666544444444445555543211111  1234444443     2 357888888887 5677778888


Q d12asa_         264 LTGDEDRLELEWHQALLRGE--MPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       264 ~~~~~~r~~~~~h~~ll~~~--lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      ..+-.++....||+++.+=.  +|.--|-.-|..||.|.|..-..|.||-+
T Consensus       152 ~~gftd~~~q~~~h~l~af~ygvp~~~g~~~g~~r~~~~~f~~~n~re~~~  202 (394)
T tr|A0A0G1C013|  152 LIGFTDKQKQQFHHMLEAFTYGVPPHGGIAPGIDRLLMVLFNEPNLREVMA  202 (394)
Confidence            88888899999999998643  56666666699999999999999999864


No 187
>tr|X1KLZ9|X1KLZ9_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S03H2_66725 PE=4 SV=1
Probab=52.45  E-value=26  Score=32.60  Aligned_cols=72  Identities=21%  Similarity=0.360  Sum_probs=55.4  Template_Neff=3.311

Q d12asa_         110 YVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG--L-APFLPDQIHFVHSQELLSRYPDLDAKGRE  183 (327)
Q Consensus       110 yVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~--l-~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE  183 (327)
                      +.|.-+|-+|+...--+..-+.++|++|++..+.==..-..+|-  . .-.|+.-  -++.+|++..+|.|+++.+.
T Consensus         6 ~~~~~~~~~v~~Rp~~~~~si~~tV~~I~e~V~~~GD~Ai~~ytekFD~v~L~~~--~v~~eei~~a~k~l~~~lk~   80 (106)
T tr|X1KLZ9|X1KL    6 NPSEEEWLEVLKRPRQDFSSIESTVRPIFENVKYRGDKAIKEYTEKFDGVRLDNS--IVSVEEIEAADKELDTELKN   80 (106)
Confidence            45778899999999889999999999999988754433333342  2 3356665  68999999999999988764


No 188
>tr|H1XVZ9|H1XVZ9_9BACT TonB family protein OS=Caldithrix abyssi DSM 13497 GN=Calab_2161 PE=4 SV=1
Probab=52.15  E-value=32  Score=37.78  Aligned_cols=70  Identities=23%  Similarity=0.240  Sum_probs=53.1  Template_Neff=1.000

Q d12asa_           3 IAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWK   74 (327)
Q Consensus         3 ~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWK   74 (327)
                      .|||..|.+=.-++.--|.++|-+.+-.||++-.++.|+||.=|-.-|-+-=  ++..-..+||+.||+-+-
T Consensus       376 aetqhqidfdqifldvfldkkltmkkkvapvvdnpergmndhgnviirvivg--kdgrveqaevlrslnyyy  445 (478)
T tr|H1XVZ9|H1XV  376 AETQHQIDFDQIFLDVFLDKKLTMKKKVAPVVDNPERGMNDHGNVIIRVIVG--KDGRVEQAEVLRSLNYYY  445 (478)
Confidence            4788999887777778899999999999999999999999987754433211  111135789999998653


No 189
>tr|A0A0A2KRV9|A0A0A2KRV9_PENIT Lysine--tRNA ligase OS=Penicillium italicum GN=PITC_020010 PE=3 SV=1
Probab=52.05  E-value=32  Score=41.23  Aligned_cols=46  Identities=35%  Similarity=0.510  Sum_probs=36.8  Template_Neff=1.287

Q d12asa_         265 TGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       265 ~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +++.......|-+++.=| ||.|=|=|-|-.||||++-.-..|+.|-
T Consensus       464 ~~e~~~idesyl~a~ewg-lp~tggwgcgvdrl~mlft~~kri~dvl  509 (1059)
T tr|A0A0A2KRV9|  464 TNESTGVDESYLQAMEWG-LPPTGGWGCGIDRLVMLFTDSKRIGDVL  509 (1059)
Confidence            334334556677777766 9999999999999999999999999873


No 190
>tr|E9LR59|E9LR59_9EUKA Lysyl-tRNA synthetase (Fragment) OS=Entamoeba nuttalli PE=3 SV=1
Probab=51.89  E-value=33  Score=37.61  Aligned_cols=32  Identities=38%  Similarity=0.581  Sum_probs=28.8  Template_Neff=1.000

Q d12asa_         280 LRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       280 l~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      +....|.|=|=|.|..||+|+|-.-..|.||-
T Consensus       217 iehafpptggwglgidrlamlladvdnikevi  248 (463)
T tr|E9LR59|E9LR  217 IEHAFPPTGGWGLGIDRLAMLLADVDNIKEVI  248 (463)
Confidence            34568999999999999999999999999986


No 191
>tr|J2U6W9|J2U6W9_9BURK Uncharacterized protein OS=Polaromonas sp. CF318 GN=PMI15_01266 PE=4 SV=1
Probab=51.83  E-value=32  Score=30.19  Aligned_cols=39  Identities=28%  Similarity=0.619  Sum_probs=30.7  Template_Neff=1.402

Q d12asa_         218 WSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDAD  256 (327)
Q Consensus       218 W~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~  256 (327)
                      |.-++.+|..-.|-|...|....+..+.+-..|||||-+
T Consensus        41 ~ayq~gsgpaptnedf~~ws~~veqrvkmkq~g~~~~ge   79 (82)
T tr|J2U6W9|J2U6   41 WAYQTGSGPAPTNEDFEQWSHIVEQRVKMKQIGIRPDGE   79 (82)
Confidence            344445555557777779999999999999999999865


No 192
>tr|E1ZP53|E1ZP53_CHLVA Putative uncharacterized protein OS=Chlorella variabilis GN=CHLNCDRAFT_138835 PE=3 SV=1
Probab=51.26  E-value=34  Score=38.62  Aligned_cols=40  Identities=35%  Similarity=0.527  Sum_probs=34.0  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQ  311 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq  311 (327)
                      ....|..++.-| ||.|=|=|+|..|..|+|-.|..|.||-
T Consensus       412 vdenfctaleyg-lpptggwglgidrftmmltdknnikevl  451 (624)
T tr|E1ZP53|E1ZP  412 VDENFCTALEYG-LPPTGGWGLGIDRFTMMLTDKNNIKEVL  451 (624)
Confidence            344566666666 9999999999999999999999999985


No 193
>tr|K0REP7|K0REP7_THAOC Uncharacterized protein OS=Thalassiosira oceanica GN=THAOC_28583 PE=4 SV=1
Probab=51.23  E-value=34  Score=30.35  Aligned_cols=35  Identities=37%  Similarity=0.447  Sum_probs=30.1  Template_Neff=1.000

Q d12asa_         179 AKGRERAIAKDLGAVFLVGIGGKLS--DGHRHDVRAP  213 (327)
Q Consensus       179 pkeRE~~i~ke~gAvFi~gIG~~L~--~G~~Hd~Rap  213 (327)
                      .|-||.+..+|.|||-.=.=|..|-  ||+.|.+|+.
T Consensus        27 fklreealireqgavtmgaegdalfvadgevhqgrsd   63 (94)
T tr|K0REP7|K0RE   27 FKLREEALIREQGAVTMGAEGDALFVADGEVHQGRSD   63 (94)
Confidence            4679999999999998888888875  9999999873


No 194
>tr|A0A0F8Z584|A0A0F8Z584_9ZZZZ Uncharacterized protein OS=marine sediment metagenome GN=LCGC14_2816120 PE=4 SV=1
Probab=51.13  E-value=34  Score=28.55  Aligned_cols=31  Identities=16%  Similarity=0.245  Sum_probs=27.3  Template_Neff=1.000

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQA   31 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsa   31 (327)
                      |.+||+.|++.||+.+.+.=..-|+|.+|+-
T Consensus        22 sareteiavkaiknalsekerrplglkkitk   52 (66)
T tr|A0A0F8Z584|   22 SARETEIAVKAIKNALSEKERRPLGLKKITK   52 (66)
Confidence            4689999999999999998888899988764


No 195
>tr|X1GNU4|X1GNU4_9ZZZZ Uncharacterized protein OS=marine sediment metagenome GN=S03H2_14990 PE=4 SV=1
Probab=50.65  E-value=35  Score=26.63  Aligned_cols=32  Identities=22%  Similarity=0.444  Sum_probs=26.7  Template_Neff=1.000

Q d12asa_         133 TVEAIWAGIKATEAAVSEEFGLAPFLPDQIHF  164 (327)
Q Consensus       133 tV~kIy~al~~te~~v~~~y~l~~~Lp~~I~F  164 (327)
                      ||+|+|+..+..++.++.+|-....+|..|-+
T Consensus         2 tvekvyekykhldkalcdkylvchnfpgnily   33 (46)
T tr|X1GNU4|X1GN    2 TVEKVYEKYKHLDKALCDKYLVCHNFPGNILY   33 (46)
Confidence            78999999999999999999744468887754


No 196
>tr|A0A067TRT7|A0A067TRT7_9AGAR Uncharacterized protein OS=Galerina marginata CBS 339.88 GN=GALMADRAFT_132550 PE=4 SV=1
Probab=50.19  E-value=36  Score=37.28  Aligned_cols=44  Identities=27%  Similarity=0.425  Sum_probs=31.2  Template_Neff=1.000

Q d12asa_          97 RPDEDRLSP-LHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGI  141 (327)
Q Consensus        97 RrDE~~ld~-~HSiyVDQWDWEkvI~~~dRnl~~Lk~tV~kIy~al  141 (327)
                      .|||+.-|| +.|+..||-|||+..... .+...|.-.|+|=|..-
T Consensus         9 qrdedqednflssyqldqcdwekrlrra-ktvnilalsvrkeyrek   53 (459)
T tr|A0A067TRT7|    9 QRDEDQEDNFLSSYQLDQCDWEKRLRRA-KTVNILALSVRKEYREK   53 (459)
Confidence            457755665 567789999999987654 35566777788877643


No 197
>tr|G0AIS0|G0AIS0_COLFT Uncharacterized protein OS=Collimonas fungivorans (strain Ter331) GN=CFU_1021 PE=4 SV=1
Probab=49.94  E-value=37  Score=34.75  Aligned_cols=25  Identities=36%  Similarity=0.664  Sum_probs=21.2  Template_Neff=1.000

Q d12asa_          91 THMKALRPDEDRLSPLHSVYVDQWD  115 (327)
Q Consensus        91 TdMnAIRrDE~~ld~~HSiyVDQWD  115 (327)
                      .||.||..|.+.-..+||+||-|+-
T Consensus        73 vdmhaihtdkelraalhsiyvaqfk   97 (251)
T tr|G0AIS0|G0AI   73 VDMHAIHTDKELRAALHSIYVAQFK   97 (251)
Confidence            5999999998555689999999974


No 198
>tr|A0A078D0F7|A0A078D0F7_BRANA BnaA03g32860D protein OS=Brassica napus GN=BnaA03g32860D PE=3 SV=1
Probab=49.86  E-value=37  Score=36.18  Aligned_cols=30  Identities=43%  Similarity=0.596  Sum_probs=27.4  Template_Neff=1.000

Q d12asa_         283 EMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       283 ~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      -+|.+-|-|+|..||.|+|-..+.|..|-+
T Consensus       195 gmppasgmglgidrlvmlltnsasirdvia  224 (353)
T tr|A0A078D0F7|  195 GMPPASGMGLGIDRLVMLLTNSASIRDVIA  224 (353)
Confidence            389999999999999999999999998865


No 199
>tr|T1C6F5|T1C6F5_9ZZZZ Asparagine synthetase A OS=mine drainage metagenome GN=B1B_01657 PE=4 SV=1
Probab=49.63  E-value=35  Score=32.88  Aligned_cols=79  Identities=24%  Similarity=0.423  Sum_probs=51.7  Template_Neff=1.809

Q d12asa_         229 LNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIG  308 (327)
Q Consensus       229 LNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIg  308 (327)
                      +|-|+ +|-.-.+.++.=+..-.-|++  +.+..+..+.. -.+..|.-.+..| |-.+.|-|||--||.-|+-.-.||+
T Consensus        49 ~~~dl-i~p~g~~e~lsg~ere~~~~r--i~~ri~r~g~~-~~~~k~~l~~a~g-l~psag~g~gver~~r~~~~~~~i~  123 (138)
T tr|T1C6F5|T1C6   49 RDMDL-IYPEGFGEALSGGEREFEVQR--IRERIARKGQP-EEQFKWYLEVAHG-LIPSAGFGIGVERLLRYLCGLPRIE  123 (138)
Confidence            56665 466667766654444443432  33334444432 2344555566666 5567899999999999999999999


Q d12asa_         309 QVQA  312 (327)
Q Consensus       309 EVq~  312 (327)
                      .||-
T Consensus       124 ~~~~  127 (138)
T tr|T1C6F5|T1C6  124 DVQP  127 (138)
Confidence            9984


No 200
>tr|T0REB3|T0REB3_9DELT WYL domain protein OS=Bacteriovorax sp. BSW11_IV GN=M899_1042 PE=4 SV=1
Probab=48.91  E-value=39  Score=35.82  Aligned_cols=30  Identities=23%  Similarity=0.458  Sum_probs=25.6  Template_Neff=1.000

Q d12asa_          84 SAGEGLYTHMKALRPDEDRLSPLHSVYVDQ  113 (327)
Q Consensus        84 ~~geGiyTdMnAIRrDE~~ld~~HSiyVDQ  113 (327)
                      -.+.|||.|.|.||||...++.+|-+.+.+
T Consensus        35 lenqgiyvdtntirrditslstthglvcse   64 (336)
T tr|T0REB3|T0RE   35 LENQGIYVDTNTIRRDITSLSTTHGLVCSE   64 (336)
Confidence            367899999999999997799999877654


No 201
>tr|V6LXH9|V6LXH9_9EUKA Uncharacterized protein OS=Spironucleus salmonicida GN=SS50377_10566 PE=4 SV=1
Probab=48.67  E-value=40  Score=32.64  Aligned_cols=35  Identities=23%  Similarity=0.313  Sum_probs=29.5  Template_Neff=1.000

Q d12asa_           9 ISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD   43 (327)
Q Consensus         9 I~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlND   43 (327)
                      |-..--.++..|+.+|.|.+.|-||+.++++|-|-
T Consensus        47 iyqadcaletalcaelqlaelsfplvmrketgana   81 (164)
T tr|V6LXH9|V6LX   47 IYQADCALETALCAELQLAELSFPLVMRKETGANA   81 (164)
Confidence            33344567889999999999999999999999874


No 202
>tr|A0A0G2FYK3|A0A0G2FYK3_9PEZI Putative lysyl-trna synthetase OS=Diaporthe ampelina GN=UCDDA912_g00790 PE=3 SV=1
Probab=48.35  E-value=41  Score=37.21  Aligned_cols=105  Identities=27%  Similarity=0.452  Sum_probs=63.4  Template_Neff=1.000

Q d12asa_           4 AKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL--PDAQFEVVHSLAKW-KRQTLGQ   80 (327)
Q Consensus         4 eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG~ErpV~F~v~~~--~~~~~EIVhSLAKW-KR~aL~k   80 (327)
                      +|....+.+-..+.+-|...|.-.+|..|++....-|      .+-||  |.....  +....-.--.-.-| ||.... 
T Consensus        16 dtlrvrscvvrklrqwletnlecvevhtpilaanagg------avarp--fttnatefpekqlalriapelwlkrlvva-   86 (490)
T tr|A0A0G2FYK3|   16 DTLRVRSCVVRKLRQWLETNLECVEVHTPILAANAGG------AVARP--FTTNATEFPEKQLALRIAPELWLKRLVVA-   86 (490)
Confidence            4555566677778888888999999999987654422      34556  433322  22222222222345 333332 


Q d12asa_          81 HDFSAGEGLYTHMKALRPDEDRLSPLHS-------VYVDQWDWERVMGDG  123 (327)
Q Consensus        81 y~~~~geGiyTdMnAIRrDE~~ld~~HS-------iyVDQWDWEkvI~~~  123 (327)
                       ||   .-||.=--|.|..  -+|.+|.       .||..+..|+++.+.
T Consensus        87 -gf---driyeigpafrne--gldsthnpeftmcefyvshinleqimqet  130 (490)
T tr|A0A0G2FYK3|   87 -GF---DRIYEIGPAFRNE--GLDSTHNPEFTMCEFYVSHINLEQIMQET  130 (490)
Confidence             33   4455555666644  4887774       799999999988654


No 203
>tr|M7XA29|M7XA29_9BACT Uncharacterized protein OS=Mariniradius saccharolyticus AK6 GN=C943_02143 PE=4 SV=1
Probab=48.23  E-value=41  Score=29.36  Aligned_cols=32  Identities=34%  Similarity=0.425  Sum_probs=27.6  Template_Neff=1.000

Q d12asa_         274 EWHQALLRGEMPQTIGGGIGQSRLTMLLLQLP  305 (327)
Q Consensus       274 ~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~  305 (327)
                      +--+-+.+|..|++.--||..|..+|-||+-.
T Consensus        46 ppprglrsgsfpqsahkgifksktamqlleld   77 (84)
T tr|M7XA29|M7XA   46 PPPRGLRSGSFPQSAHKGIFKSKTAMQLLELD   77 (84)
Confidence            33467889999999999999999999999853


No 204
>tr|I7ATP4|I7ATP4_9CAUD Uncharacterized protein OS=Escherichia phage ECML-117 PE=4 SV=1
Probab=48.17  E-value=34  Score=33.92  Aligned_cols=55  Identities=24%  Similarity=0.383  Sum_probs=42.4  Template_Neff=3.100

Q d12asa_         190 LGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSEL---GHAGLNGDILVWNPVLEDAF  244 (327)
Q Consensus       190 ~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~---~~~gLNGDilv~n~~l~~a~  244 (327)
                      .|+.+.+.--+..++|-.-..+|-|-|....|.-+   ..-|+|||+++|+...-..+
T Consensus         7 fG~~~~i~As~TFP~G~~it~FADD~DPlD~p~v~i~~~~m~~NG~lv~ws~~~p~~v   64 (155)
T tr|I7ATP4|I7AT    7 FGLSIVIVASKTFPNGFTITEFADDADPLDSPAVDIADTGMGLNGDLVVWSKATPLEV   64 (155)
Confidence            46666777777888888888899888888766654   55689999999997654433


No 205
>tr|W1SLJ3|W1SLJ3_9BACI Uncharacterized protein OS=Bacillus vireti LMG 21834 GN=BAVI_04664 PE=4 SV=1
Probab=48.12  E-value=41  Score=34.37  Aligned_cols=44  Identities=36%  Similarity=0.474  Sum_probs=39.2  Template_Neff=1.000

Q d12asa_         273 LEWHQALLRGEMPQTIGGGIGQSR-LTMLLLQLPHIGQVQAGVWP  316 (327)
Q Consensus       273 ~~~h~~ll~~~lP~TIgGGIgqSR-l~M~lL~k~HIgEVq~svW~  316 (327)
                      .+|...+.-|..|..|+|=|=||- +.++||+.+..|||-|+++-
T Consensus       184 snyydllqvgavplhisgiiiqsdslslylledaymgevfagvfi  228 (246)
T tr|W1SLJ3|W1SL  184 SNYYDLLQVGAVPLHISGIIIQSDSLSLYLLEDAYMGEVFAGVFI  228 (246)
Confidence            567888999999999999999984 88999999999999998764


No 206
>tr|A0A0C3PMV5|A0A0C3PMV5_PHLGI Uncharacterized protein OS=Phlebiopsis gigantea 11061_1 CR5-6 GN=PHLGIDRAFT_393685 PE=4 SV=1
Probab=48.02  E-value=41  Score=34.00  Aligned_cols=31  Identities=35%  Similarity=0.582  Sum_probs=27.6  Template_Neff=1.000

Q d12asa_          70 LAKWKRQTLGQHDFSAGEGLYTHMKALRPDE  100 (327)
Q Consensus        70 LAKWKR~aL~ky~~~~geGiyTdMnAIRrDE  100 (327)
                      ---|+|.+|..-....++|++|+.-|+|-.-
T Consensus       130 tlvwrrralgtqrvgahegfitnapalrggv  160 (227)
T tr|A0A0C3PMV5|  130 TLVWRRRALGTQRVGAHEGFITNAPALRGGV  160 (227)
Confidence            3469999999999999999999999999653


No 207
>tr|A0A0D7X2X4|A0A0D7X2X4_9BACL Uncharacterized protein OS=Paenibacillus terrae GN=QD47_09715 PE=4 SV=1
Probab=47.66  E-value=42  Score=35.01  Aligned_cols=23  Identities=30%  Similarity=0.689  Sum_probs=18.4  Template_Neff=1.000

Q d12asa_         221 PSELGHAGLNGDILVWNPVLEDA  243 (327)
Q Consensus       221 ~~~~~~~gLNGDilv~n~~l~~a  243 (327)
                      -..-.|+.|+||++|||...++|
T Consensus        50 ygdfryqqlygdlivwnqmidra   72 (291)
T tr|A0A0D7X2X4|   50 YGDFRYQQLYGDLIVWNQMIDRA   72 (291)
Confidence            33445788999999999988876


No 208
>tr|A0A0R2AJ50|A0A0R2AJ50_LACPA Uncharacterized protein OS=Lactobacillus paracasei subsp. paracasei DSM 5622 GN=FC74_GL002357 PE=4 SV=1
Probab=47.43  E-value=41  Score=28.33  Aligned_cols=32  Identities=25%  Similarity=0.279  Sum_probs=27.7  Template_Neff=1.437

Q d12asa_         166 HSQELLSRYPDLDAKGRERAIAKDLGAVFLVG  197 (327)
Q Consensus       166 ~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~g  197 (327)
                      |...+....-.||-.+|...+||.||-.+|-|
T Consensus        25 ~~~~~kq~igel~ns~rl~~fakahg~~vi~g   56 (64)
T tr|A0A0R2AJ50|   25 TTANLKQEIGELSNSERLDAFAKAHGFKVING   56 (64)
Confidence            55678888899999999999999999887754


No 209
>tr|D0NTV5|D0NTV5_PHYIT Aspartyl-tRNA synthetase, putative OS=Phytophthora infestans (strain T30-4) GN=PITG_16530 PE=3 SV=1
Probab=46.39  E-value=46  Score=35.97  Aligned_cols=62  Identities=26%  Similarity=0.433  Sum_probs=41.4  Template_Neff=1.000

Q d12asa_         254 DADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAA  318 (327)
Q Consensus       254 d~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~svW~~~  318 (327)
                      |+.-|.+.....++....-..|---+.-+.||-. |||||.-|+.|+..---.|.  .+|.+|.+
T Consensus       321 dpkllmermdelgapqesirnyiyylrldslphg-gggiglervvmlyfclgnic--kssmfprd  382 (388)
T tr|D0NTV5|D0NT  321 DPKLLMERMDELGAPQESIRNYIYYLRLDSLPHG-GGGIGLERVVMLYFCLGNIC--KSSMFPRD  382 (388)
Confidence            6666777777777665544445445555678764 89999999999887665553  35555544


No 210
>tr|E9EJY1|E9EJY1_METRA Methyl transferase OS=Metarhizium robertsii (strain ARSEF 23 / ATCC MYA-3075) GN=MAA_01386 PE=4 SV=2
Probab=46.21  E-value=42  Score=32.81  Aligned_cols=54  Identities=24%  Similarity=0.346  Sum_probs=46.0  Template_Neff=2.085

Q d12asa_         217 DWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDR  270 (327)
Q Consensus       217 DW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r  270 (327)
                      ||+..-.++..-.|.-+.=|...+=.|+..+-...||.|+...+||+.+|-.+-
T Consensus        18 dw~pr~~~~~~p~~~a~~~w~~~~l~amd~~~r~~rv~p~k~~~~l~~agf~di   71 (145)
T tr|E9EJY1|E9EJ   18 DWTPRWQGDDQPENSALQRWSEKLLSAMDKSGRSMRVVPEKTRRQLKAAGFTDI   71 (145)
Confidence            677555566666788888999999999999999999999999999999987653


No 211
>tr|A0A0Q9GYJ6|A0A0Q9GYJ6_9BACI Uncharacterized protein OS=Bacillus sp. Root131 GN=ASE54_06390 PE=4 SV=1
Probab=45.63  E-value=47  Score=32.73  Aligned_cols=36  Identities=36%  Similarity=0.405  Sum_probs=29.2  Template_Neff=1.300

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPH  306 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~H  306 (327)
                      +......|-++|+...--|||||..|-..+||-++.
T Consensus         5 l~egvm~mklagkvavvtgggigigrstallla~qg   40 (169)
T tr|A0A0Q9GYJ6|    5 LIEGVMKMKLAGKVAIVTGGGIGIGRNTALLLAKQG   40 (169)
Confidence            334445667788899999999999999999998765


No 212
>tr|A0A0G0YPD7|A0A0G0YPD7_9BACT Uncharacterized protein OS=Parcubacteria bacterium GW2011_GWA2_42_14 GN=UV01_C0003G0026 PE=4 SV=1
Probab=45.46  E-value=48  Score=34.09  Aligned_cols=55  Identities=25%  Similarity=0.408  Sum_probs=50.5  Template_Neff=1.000

Q d12asa_         154 LAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRH  208 (327)
Q Consensus       154 l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~H  208 (327)
                      +..+-|-.|.||-.-||-+..|.-.-.|...++||+..++|+.++|..|.=|+..
T Consensus        66 lqekspltiefiidpeltksipqenieenikevakqassifltqvgrrleigkvt  120 (255)
T tr|A0A0G0YPD7|   66 LQEKSPLTIEFIIDPELTKSIPQENIEENIKEVAKQASSIFLTQVGRRLEIGKVT  120 (255)
Confidence            5567888999999999999999999999999999999999999999999988754


No 213
>tr|A0A0Q7F9L2|A0A0Q7F9L2_9CAUL Uncharacterized protein OS=Phenylobacterium sp. Root1290 GN=ASC79_11505 PE=4 SV=1
Probab=44.99  E-value=50  Score=28.74  Aligned_cols=34  Identities=35%  Similarity=0.469  Sum_probs=30.6  Template_Neff=1.000

Q d12asa_         235 VWNPVLEDAFELSSMGIRVDADTLKHQLALTGDE  268 (327)
Q Consensus       235 v~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~  268 (327)
                      .--|..+++|.++..|---|-.+|++||+..+|.
T Consensus        15 mtiptierafqiarsgqardiatlkrqleadgcr   48 (81)
T tr|A0A0Q7F9L2|   15 MTIPTIERAFQIARSGQARDIATLKRQLEADGCR   48 (81)
Confidence            3458899999999999999999999999999984


No 214
>tr|B2ACG1|B2ACG1_PODAN Podospora anserina S mat+ genomic DNA chromosome 3, supercontig 1 OS=Podospora anserina (strain S / ATCC MYA-4624 / DSM 980 / FGSC 10383) GN=PODANS_3_923 PE=4 SV=1
Probab=44.29  E-value=52  Score=33.35  Aligned_cols=19  Identities=37%  Similarity=0.814  Sum_probs=17.1  Template_Neff=1.000

Q d12asa_         203 SDGHRHDVRAPDYDDWSTP  221 (327)
Q Consensus       203 ~~G~~Hd~RapDYDDW~t~  221 (327)
                      +||+||..|..-||.|..+
T Consensus        97 pdgqphsnreftydeweke  115 (224)
T tr|B2ACG1|B2AC   97 PDGQPHSNREFTYDEWEKE  115 (224)
Confidence            5999999999999999933


No 215
>tr|A0A0G0L7W7|A0A0G0L7W7_9BACT Uncharacterized protein OS=candidate division TM6 bacterium GW2011_GWF2_38_10 GN=US69_C0014G0002 PE=4 SV=1
Probab=43.60  E-value=52  Score=32.81  Aligned_cols=37  Identities=32%  Similarity=0.442  Sum_probs=32.9  Template_Neff=1.389

Q d12asa_         229 LNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALT  265 (327)
Q Consensus       229 LNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~  265 (327)
                      .-||.++||..-.+--.-|+.||.||-.++.-||...
T Consensus        53 trgdff~ynakkkrf~~ps~~gi~v~~~~~vg~l~~d   89 (180)
T tr|A0A0G0L7W7|   53 TRGDFFVYNAKKKRFSSPSSVGIKVDMNAAVGQLRMD   89 (180)
Confidence            6799999999988888899999999999998887654


No 216
>tr|A6DL15|A6DL15_9BACT Uncharacterized protein OS=Lentisphaera araneosa HTCC2155 GN=LNTAR_20463 PE=4 SV=1
Probab=41.85  E-value=60  Score=30.88  Aligned_cols=30  Identities=27%  Similarity=0.420  Sum_probs=25.9  Template_Neff=1.000

Q d12asa_         169 ELLSRYPDLDAKGRERAIAKDLGAVFLVGI  198 (327)
Q Consensus       169 eL~~~YP~LtpkeRE~~i~ke~gAvFi~gI  198 (327)
                      .+.+.+|.=||+|-|+.|+-+.||||---+
T Consensus        44 aikekwpegtprelensiaiqwgavfsrml   73 (140)
T tr|A6DL15|A6DL   44 AIKEKWPEGTPRELENSIAIQWGAVFSRML   73 (140)
Confidence            356789999999999999999999996443


No 217
>tr|E4XWZ8|E4XWZ8_OIKDI Uncharacterized protein OS=Oikopleura dioica GN=GSOID_T00007176001 PE=4 SV=1
Probab=41.61  E-value=61  Score=37.07  Aligned_cols=36  Identities=31%  Similarity=0.559  Sum_probs=31.5  Template_Neff=1.000

Q d12asa_         263 ALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLT  298 (327)
Q Consensus       263 ~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqSRl~  298 (327)
                      ++.++++-.+..|..++.+|-||.+=.||-|.||.-
T Consensus        46 eieneedviknqydtalangllppagsggrgrsrtn   81 (652)
T tr|E4XWZ8|E4XW   46 EIENEEDVIKNQYDTALANGLLPPAGSGGRGRSRTN   81 (652)
Confidence            456677778899999999999999999999999963


No 218
>tr|A0A067PEU0|A0A067PEU0_9HOMO Uncharacterized protein OS=Jaapia argillacea MUCL 33604 GN=JAAARDRAFT_201251 PE=4 SV=1
Probab=41.10  E-value=63  Score=31.69  Aligned_cols=30  Identities=23%  Similarity=0.490  Sum_probs=26.3  Template_Neff=1.000

Q d12asa_          39 DGTQDNLSGAEKAVQVKVKALPDAQFEVVH   68 (327)
Q Consensus        39 sGlNDdLnG~ErpV~F~v~~~~~~~~EIVh   68 (327)
                      .-|.|||||-..-|.|.+.+......-|+|
T Consensus       125 qqlrddlnglqhlveftiadpqpatfpiih  154 (172)
T tr|A0A067PEU0|  125 QQLRDDLNGLQHLVEFTIADPQPATFPIIH  154 (172)
Confidence            358899999999999999987778888887


No 219
>tr|A0A0N0Y7P1|A0A0N0Y7P1_9ACTN rRNA cytosine-C5-methyltransferase (Fragment) OS=Streptomyces sp. NRRL F-6602 GN=ADL27_53365 PE=4 SV=1
Probab=40.74  E-value=64  Score=29.25  Aligned_cols=48  Identities=25%  Similarity=0.344  Sum_probs=34.4  Template_Neff=1.000

Q d12asa_          42 QDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGL   89 (327)
Q Consensus        42 NDdLnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGi   89 (327)
                      .|||.|.-.-|.=.....++...-||||...|--.||...==.+..||
T Consensus        10 addldgwlervappydedpedhlavvhshprwvvsalwdslggpragi   57 (102)
T tr|A0A0N0Y7P1|   10 ADDLDGWLERVAPPYDEDPEDHLAVVHSHPRWVVSALWDSLGGPRAGI   57 (102)
Confidence            478888666555445555678999999999999999876543444444


No 220
>tr|A0A0F9V184|A0A0F9V184_9ZZZZ Uncharacterized protein OS=marine sediment metagenome GN=LCGC14_0142350 PE=4 SV=1
Probab=40.51  E-value=62  Score=31.18  Aligned_cols=43  Identities=23%  Similarity=0.398  Sum_probs=34.1  Template_Neff=1.501

Q d12asa_          62 AQFEVVHSLAKWKR---QTLGQHDFSAGEGLYTHMKALRPDEDRLS  104 (327)
Q Consensus        62 ~~~EIVhSLAKWKR---~aL~ky~~~~geGiyTdMnAIRrDE~~ld  104 (327)
                      ...-.-|||.||--   .-+++.+...++|.+-|-+.+-.||++.|
T Consensus        20 k~~a~rhs~~kw~gl~~~n~~kh~v~l~~g~~~~~~dvtddeddcd   65 (137)
T tr|A0A0F9V184|   20 KGWAMRHSIDKWTGLKCRNRRKHKVNLDEGNVYDNNDVTDDEDDCD   65 (137)
Confidence            34557899999954   45778888899999999999988886544


No 221
>tr|J9E6C3|J9E6C3_WUCBA Uncharacterized protein (Fragment) OS=Wuchereria bancrofti GN=WUBG_18181 PE=4 SV=1
Probab=39.70  E-value=68  Score=29.89  Aligned_cols=48  Identities=23%  Similarity=0.402  Sum_probs=40.5  Template_Neff=1.000

Q d12asa_          45 LSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMK   94 (327)
Q Consensus        45 LnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~~geGiyTdMn   94 (327)
                      ..|+--||...|.+.+.....|+||-.  +|..|.+.|+..+.-|+-++.
T Consensus        50 iggvlvpvqltvqslnapetrilhssd--rrktldelgidanshiiinfq   97 (121)
T tr|J9E6C3|J9E6   50 IGGVLVPVQLTVQSLNAPETRILHSSD--RRKTLDELGIDANSHIIINFQ   97 (121)
Confidence            458999999999999878889999976  688999999998887776653


No 222
>tr|A0A0J6EHQ9|A0A0J6EHQ9_9BACI Uncharacterized protein OS=Bacillus glycinifermentans GN=AB447_05940 PE=4 SV=1
Probab=39.12  E-value=70  Score=27.75  Aligned_cols=28  Identities=29%  Similarity=0.526  Sum_probs=22.0  Template_Neff=1.000

Q d12asa_         196 VGIGGKLSDGHRHDVRAPDYDDWSTPSE  223 (327)
Q Consensus       196 ~gIG~~L~~G~~Hd~RapDYDDW~t~~~  223 (327)
                      ..-|-+|+..+.||.-+|||..=+..+.
T Consensus        34 fnkgtklkenkahdlaspdyrgstkmtp   61 (78)
T tr|A0A0J6EHQ9|   34 FNKGTKLKENKAHDLASPDYRGSTKMTP   61 (78)
Confidence            3458899999999999999987654333


No 223
>tr|A0A0P9J0K9|A0A0P9J0K9_9PSED Uncharacterized protein OS=Pseudomonas coronafaciens pv. atropurpurea GN=ALO66_00479 PE=4 SV=1
Probab=38.53  E-value=73  Score=32.73  Aligned_cols=28  Identities=39%  Similarity=0.541  Sum_probs=25.7  Template_Neff=1.000

Q d12asa_         172 SRYPDLDAKGRERAIAKDLGAVFLVGIG  199 (327)
Q Consensus       172 ~~YP~LtpkeRE~~i~ke~gAvFi~gIG  199 (327)
                      .++|+=+|.-|..+|+|.+|+.|-.||=
T Consensus       200 arfpdgppsvrldeitkdngstfelgif  227 (241)
T tr|A0A0P9J0K9|  200 ARFPDGPPSVRLDEITKDNGSTFELGIF  227 (241)
Confidence            5789999999999999999999988873


No 224
>tr|E2NUA3|E2NUA3_9FIRM Uncharacterized protein OS=Catenibacterium mitsuokai DSM 15897 GN=CATMIT_02178 PE=4 SV=1
Probab=38.33  E-value=74  Score=25.96  Aligned_cols=22  Identities=27%  Similarity=0.490  Sum_probs=18.8  Template_Neff=1.000

Q d12asa_          87 EGLYTHMKALRPDEDRLSPLHS  108 (327)
Q Consensus        87 eGiyTdMnAIRrDE~~ld~~HS  108 (327)
                      -|+...+||||+-.+.+|-+|-
T Consensus        21 igilasfnairkyqesldilhe   42 (55)
T tr|E2NUA3|E2NU   21 IGILASFNAIRKYQESLDILHE   42 (55)
Confidence            5899999999998877887774


No 225
>tr|L0L8H6|L0L8H6_9CAUD Uncharacterized protein OS=Bacillus phage phiAGATE PE=4 SV=1
Probab=37.94  E-value=69  Score=31.53  Aligned_cols=48  Identities=33%  Similarity=0.414  Sum_probs=38.5  Template_Neff=1.936

Q d12asa_         133 TVEAIWAGIKATEAAVSEEFG-LAPFLPDQIHFVHSQELLSRYPDLDAKG  181 (327)
Q Consensus       133 tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~FI~sqeL~~~YP~Ltpke  181 (327)
                      |-++|.+-|+++|+.+.+.-. ++ .|-.++.=..+-+|..+||.||...
T Consensus         3 tke~il~el~~~ek~~k~~~~~~k-~l~~el~~~~~~~l~srfphlt~~~   51 (148)
T tr|L0L8H6|L0L8    3 TKEEILAELAQAEKEVKETTNQLK-ALQQELSEASSSDLRSRFPHLTADT   51 (148)
Confidence            456888999999999977665 44 5667777788899999999999753


No 226
>tr|A0A0F0HS90|A0A0F0HS90_9PSEU Uncharacterized protein OS=Saccharothrix sp. ST-888 GN=UK12_17830 PE=4 SV=1
Probab=37.85  E-value=76  Score=28.69  Aligned_cols=60  Identities=25%  Similarity=0.508  Sum_probs=48.6  Template_Neff=1.000

Q d12asa_         226 HAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIG  289 (327)
Q Consensus       226 ~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIg  289 (327)
                      -.|..|.--+|...-    -|+...-|.|+++|.+-........|.-..||+++....||+.-|
T Consensus        29 dtgvegkaalwralr----aiadqdprldpdaldrlaqraddqartlnrwhravsaqvlpqdeg   88 (99)
T tr|A0A0F0HS90|   29 DTGVEGKAALWRALR----AIADQDPRLDPDALDRLAQRADDQARTLNRWHRAVSAQVLPQDEG   88 (99)
Confidence            346778888886532    356777899999999887778887899999999999999998765


No 227
>tr|S8DGJ7|S8DGJ7_9LAMI Uncharacterized protein (Fragment) OS=Genlisea aurea GN=M569_16328 PE=4 SV=1
Probab=37.75  E-value=77  Score=32.64  Aligned_cols=52  Identities=27%  Similarity=0.417  Sum_probs=37.0  Template_Neff=1.000

Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAG--IKATEAAVSEEFGLAPFLPDQIHFVHSQE  169 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~a--l~~te~~v~~~y~l~~~Lp~~I~FI~sqe  169 (327)
                      ||.++..+.-++.|-|...+.+.+.  |+.+--.+++.  +.-.||++|.|+.+|.
T Consensus        60 wellvcdetlslqytkyflnnvinsvtlknalisicee--lgvplpekirffrsqm  113 (243)
T tr|S8DGJ7|S8DG   60 WELLVCDETLSLQYTKYFLNNVINSVTLKNALISICEE--LGVPLPEKIRFFRSQM  113 (243)
Confidence            8888888888888888777776654  34444444443  3346999999998874


No 228
>tr|A0A0B0EA93|A0A0B0EA93_9BACT Putative geranylgeranyl reductase OS=Candidatus Scalindua brodae GN=SCABRO_04027 PE=4 SV=1
Probab=37.33  E-value=78  Score=31.60  Aligned_cols=53  Identities=21%  Similarity=0.273  Sum_probs=38.7  Template_Neff=1.000

Q d12asa_         271 LELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHI--GQVQAGVWPAAVRESV  323 (327)
Q Consensus       271 ~~~~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HI--gEVq~svW~~~~~~~~  323 (327)
                      .+..-|-.+.||.--..||||-.-|=.|+|||+-+.=  ..|.-.++++.-...|
T Consensus        18 pkdnhhlqlrngsrvavigggpagsffclfllrfakrvgkdvsidiyddkdfskc   72 (193)
T tr|A0A0B0EA93|   18 PKDNHHLQLRNGSRVAVIGGGPAGSFFCLFLLRFAKRVGKDVSIDIYDDKDFSKC   72 (193)
Confidence            4456788899999999999999999999999985432  2344445555444333


No 229
>tr|A0A0R0HVZ6|A0A0R0HVZ6_SOYBN Uncharacterized protein OS=Glycine max GN=GLYMA_10G089700 PE=4 SV=1
Probab=37.26  E-value=76  Score=33.54  Aligned_cols=54  Identities=28%  Similarity=0.305  Sum_probs=41.4  Template_Neff=1.334

Q d12asa_         116 WERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQE  169 (327)
Q Consensus       116 WEkvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~l~~~Lp~~I~FI~sqe  169 (327)
                      ||.|+-...-++.|-|..-+.+.+.+---|..++---.|.-.+|++|.|+.+|.
T Consensus       115 welvvcd~~lslqytkyfpnnvinsitlk~aiva~~d~lgvp~p~~irffrsqm  168 (272)
T tr|A0A0R0HVZ6|  115 WELVVCDKTLSLQYTKYFPNNVINSITLKDAIVAVSDQLGVPLPRNIRFFRSQM  168 (272)
Confidence            999999999999999988888887765555554422225557999999999874


No 230
>tr|G9MSI5|G9MSI5_HYPVG Uncharacterized protein OS=Hypocrea virens (strain Gv29-8 / FGSC 10586) GN=TRIVIDRAFT_201253 PE=4 SV=1
Probab=37.19  E-value=79  Score=35.01  Aligned_cols=60  Identities=25%  Similarity=0.408  Sum_probs=50.7  Template_Neff=1.000

Q d12asa_         119 VMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEF-G-L---APFLPDQIHFVH--SQELLSRYPDLDA  179 (327)
Q Consensus       119 vI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y-~-l---~~~Lp~~I~FI~--sqeL~~~YP~Ltp  179 (327)
                      +..+.-|.++.|..+++|+.+.|..+.+.+.+ | | +   +..|-.+|.|-|  +||....-|.++.
T Consensus        55 llekqprslerlnsiiekvhnglldasrivek-yrpelhsgkkslqtqiawahrdsqefqalgptmsq  121 (456)
T tr|G9MSI5|G9MS   55 LLEKQPRSLERLNSIIEKVHNGLLDASRIVEK-YRPELHSGKKSLQTQIAWAHRDSQEFQALGPTMSQ  121 (456)
Confidence            56788999999999999999999999887764 6 6 6   568889999876  7888888888764


No 231
>tr|A0A0F9Y4R0|A0A0F9Y4R0_TRIHA Uncharacterized protein OS=Trichoderma harzianum GN=THAR02_00902 PE=4 SV=1
Probab=37.01  E-value=80  Score=34.40  Aligned_cols=62  Identities=24%  Similarity=0.409  Sum_probs=51.8  Template_Neff=1.000

Q d12asa_         118 RVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEF-G-L---APFLPDQIHFVH--SQELLSRYPDLDAK  180 (327)
Q Consensus       118 kvI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y-~-l---~~~Lp~~I~FI~--sqeL~~~YP~Ltpk  180 (327)
                      ++..+.-.+++.|..+++|+.+.|.++.+.+.+ | | +   +..|-.+|.|-|  +||...+-|.++..
T Consensus        54 rllekqpqtlerlnsiiekvhnglleasrivek-yrpelhsgkkslqtqiawahrdsqefqamgptmsqh  122 (390)
T tr|A0A0F9Y4R0|   54 RLLEKQPQTLERLNSIIEKVHNGLLEASRIVEK-YRPELHSGKKSLQTQIAWAHRDSQEFQAMGPTMSQH  122 (390)
Confidence            456777889999999999999999999887764 6 6 6   568889999875  78999999988743


No 232
>tr|U1X922|U1X922_9BURK Methyltransferase domain protein (Fragment) OS=Burkholderia cenocepacia BC7 GN=BURCENBC7_AP4763 PE=4 SV=1
Probab=36.91  E-value=81  Score=29.08  Aligned_cols=43  Identities=23%  Similarity=0.458  Sum_probs=38.6  Template_Neff=1.000

Q d12asa_           5 KQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG   47 (327)
Q Consensus         5 Tq~aI~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlNDdLnG   47 (327)
                      |-.|+.++-.+..+++.++|..++||+.-++..+-|-.|||-.
T Consensus        31 tfdavaflpreiaqrmnerleyikvspaavldagcgpgddlpa   73 (111)
T tr|U1X922|U1X9   31 TFDAVAFLPREIAQRMNERLEYIKVSPAAVLDAGCGPGDDLPA   73 (111)
Confidence            4567788888999999999999999999999999999999964


No 233
>tr|G8M0B5|G8M0B5_CLOCD Uncharacterized protein OS=Clostridium clariflavum (strain DSM 19732 / NBRC 101661 / EBR45) GN=Clocl_0176 PE=4 SV=1
Probab=35.54  E-value=83  Score=28.74  Aligned_cols=55  Identities=13%  Similarity=0.252  Sum_probs=44.0  Template_Neff=1.545

Q d12asa_         125 RQFSTLKSTVEAIWAGIKATEAAVSEEFG-LAPFLPDQIHFVHSQELLSRYPDLDA  179 (327)
Q Consensus       125 Rnl~~Lk~tV~kIy~al~~te~~v~~~y~-l~~~Lp~~I~FI~sqeL~~~YP~Ltp  179 (327)
                      .+++|--.+...+|..=+.-++.+..... |-.+|.++|.=.+.+...+.||+.|.
T Consensus        29 q~~~y~~dil~tvy~yn~~~~~dll~elk~llnkls~k~pel~pe~f~~~yp~vse   84 (95)
T tr|G8M0B5|G8M0   29 QNLDYYFDILVTVYEYNYDSEKDLLGELKDLINKLSEKIPELTPEQFKEKYPNISE   84 (95)
Confidence            35677778888888888888877766665 55588999998999999999999874


No 234
>tr|A0A0J6GXZ2|A0A0J6GXZ2_9BACI Aspartate--tRNA ligase (Fragment) OS=Bacillus paralicheniformis GN=aspS PE=4 SV=1
Probab=34.91  E-value=75  Score=30.05  Aligned_cols=29  Identities=38%  Similarity=0.413  Sum_probs=24.8  Template_Neff=3.129

Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      -|.-=|-..|-.||.|+|.....|.||-|
T Consensus        45 aPPHGGiA~GlDRlvmllag~~siRd~IA   73 (109)
T tr|A0A0J6GXZ2|   45 APPHGGIALGLDRLVMLLAGEKSIRDVIA   73 (109)
Confidence            46666777799999999999999999875


No 235
>tr|C1EBJ0|C1EBJ0_MICSR Uncharacterized protein OS=Micromonas sp. (strain RCC299 / NOUM17) GN=MICPUN_108852 PE=4 SV=1
Probab=34.60  E-value=93  Score=33.82  Aligned_cols=56  Identities=29%  Similarity=0.426  Sum_probs=41.3  Template_Neff=1.000

Q d12asa_         238 PVLEDAFELSS-MGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQS  295 (327)
Q Consensus       238 ~~l~~a~ElSS-mGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqS  295 (327)
                      |.-.-.|||+| .|.|.|++.|.--+..  .+.---..|-+.-..|.-|-|||-|.|..
T Consensus       256 ppghfrielasvqglrldreelavairr--aeagwgaafgrsasagsspgtigegmgla  312 (374)
T tr|C1EBJ0|C1EB  256 PPGHFRIELASVQGLRLDREELAVAIRR--AEAGWGAAFGRSASAGSSPGTIGEGMGLA  312 (374)
Confidence            33445678877 5999999998754433  33345567888889999999999999853


No 236
>tr|M2NJW0|M2NJW0_BAUCO Lysine--tRNA ligase OS=Baudoinia compniacensis (strain UAMH 10762) GN=BAUCODRAFT_571914 PE=3 SV=1
Probab=34.38  E-value=94  Score=38.32  Aligned_cols=39  Identities=31%  Similarity=0.393  Sum_probs=31.4  Template_Neff=1.000

Q d12asa_         274 EWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAG  313 (327)
Q Consensus       274 ~~h~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~s  313 (327)
                      .|-.++.-| ||.|=|=|.|..|++|||-..-.|.||-.-
T Consensus      1418 gfltameyg-lpptggwgmgidrmvmfltnnysikevltf 1456 (1487)
T tr|M2NJW0|M2NJ 1418 GFLTAMEYG-LPPTGGWGMGIDRMVMFLTNNYSIKEVLTF 1456 (1487)
Confidence            344443334 899999999999999999999999998653


No 237
>tr|A0A0D2IXC4|A0A0D2IXC4_9CHLO Uncharacterized protein OS=Monoraphidium neglectum GN=MNEG_15366 PE=4 SV=1
Probab=34.37  E-value=94  Score=29.74  Aligned_cols=34  Identities=32%  Similarity=0.499  Sum_probs=28.4  Template_Neff=1.000

Q d12asa_         170 LLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSD  204 (327)
Q Consensus       170 L~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~  204 (327)
                      =-.+||.|||.||.+.| ..+|+++.-|-|..-+|
T Consensus        92 gprrfprltpderqrli-esngcfycrgtghvasd  125 (140)
T tr|A0A0D2IXC4|   92 GPRRFPRLTPDERQRLI-ESNGCFYCRGTGHVASD  125 (140)
Confidence            34689999999998765 67999999999988765


No 238
>tr|X1J2A8|X1J2A8_9ZZZZ Uncharacterized protein (Fragment) OS=marine sediment metagenome GN=S03H2_62063 PE=4 SV=1
Probab=34.03  E-value=96  Score=29.12  Aligned_cols=33  Identities=15%  Similarity=0.171  Sum_probs=28.5  Template_Neff=1.000

Q d12asa_           2 YIAKQRQISFVKSHFSRQLEERLGLIEVQAPIL   34 (327)
Q Consensus         2 ~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPLf   34 (327)
                      +..|.++|+.||+.+...|+-.-||+.|+....
T Consensus         4 vfntdegikaikdqltnllsidsnltpvsnsyw   36 (123)
T tr|X1J2A8|X1J2    4 VFNTDEGIKAIKDQLTNLLSIDSNLTPVSNSYW   36 (123)
Confidence            456889999999999999999999999886544


No 239
>tr|A0A0M2WTF0|A0A0M2WTF0_9BURK Tyrosine recombinase XerC OS=Janthinobacterium sp. KBS0711 GN=xerC_2 PE=4 SV=1
Probab=34.00  E-value=91  Score=35.02  Aligned_cols=77  Identities=25%  Similarity=0.365  Sum_probs=51.3  Template_Neff=1.593

Q d12asa_         216 DDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQS  295 (327)
Q Consensus       216 DDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgqS  295 (327)
                      -.|.+.-+.....|--.|--||..+      +|.||.-      .-...++---|.+..=+++++.|-+|.|+|||-||=
T Consensus       216 ~gw~~~~~g~~a~le~~~~~~~~~~------~~~g~tk------~~~~~tghglraq~aen~al~~~mip~t~gg~~gq~  283 (413)
T tr|A0A0M2WTF0|  216 MGWEYQCSGQIASLEQNIRRWNNMM------TSFGFTK------QDSNCTGHGLRAQFAENHALIADMIPPTMGGSAGQM  283 (413)
Confidence            3577766666666777777777653      5666642      222344555577777788999999999999999994


Q d12asa_         296 RLTMLLLQL  304 (327)
Q Consensus       296 Rl~M~lL~k  304 (327)
                      -=.-.=|.|
T Consensus       284 ~~a~~~~~~  292 (413)
T tr|A0A0M2WTF0|  284 DKADSDVKK  292 (413)
Confidence            433333333


No 240
>tr|E8YJA2|E8YJA2_9BURK Putative uncharacterized protein OS=Burkholderia sp. CCGE1001 GN=BC1001_0009 PE=4 SV=1
Probab=33.80  E-value=97  Score=27.51  Aligned_cols=40  Identities=28%  Similarity=0.394  Sum_probs=35.1  Template_Neff=1.000

Q d12asa_          45 LSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS   84 (327)
Q Consensus        45 LnG~ErpV~F~v~~~~~~~~EIVhSLAKWKR~aL~ky~~~   84 (327)
                      |.|+|..|.+...+..+-..|..|--||---..|.+|+..
T Consensus         7 legsegtvtvqlfpstnvhiemahdaakhlvqllerynwt   46 (87)
T tr|E8YJA2|E8YJ    7 LEGSEGTVTVQLFPSTNVHIEMAHDAAKHLVQLLERYNWT   46 (87)
Confidence            5677888888888877899999999999999999999874


No 241
>tr|I2G2J9|I2G2J9_USTH4 Uncharacterized protein OS=Ustilago hordei (strain Uh4875-4) GN=UHOR_15983 PE=4 SV=1
Probab=33.56  E-value=99  Score=26.91  Aligned_cols=24  Identities=29%  Similarity=0.354  Sum_probs=21.1  Template_Neff=1.000

Q d12asa_          45 LSGAEKAVQVKVKALPDAQFEVVH   68 (327)
Q Consensus        45 LnG~ErpV~F~v~~~~~~~~EIVh   68 (327)
                      |-|+|.||.++-.+.++...||+|
T Consensus        32 ligtetpvrydsdnndndvsevmh   55 (77)
T tr|I2G2J9|I2G2   32 LIGTETPVRYDSDNNDNDVSEVMH   55 (77)
Confidence            679999999999888788888888


No 242
>tr|T0CYU0|T0CYU0_9BACL Uncharacterized protein OS=Alicyclobacillus acidoterrestris ATCC 49025 GN=N007_14340 PE=4 SV=1
Probab=32.91  E-value=1e+02  Score=35.08  Aligned_cols=62  Identities=29%  Similarity=0.373  Sum_probs=51.2  Template_Neff=1.000

Q d12asa_         173 RYPDLDAKG----RERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFEL  246 (327)
Q Consensus       173 ~YP~Ltpke----RE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~El  246 (327)
                      -.|+-|.|+    .|+..+-.|||+|.+.=-+.|-+|+.--..|.|.|-+.|.            .-|.|+.+|.+-|
T Consensus       428 flpdktakqvqwmqenqltisygaffymnergelisgkvgtafasdvdfyrtm------------sqwcpvngctvmi  493 (576)
T tr|T0CYU0|T0CY  428 FLPDKTAKQVQWMQENQLTISYGAFFYMNERGELISGKVGTAFASDVDFYRTM------------SQWCPVNGCTVMI  493 (576)
Confidence            458888886    5889999999999999999999999999999999988743            3577777776544


No 243
>tr|A0A0E0LIU8|A0A0E0LIU8_ORYPU Uncharacterized protein OS=Oryza punctata PE=4 SV=1
Probab=32.77  E-value=1e+02  Score=33.10  Aligned_cols=66  Identities=26%  Similarity=0.443  Sum_probs=40.6  Template_Neff=1.000

Q d12asa_         177 LDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDIL-VWNPVLEDAFELSSMGIRVDA  255 (327)
Q Consensus       177 LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDil-v~n~~l~~a~ElSSmGirVd~  255 (327)
                      ++.+||.-..--|.--|+++     |.||...      .--|.|.-+   ...-|||. ||+|...-.+|++.||+|..+
T Consensus       252 mpyrerlwrlegestevvvv-----lddgrvv------fwvwdtgfs---ratagdimrvydpatggqmevaamgvrpcr  317 (337)
T tr|A0A0E0LIU8|  252 MPYRERLWRLEGESTEVVVV-----LDDGRVV------FWVWDTGFS---RATAGDIMRVYDPATGGQMEVAAMGVRPCR  317 (337)
Confidence            45667766665554444332     3344321      223442211   23568886 999999999999999999765


Q d12asa_         256 D  256 (327)
Q Consensus       256 ~  256 (327)
                      .
T Consensus       318 r  318 (337)
T tr|A0A0E0LIU8|  318 R  318 (337)
Confidence            3


No 244
>tr|A0A0C3C5V9|A0A0C3C5V9_HEBCY Uncharacterized protein OS=Hebeloma cylindrosporum h7 GN=M413DRAFT_29574 PE=4 SV=1
Probab=32.70  E-value=1e+02  Score=34.27  Aligned_cols=53  Identities=32%  Similarity=0.484  Sum_probs=42.5  Template_Neff=1.000

Q d12asa_         154 LAPFLPDQIHFVHSQELLSRYPDLDAKG--------RERAIAKD-LGAVFLVGIGGKLSDGHR  207 (327)
Q Consensus       154 l~~~Lp~~I~FI~sqeL~~~YP~Ltpke--------RE~~i~ke-~gAvFi~gIG~~L~~G~~  207 (327)
                      +...||-.+..+++-+|.+..|-|+..|        ||+.+++| .|||-+++|..- +.|.|
T Consensus       169 laevlpihlelfttgdladlipflpsseyvpsaqfsrensvtredlgavalvqiatc-kggrp  230 (463)
T tr|A0A0C3C5V9|  169 LAEVLPIHLELFTTGDLADLIPFLPSSEYVPSAQFSRENSVTREDLGAVALVQIATC-KGGRP  230 (463)
Confidence            3457888888999999999999998775        89999987 799999999643 24444


No 245
>tr|A0A061IL71|A0A061IL71_CRIGR Membrane protein (Fragment) OS=Cricetulus griseus GN=H671_2g5923 PE=4 SV=1
Probab=32.03  E-value=1.1e+02  Score=30.68  Aligned_cols=37  Identities=32%  Similarity=0.567  Sum_probs=31.3  Template_Neff=1.000

Q d12asa_         203 SDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLE  241 (327)
Q Consensus       203 ~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~  241 (327)
                      +||.||  |.|..|--.|....|--.|.|||++-.|..+
T Consensus        77 angdph--rkpqvdtvqtsgahgdlslkgdiyitspaqg  113 (190)
T tr|A0A061IL71|   77 ANGDPH--RKPQVDTVQTSGAHGDLSLKGDIYITSPAQG  113 (190)
Confidence            588888  7899998888888888889999999888644


No 246
>tr|A0A0P4VRZ8|A0A0P4VRZ8_9EUCA Uncharacterized protein OS=Scylla olivacea PE=4 SV=1
Probab=31.72  E-value=1.1e+02  Score=33.45  Aligned_cols=29  Identities=28%  Similarity=0.504  Sum_probs=24.1  Template_Neff=1.000

Q d12asa_          71 AKWKRQTLGQHDFSAGEGLYTHMKALRPDE  100 (327)
Q Consensus        71 AKWKR~aL~ky~~~~geGiyTdMnAIRrDE  100 (327)
                      .-|+|-|-..|=. ..-||||.||.|.|--
T Consensus       146 plwrrdatghylc-nacglytkmnginrpl  174 (388)
T tr|A0A0P4VRZ8|  146 PLWRRDATGHYLC-NACGLYTKMNGINRPL  174 (388)
Confidence            3699999998854 5689999999998764


No 247
>tr|A0A0G0CFV7|A0A0G0CFV7_9BACT Uncharacterized protein OS=Parcubacteria bacterium GW2011_GWA2_33_14 GN=UR31_C0016G0012 PE=4 SV=1
Probab=31.16  E-value=1.2e+02  Score=29.11  Aligned_cols=27  Identities=33%  Similarity=0.393  Sum_probs=22.9  Template_Neff=1.000

Q d12asa_         183 ERAIAKDLGAVFLVGIGGKLSDGHRHD  209 (327)
Q Consensus       183 E~~i~ke~gAvFi~gIG~~L~~G~~Hd  209 (327)
                      .-++||+.|+||-+.|.++++|-..-|
T Consensus        93 lleiakkkgvvfainiakkmanpyild  119 (136)
T tr|A0A0G0CFV7|   93 LLEIAKKKGVVFAINIAKKMANPYILD  119 (136)
Confidence            457999999999999999998866544


No 248
>tr|T1BE53|T1BE53_9ZZZZ Aspartyl-tRNA synthetase (Fragment) OS=mine drainage metagenome GN=B2A_01004 PE=4 SV=1
Probab=31.06  E-value=86  Score=35.54  Aligned_cols=117  Identities=27%  Similarity=0.346  Sum_probs=71.1  Template_Neff=4.300

Q d12asa_         173 RYPDLDAKGRERAIAKDLGAVFLVGIGG-KLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGI  251 (327)
Q Consensus       173 ~YP~LtpkeRE~~i~ke~gAvFi~gIG~-~L~~G~~Hd~RapDYDDW~t~~~~~~~gLNGDilv~n~~l~~a~ElSSmGi  251 (327)
                      -||-+.+.|-+.++.-.|..+=.-.=+. .+-+-.|...||--||                     =|++ -.||+.=.|
T Consensus       215 DFPlFe~~ee~~r~~a~HHPFTaP~~ed~~~l~tdP~~vra~aYD---------------------lVlN-G~EiGGGSi  272 (376)
T tr|T1BE53|T1BE  215 DFPLFEPDEEEGRWEAAHHPFTAPKPEDIELLDTDPLKVRAQAYD---------------------LVLN-GVEIGGGSI  272 (376)
Confidence            3677777777666666665543222211 0112245555554443                     2232 368999999


Q d12asa_         252 RVDADTLKHQL-ALTGD-EDRLELEWH---QALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       252 rVd~~~L~~Ql-~~~~~-~~r~~~~~h---~~ll~~~lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      |+-...|.+.. +..+- .+..+..|-   .++.-| -|.--|=.+|--||+|+|....-|.+|-|
T Consensus       273 RIH~~elQ~~vf~~lg~~~eea~~kFG~LLeAl~yG-aPPHGGiA~GlDRlvmLL~g~~sIRdVIA  337 (376)
T tr|T1BE53|T1BE  273 RIHDAELQEKVFELLGISEEEAQEKFGFLLEALKYG-APPHGGIALGLDRLVMLLAGAESIRDVIA  337 (376)
Confidence            99877776653 22222 222223333   344456 56666777899999999999999999976


No 249
>tr|Q14KG0|Q14KG0_SPICI Putative uncharacterized protein OS=Spiroplasma citri GN=SPICINP12_014 PE=4 SV=1
Probab=30.42  E-value=1.2e+02  Score=29.95  Aligned_cols=29  Identities=24%  Similarity=0.299  Sum_probs=23.0  Template_Neff=1.174

Q d12asa_         164 FVHSQELLSRYPDLDAKGRERAIAKDLGA  192 (327)
Q Consensus       164 FI~sqeL~~~YP~LtpkeRE~~i~ke~gA  192 (327)
                      .|++===-+.|.+||..||.++||.+|..
T Consensus        81 iissli~n~~y~dlskqerinkiaeeyni  109 (161)
T tr|Q14KG0|Q14K   81 IISSLIDNEKYNDLSKQERINKIAEEYNI  109 (161)
Confidence            34444445689999999999999999974


No 250
>tr|R6KA16|R6KA16_9FIRM Uncharacterized protein OS=Eubacterium sp. CAG:252 GN=BN564_01583 PE=4 SV=1
Probab=30.37  E-value=1.2e+02  Score=27.35  Aligned_cols=29  Identities=28%  Similarity=0.438  Sum_probs=22.9  Template_Neff=1.000

Q d12asa_         221 PSELGHAGLNGDILVWNPVLEDAFELSSM  249 (327)
Q Consensus       221 ~~~~~~~gLNGDilv~n~~l~~a~ElSSm  249 (327)
                      +...+..|-|.|++||.-.+-.|+.+-|.
T Consensus        58 eattpktgdnsdlyvwyiilvaaisvksh   86 (94)
T tr|R6KA16|R6KA   58 EATTPKTGDNSDLYVWYIILVAAISVKSH   86 (94)
Confidence            34447789999999999999888876653


No 251
>tr|M7C583|M7C583_CHEMY Uncharacterized protein OS=Chelonia mydas GN=UY3_07148 PE=4 SV=1
Probab=30.11  E-value=1.2e+02  Score=27.76  Aligned_cols=31  Identities=32%  Similarity=0.557  Sum_probs=28.0  Template_Neff=1.000

Q d12asa_         216 DDWSTPSELGHAGLNGDILVWNPVLEDAFEL  246 (327)
Q Consensus       216 DDW~t~~~~~~~gLNGDilv~n~~l~~a~El  246 (327)
                      .-|+|++.-....||.|+++-||.+...||-
T Consensus        50 erwktqttcfrgalnndlllsnpplktviet   80 (104)
T tr|M7C583|M7C5   50 ERWKTQTTCFRGALNNDLLLSNPPLKTVIET   80 (104)
Confidence            5799999888889999999999999998884


No 252
>tr|V6LT09|V6LT09_9EUKA Uncharacterized protein OS=Spironucleus salmonicida GN=SS50377_12380 PE=4 SV=1
Probab=29.64  E-value=1.3e+02  Score=32.47  Aligned_cols=35  Identities=23%  Similarity=0.313  Sum_probs=29.8  Template_Neff=1.000

Q d12asa_           9 ISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD   43 (327)
Q Consensus         9 I~~iK~~F~~~L~~~LnL~rVsaPLfv~~~sGlND   43 (327)
                      |-..--.++..|+.+|.|.+.|-||+.++++|-|-
T Consensus        47 iyqadcaletalcaelqlaelsfplvmrketgana   81 (331)
T tr|V6LT09|V6LT   47 IYQADCALETALCAELQLAELSFPLVMRKETGANA   81 (331)
Confidence            33444567889999999999999999999999884


No 253
>tr|A8SD86|A8SD86_9FIRM Uncharacterized protein OS=Faecalibacterium prausnitzii M21/2 GN=FAEPRAM212_02149 PE=4 SV=1
Probab=29.61  E-value=1.3e+02  Score=28.20  Aligned_cols=35  Identities=34%  Similarity=0.574  Sum_probs=29.5  Template_Neff=1.000

Q d12asa_         270 RLELEWHQALLRGEMPQT--IGGGIGQSRLTMLLLQL  304 (327)
Q Consensus       270 r~~~~~h~~ll~~~lP~T--IgGGIgqSRl~M~lL~k  304 (327)
                      ...+-||.+++.|+|-+-  .+||+|-.-||+-+|+.
T Consensus        66 figyifhnallagklvqitvysggvgacalclqmlqn  102 (117)
T tr|A8SD86|A8SD   66 FIGYIFHNALLAGKLVQITVYSGGVGACALCLQMLQN  102 (117)
Confidence            356779999999999763  58999999999988875


No 254
>tr|A5BX84|A5BX84_VITVI Putative uncharacterized protein OS=Vitis vinifera GN=VITISV_011037 PE=4 SV=1
Probab=29.04  E-value=1.3e+02  Score=32.57  Aligned_cols=28  Identities=36%  Similarity=0.462  Sum_probs=25.9  Template_Neff=1.000

Q d12asa_         226 HAGLNGDILVWNPVLEDAFELSSMGIRV  253 (327)
Q Consensus       226 ~~gLNGDilv~n~~l~~a~ElSSmGirV  253 (327)
                      -+|.|-|++.||---..|.+|||.|.|-
T Consensus       298 chgmnddvvaynygersahclssagfry  325 (350)
T tr|A5BX84|A5BX  298 CHGMNDDVVAYNYGERSAHCLSSAGFRY  325 (350)
Confidence            5789999999999999999999999985


No 255
>tr|R7E9K0|R7E9K0_9BACE Uncharacterized protein OS=Bacteroides uniformis CAG:3 GN=BN594_01607 PE=4 SV=1
Probab=28.67  E-value=1.3e+02  Score=26.95  Aligned_cols=34  Identities=18%  Similarity=0.453  Sum_probs=30.1  Template_Neff=1.556

Q d12asa_         168 QELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGK  201 (327)
Q Consensus       168 qeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~  201 (327)
                      .+.++.+...+|.|..+-.|++||.|..+.-|.+
T Consensus        45 n~ve~~f~~v~~~e~~~~~c~~yg~~y~~r~gek   78 (81)
T tr|R7E9K0|R7E9   45 NQVESYFSSVKPNEKRRYMCRDYGIIYYFRRGEK   78 (81)
Confidence            5677888899999999999999999999988875


No 256
>tr|K8XZW9|K8XZW9_RHOOP Uncharacterized protein (Fragment) OS=Rhodococcus opacus M213 GN=WSS_A10202 PE=4 SV=1
Probab=27.88  E-value=1.4e+02  Score=28.93  Aligned_cols=67  Identities=31%  Similarity=0.470  Sum_probs=44.3  Template_Neff=1.000

Q d12asa_         192 AVFLVGIGGKLSDGHRHDVRA----PDYDDWSTPSELGHAGLNGDILV----WN--PVLEDAFELSSM-GIRVDADTLKH  260 (327)
Q Consensus       192 AvFi~gIG~~L~~G~~Hd~Ra----pDYDDW~t~~~~~~~gLNGDilv----~n--~~l~~a~ElSSm-GirVd~~~L~~  260 (327)
                      ..+|.|+|.--..|.|||-..    ...-||+|.+.+|+.   |-+-|    |.  .-...+-|-|.. -|||-+..|+.
T Consensus        21 vtvllglgagtaggaphdweevaecesggdwttdtgngfy---gglqftpetwkanggtgmaseasqaeqirvaenvlkt   97 (148)
T tr|K8XZW9|K8XZ   21 VTVLLGLGAGTAGGAPHDWEEVAECESGGDWTTDTGNGFY---GGLQFTPETWKANGGTGMASEASQAEQIRVAENVLKT   97 (148)
Confidence            346788999999999999653    456799988887543   33332    32  223445555543 58888888877


Q d12asa_         261 Q  261 (327)
Q Consensus       261 Q  261 (327)
                      |
T Consensus        98 q   98 (148)
T tr|K8XZW9|K8XZ   98 Q   98 (148)
Confidence            6


No 257
>tr|T1DKS2|T1DKS2_ANOAQ Putative secreted protein (Fragment) OS=Anopheles aquasalis PE=2 SV=1
Probab=27.83  E-value=1.3e+02  Score=28.91  Aligned_cols=31  Identities=35%  Similarity=0.543  Sum_probs=23.7  Template_Neff=1.761

Q d12asa_         261 QLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQ  294 (327)
Q Consensus       261 Ql~~~~~~~r~~~~~h~~ll~~~lP~TIgGGIgq  294 (327)
                      ||+..+.   .-+.|.+.+++..-|+.+.||||.
T Consensus        57 ~lk~l~d---eyly~nr~l~n~~~~~~~~gg~~~   87 (123)
T tr|T1DKS2|T1DK   57 QLKQLND---EYLYFNRLLLNGGSPLELAGGISR   87 (123)
Confidence            4444444   345699999999999999999985


No 258
>tr|L1JXI5|L1JXI5_GUITH Uncharacterized protein OS=Guillardia theta CCMP2712 GN=GUITHDRAFT_133463 PE=4 SV=1
Probab=27.35  E-value=1.5e+02  Score=34.65  Aligned_cols=51  Identities=31%  Similarity=0.485  Sum_probs=36.6  Template_Neff=1.000

Q d12asa_         226 HAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQ  277 (327)
Q Consensus       226 ~~gLNGDilv~n~~l~~a~ElSSmGirVd~~~L~~Ql~~~~~~~r~~~~~h~  277 (327)
                      ..|--||++|-.|+.+||| |-.||||--+.+|-.-........+.-..|.+
T Consensus       615 argepgdvlvhgpvtqcaf-lhamgireraqalysasqeeeqqqkiikqyer  665 (697)
T tr|L1JXI5|L1JX  615 ARGEPGDVLVHGPVTQCAF-LHAMGIRERAQALYSASQEEEQQQKIIKQYER  665 (697)
Confidence            4567899999999999998 78899999888886544433333344344443


No 259
>tr|A0A0R1QJV5|A0A0R1QJV5_9LACO Uncharacterized protein OS=Lactobacillus manihotivorans DSM 13343 = JCM 12514 GN=FD01_GL000860 PE=4 SV=1
Probab=27.07  E-value=1.5e+02  Score=25.86  Aligned_cols=37  Identities=30%  Similarity=0.460  Sum_probs=0.0  Template_Neff=1.343

Q d12asa_         178 DAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDY  215 (327)
Q Consensus       178 tpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDY  215 (327)
                      +|.|----+|.+..+.|++ ||....||..|--.+||.
T Consensus        11 dpsedfayfaq~~p~~fly-~gc~~~dg~~hphhspdf   47 (70)
T tr|A0A0R1QJV5|   11 DPSEDFAYFAQSRPSSFLY-IGCDKDDGQDHPHHSPDF   47 (70)


No 260
>tr|M9M860|M9M860_PSEA3 Uncharacterized protein OS=Pseudozyma antarctica (strain T-34) GN=PANT_26d00073 PE=4 SV=1
Probab=25.99  E-value=1.6e+02  Score=36.68  Aligned_cols=170  Identities=22%  Similarity=0.319  Sum_probs=90.8  Template_Neff=1.000

Q d12asa_          64 FEVVHSLAKWKRQTLGQHDFS-AGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGE---RQFSTLKSTVEAIWA  139 (327)
Q Consensus        64 ~EIVhSLAKWKR~aL~ky~~~-~geGiyTdMnAIRrDE~~ld~~HSiyVDQWDWEkvI~~~d---Rnl~~Lk~tV~kIy~  139 (327)
                      .|-+|-.-|=.-.++--||-. .|.||-     +| |. -+.-+--+--|.-.|...|...-   .|+--+.+..+|.--
T Consensus      1183 leklqdyfknetfailgygsqghgqgln-----lr-dn-glnvivgvrkdgaswkeaiedgwqpgknlfpieeaaekati 1255 (1514)
T tr|M9M860|M9M8 1183 LEKLQDYFKNETFAILGYGSQGHGQGLN-----LR-DN-GLNVIVGVRKDGASWKEAIEDGWQPGKNLFPIEEAAEKATI 1255 (1514)
Confidence            344444444445555556654 566653     34 32 35566666677778888777543   455555555555443


Q d12asa_         140 GIKAT-EAAVSEEFG-LAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDD  217 (327)
Q Consensus       140 al~~t-e~~v~~~y~-l~~~Lp~~I~FI~sqeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDD  217 (327)
                      ++.-. +.....-|+ ++.+.|..-+..-+.-..-.|-+-|--|    ....... ++.     -+.|.-.-.|+-=-  
T Consensus      1256 amlllsdaactsvynsikdkmpnthtlyfshgfqvtyknetgfe----apsdkdi-ilc-----apkgsgrtvrslfl-- 1323 (1514)
T tr|M9M860|M9M8 1256 AMLLLSDAACTSVYNSIKDKMPNTHTLYFSHGFQVTYKNETGFE----APSDKDI-ILC-----APKGSGRTVRSLFL-- 1323 (1514)
Confidence            43333 334445577 7778887655444444444444322111    0000110 100     00111111111100  


Q d12asa_         218 WSTPSELGHAGLNGDILVWNPVLEDAFE-LSSMGIRVDADTLK  259 (327)
Q Consensus       218 W~t~~~~~~~gLNGDilv~n~~l~~a~E-lSSmGirVd~~~L~  259 (327)
                             ...|.||.|-||..+.++++| +|.||+-|.-.-|-
T Consensus      1324 -------esrgingsiavyqdvsgkalekvsamgvavgcgyly 1359 (1514)
T tr|M9M860|M9M8 1324 -------ESRGINGSIAVYQDVSGKALEKVSAMGVAVGCGYLY 1359 (1514)
Confidence                   246799999999999999998 78899988655443


No 261
>tr|A0A0J9SNC8|A0A0J9SNC8_PLAVI Uncharacterized protein (Fragment) OS=Plasmodium vivax Brazil I GN=PVBG_00579 PE=4 SV=1
Probab=25.87  E-value=1.6e+02  Score=27.51  Aligned_cols=24  Identities=25%  Similarity=0.577  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         151 EFGLAPFLPDQIHFVHSQELLSRY  174 (327)
Q Consensus       151 ~y~l~~~Lp~~I~FI~sqeL~~~Y  174 (327)
                      +|.-..+|..+|+||+-++..++|
T Consensus         2 kfnaeqklknkihfiteedfkkmy   25 (115)
T tr|A0A0J9SNC8|    2 KFNAEQKLKNKIHFITEEDFKKMY   25 (115)


No 262
>tr|F8NHN4|F8NHN4_SERL9 Putative uncharacterized protein (Fragment) OS=Serpula lacrymans var. lacrymans (strain S7.9) GN=SERLADRAFT_456645 PE=4 SV=1
Probab=25.73  E-value=1.7e+02  Score=31.57  Aligned_cols=29  Identities=28%  Similarity=0.503  Sum_probs=25.6  Template_Neff=1.000

Q d12asa_         297 LTMLLLQLPHIGQVQAGVWPAAVRESVPS  325 (327)
Q Consensus       297 l~M~lL~k~HIgEVq~svW~~~~~~~~~~  325 (327)
                      -|||||-+-|+|....-+=|.+++...++
T Consensus         4 acmfllcrlhvgsletkitpkevrnalen   32 (317)
T tr|F8NHN4|F8NH    4 ACMFLLCRLHVGSLETKITPKEVRNALEN   32 (317)
Confidence            49999999999999999999998877654


No 263
>tr|L8WKL1|L8WKL1_THACA CK1/CK1/CK1-G protein kinase OS=Thanatephorus cucumeris (strain AG1-IA) GN=AG1IA_08637 PE=3 SV=1
Probab=25.27  E-value=1.7e+02  Score=34.64  Aligned_cols=37  Identities=35%  Similarity=0.610  Sum_probs=28.3  Template_Neff=1.000

Q d12asa_         183 ERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSEL  224 (327)
Q Consensus       183 E~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~  224 (327)
                      ..++..|.|.|||.. |..    ..||+|.|.|..++|.+..
T Consensus        50 deaageeegsvflcd-gya----mshdgrppeyspystqtrp   86 (797)
T tr|L8WKL1|L8WK   50 DEAAGEEEGSVFLCD-GYA----MSHDGRPPEYSPYSTQTRP   86 (797)
Confidence            346677889999863 222    2799999999999988765


No 264
>tr|H9J8U5|H9J8U5_BOMMO Uncharacterized protein OS=Bombyx mori PE=4 SV=1
Probab=25.21  E-value=1.7e+02  Score=28.43  Aligned_cols=22  Identities=27%  Similarity=0.537  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         175 PDLDAKGRERAIAKDLGAVFLV  196 (327)
Q Consensus       175 P~LtpkeRE~~i~ke~gAvFi~  196 (327)
                      |..++||...+.+||||+.+++
T Consensus        33 ppvsakeklkkavkeygstviv   54 (147)
T tr|H9J8U5|H9J8   33 PPVSAKEKLKKAVKEYGSTVIV   54 (147)


No 265
>tr|Q22DK8|Q22DK8_TETTS Uncharacterized protein OS=Tetrahymena thermophila (strain SB210) GN=TTHERM_00942710 PE=4 SV=1
Probab=23.22  E-value=2e+02  Score=32.66  Aligned_cols=32  Identities=28%  Similarity=0.400  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_           1 AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPI   33 (327)
Q Consensus         1 ~~~eTq~aI~~iK~~F~~~L~~~LnL~rVsaPL   33 (327)
                      |++|.|..|+.+|+.|++ ||.-..--+|.|||
T Consensus       225 sfkelqssikdlknqfqq-lcqyvkeneiqapl  256 (497)
T tr|Q22DK8|Q22D  225 SFKELQSSIKDLKNQFQQ-LCQYVKENEIQAPL  256 (497)


No 266
>tr|C0JP56|C0JP56_RHYSE MAT1-1-3 (Fragment) OS=Rhynchosporium secalis GN=MAT1-1 PE=4 SV=1
Probab=22.56  E-value=2.1e+02  Score=29.26  Aligned_cols=28  Identities=29%  Similarity=0.506  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         208 HDVRAPDYDDWSTPSELGHAGLNGDILV  235 (327)
Q Consensus       208 Hd~RapDYDDW~t~~~~~~~gLNGDilv  235 (327)
                      |-..+.-.|.|..++++...|||.|+.|
T Consensus       131 hvnfaaafdswdaetsstmnglnsdlsf  158 (204)
T tr|C0JP56|C0JP  131 HVNFAAAFDSWDAETSSTMNGLNSDLSF  158 (204)


No 267
>tr|W9DTD3|W9DTD3_METTI PDK repeat-containing protein OS=Methanolobus tindarius DSM 2278 GN=MettiDRAFT_0355 PE=4 SV=1
Probab=22.35  E-value=2.1e+02  Score=34.10  Aligned_cols=76  Identities=21%  Similarity=0.396  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         213 PDYDDWSTPSELGHAGLNGDILVWNPVLEDA--FELS---SMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQT  287 (327)
Q Consensus       213 pDYDDW~t~~~~~~~gLNGDilv~n~~l~~a--~ElS---SmGirVd~~~L~~Ql~~~~~~~r~~~~~h~~ll~~~lP~T  287 (327)
                      .|.|.-+-.+...+..|+|.++-||..--.|  .|+|   |-|+..|.+++.--.+.....-..-..|...+.+|++|+|
T Consensus       189 adndkltysttaefgtlygnvfewnttgieagkyeisflvsdgystdsetimitikkeandvfpvvnfttnvtngkiplt  268 (822)
T tr|W9DTD3|W9DT  189 ADNDKLTYSTTAEFGTLYGNVFEWNTTGIEAGKYEISFLVSDGYSTDSETIMITIKKEANDVFPVVNFTTNVTNGKIPLT  268 (822)


Q d12asa_         288 I  288 (327)
Q Consensus       288 I  288 (327)
                      +
T Consensus       269 v  269 (822)
T tr|W9DTD3|W9DT  269 V  269 (822)


No 268
>tr|A0A0C1D7Y4|A0A0C1D7Y4_9FLAO Uncharacterized protein OS=Chryseobacterium jeonii GN=OA86_05190 PE=4 SV=1
Probab=22.17  E-value=2.1e+02  Score=25.81  Aligned_cols=35  Identities=14%  Similarity=0.283  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         119 VMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFG  153 (327)
Q Consensus       119 vI~~~dRnl~~Lk~tV~kIy~al~~te~~v~~~y~  153 (327)
                      +..+...+++.||.|+.--|+.-+++|-.+.-.||
T Consensus        37 lfsedkisleklkktlkltynsykeiefyiayhyp   71 (90)
T tr|A0A0C1D7Y4|   37 LFSEDKISLEKLKKTLKLTYNSYKEIEFYIAYHYP   71 (90)


No 269
>tr|R9NAC4|R9NAC4_9FIRM Uncharacterized protein OS=Dorea sp. 5-2 GN=C817_04031 PE=4 SV=1
Probab=22.15  E-value=2.1e+02  Score=27.16  Aligned_cols=30  Identities=27%  Similarity=0.360  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_          12 VKSHFSRQLEER---LGLIEVQAPILSRVGDGT   41 (327)
Q Consensus        12 iK~~F~~~L~~~---LnL~rVsaPLfv~~~sGl   41 (327)
                      ||..|++.|++.   +.|.+-+|-||+...+||
T Consensus         3 iktyfekslcesaglfslqkkqadlfllsqtgl   35 (124)
T tr|R9NAC4|R9NA    3 IKTYFEKSLCESAGLFSLQKKQADLFLLSQTGL   35 (124)


No 270
>tr|G7UUD7|G7UUD7_PSEUP Elongation factor P--(R)-beta-lysine ligase OS=Pseudoxanthomonas spadix (strain BD-a59) GN=DSC_07885 PE=3 SV=1
Probab=21.97  E-value=2.1e+02  Score=33.17  Aligned_cols=29  Identities=38%  Similarity=0.589  Sum_probs=0.0  Template_Neff=1.324

Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIGQVQA  312 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIgEVq~  312 (327)
                      +|-+.|=.+|--||.|.+|+...|.||-+
T Consensus       495 ~p~cagva~g~~rllmaml~t~~ia~vla  523 (529)
T tr|G7UUD7|G7UU  495 MPHCAGVAVGMDRLLMAMLGTSRIAEVLA  523 (529)


No 271
>tr|A0A061D193|A0A061D193_BABBI Uncharacterized protein OS=Babesia bigemina GN=BBBOND_0102010 PE=4 SV=1
Probab=21.22  E-value=2.3e+02  Score=27.02  Aligned_cols=49  Identities=27%  Similarity=0.338  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         181 GRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGL  229 (327)
Q Consensus       181 eRE~~i~ke~gAvFi~gIG~~L~~G~~Hd~RapDYDDW~t~~~~~~~gL  229 (327)
                      |.-.+-+|-.-.||+|.|=-+|..+.....|+|-.|.++||+.....|.
T Consensus        25 eqadkdakvqmivfmitiifplagnnkrgqrspaqdafstptaatlggi   73 (125)
T tr|A0A061D193|   25 EQADKDAKVQMIVFMITIIFPLAGNNKRGQRSPAQDAFSTPTAATLGGI   73 (125)


No 272
>tr|A0A0C9N3C5|A0A0C9N3C5_9FUNG Lysine--tRNA ligase OS=Mucor ambiguus GN=MAM1_0243c08615 PE=3 SV=1
Probab=20.23  E-value=2.5e+02  Score=33.27  Aligned_cols=26  Identities=35%  Similarity=0.621  Sum_probs=0.0  Template_Neff=1.000

Q d12asa_         284 MPQTIGGGIGQSRLTMLLLQLPHIGQ  309 (327)
Q Consensus       284 lP~TIgGGIgqSRl~M~lL~k~HIgE  309 (327)
                      ||.|.|=|.|..|++-++-.-.||.|
T Consensus       507 lpptagwgmgidrvvqlmtgathire  532 (738)
T tr|A0A0C9N3C5|  507 LPPTAGWGMGIDRVVQLMTGATHIRE  532 (738)


No 273
>tr|R5MDR2|R5MDR2_9BACE Uncharacterized protein OS=Bacteroides intestinalis CAG:564 GN=BN711_00177 PE=4 SV=1
Probab=20.06  E-value=2e+02  Score=31.09  Aligned_cols=34  Identities=18%  Similarity=0.374  Sum_probs=0.0  Template_Neff=3.700

Q d12asa_         168 QELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGK  201 (327)
Q Consensus       168 qeL~~~YP~LtpkeRE~~i~ke~gAvFi~gIG~~  201 (327)
                      .+.++.|-+++|.|...+-|++||.|..+.=|.+
T Consensus       197 NqVE~~f~~~pp~e~R~~~Cr~yGvVYyyR~~Ek  230 (237)
T tr|R5MDR2|R5MD  197 NQVEAKFQDEPPNEKRLEKCRDYGVVYYYRRGEK  230 (237)

"""


atab_text = """\
>3m4p_A Ehasnrs, asparaginyl-tRNA synthetase, putative; aminoacyl-tRNA synthetase, tRNA ligase, AARS, translation, ATP-binding, nucleotide-binding; HET: 4AD; 2.83A {Entamoeba histolytica} PDB: 3m4q_A
    i     j  score     SS  probab  dssp
    9   159   0.10   0.00  0.3565     H
   10   160   0.35   0.00  0.3835     H
   11   161   0.15   0.00  0.4059     H
   12   162   1.39   0.00  0.4281     H
   13   163  -0.19   0.00  0.4436     H
>1h4v_B Histidyl-tRNA synthetase; class IIA aminoacyl-tRNA synthetase, ATP + L-histidine tRNA(His)-> AMP + PPI + L-histidyl-tRNA(His); 2.4A {Thermus thermophilus} SCOP: c.51.1.1 d.104.1.1 PDB: 1ady_A* 1adj_A 4rdx_A*
missing dssp
    i     j  score     SS  probab
    6    17   0.23   0.00  0.4388
    7    18   0.55   0.00  0.5639
    8    19   0.37   0.00  0.6517
    9    20   0.64   0.00  0.7241
   10    21   0.43   0.00  0.7749
   11    22   0.37   0.00  0.8221
"""

if __name__ == '__main__':
    unittest.main()
