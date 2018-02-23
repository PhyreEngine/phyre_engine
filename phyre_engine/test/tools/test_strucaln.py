import os
from pathlib import Path
import unittest
import tempfile

import phyre_engine.test
import phyre_engine.tools.strucaln as sa
import textwrap

DATA_DIR = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')


class TestTMAlign(unittest.TestCase):
    """Test common methods for TMAlign tools."""

    @phyre_engine.test.requireFields(["bin_dir"], ["tools", "tmalign"])
    def test_runner(self):
        """Try running TMalign"""
        # pylint: disable=unsubscriptable-object
        config = phyre_engine.test.config["tools"]["tmalign"]
        tmalign_dir = config["bin_dir"]
        tmalign_exec = config.get("executable", "TMalign")
        pdb = str(Path(DATA_DIR, "pdb_chains/2a/12as_A.pdb"))
        tmalign = sa.TMAlign(bin_dir=tmalign_dir, executable=tmalign_exec)

        with tempfile.TemporaryDirectory() as tmpdir:
            result = tmalign.align(pdb, pdb, str(Path(tmpdir, "sup")))
        self.assertTupleEqual(result.length, (327, 327), "Parsed length")
        self.assertAlmostEqual(
            result.tm[0], 1.0, 6, "Parsed first TM")
        self.assertAlmostEqual(
            result.tm[1], 1.0, 6, "Parsed second TM")
        self.assertAlmostEqual(
            result.rmsd, 0.0, 3, "Parsed RMSD")
        self.assertAlmostEqual(
            result.seqid, 1.0, 4, "Parsed sequence identity")
        self.assertEqual(
            result.sequences[0],
            'AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL',
            "Parsed query sequence")
        self.assertEqual(
             result.sequences[1],
             ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::',
             "Parsed map sequence")
        self.assertEqual(
            result.sequences[2],
            'AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL',
            "Parsed template sequence")

    def test_parser(self):
        """Make sure we can parse known output."""
        test_output = r'''

 **************************************************************************
 *                        TM-align (Version 20160521)                     *
 * An algorithm for protein structure alignment and comparison            *
 * Based on statistics:                                                   *
 *       0.0 < TM-score < 0.30, random structural similarity              *
 *       0.5 < TM-score < 1.00, in about the same fold                    *
 * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)    *
 * Please email your comments and suggestions to: zhng@umich.edu          *
 **************************************************************************

Name of Chain_1: /home/stefans/Downloads/12asA.pdb                 
Name of Chain_2: /home/stefans/Downloads/4lns.pdb                  
Length of Chain_1:  327 residues
Length of Chain_2:  304 residues

Aligned length=  297, RMSD=   2.31, Seq_ID=n_identical/n_aligned= 0.586
TM-score= 0.85066 (if normalized by length of Chain_1)
TM-score= 0.91243 (if normalized by length of Chain_2)
(You should use TM-score normalized by length of the reference protein)

(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
-----AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGK-LSDGHRHDVRAPDYDDWSTPSELGH-AGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
     :::::::::::::::::::::::::::::::::::::::::::       :::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::...  ..              .::.      .:::::::::::::::::::::::::::::::::........:::::::::::::::::::::::::::::::::::::::::::::::::::: ::: ::
DDGYSSYVLLQEQILTVKRSFSEALEKELNLVEVRAPILFRVGDGTQD-------AVQVPVKAIPNASFEVVHSLAKWKRRTLANYKFAPGHGLYTHMTALRVDD-VLDNIHSVVVDQWDWEMVMKDDQRNLAFLKEVVCKVYAAIRKTELAVCEKYKQKPILPETIQFVHAEHLLLAYPNLTAKEREREIAREYGAVFLIGIGA-VLS--------------SLSS-----LKGLNGDILLYNPTLDDSLEVSSMGIRVNAEALRHQISLTGDDSLLKSEWHQQLLNGEFPQTVGGGIGQSRMVMFMLRKKHIGEVQCSVWPEEIR-KKH-NL

'''
        result = sa.TMAlign.Result.parse_str(test_output, False, "sup")
        self.assertTupleEqual(result.length, (327, 304), "Parsed length")
        self.assertAlmostEqual(
            result.tm[0], 0.85066, 6, "Parsed first TM")
        self.assertAlmostEqual(
            result.tm[1], 0.91243, 6, "Parsed second TM")
        self.assertAlmostEqual(
            result.rmsd, 2.31, 3, "Parsed RMSD")
        self.assertAlmostEqual(
            result.seqid, 0.586, 4, "Parsed sequence identity")
        self.assertEqual(
            result.sequences[0],
            '-----AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGK-LSDGHRHDVRAPDYDDWSTPSELGH-AGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL',
            "Parsed query sequence")
        self.assertEqual(
             result.sequences[1],
             '     :::::::::::::::::::::::::::::::::::::::::::       :::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::...  ..              .::.      .:::::::::::::::::::::::::::::::::........:::::::::::::::::::::::::::::::::::::::::::::::::::: ::: ::',
             "Parsed map sequence")
        self.assertEqual(
            result.sequences[2],
            'DDGYSSYVLLQEQILTVKRSFSEALEKELNLVEVRAPILFRVGDGTQD-------AVQVPVKAIPNASFEVVHSLAKWKRRTLANYKFAPGHGLYTHMTALRVDD-VLDNIHSVVVDQWDWEMVMKDDQRNLAFLKEVVCKVYAAIRKTELAVCEKYKQKPILPETIQFVHAEHLLLAYPNLTAKEREREIAREYGAVFLIGIGA-VLS--------------SLSS-----LKGLNGDILLYNPTLDDSLEVSSMGIRVNAEALRHQISLTGDDSLLKSEWHQQLLNGEFPQTVGGGIGQSRMVMFMLRKKHIGEVQCSVWPEEIR-KKH-NL',
            "Parsed template sequence")

        #Test Inverted
        result = sa.TMAlign.Result.parse_str(test_output, True, "sup")
        self.assertTupleEqual(result.length, (304, 327), "Parsed length")
        self.assertAlmostEqual(
            result.tm[1], 0.85066, 6, "Parsed first TM")
        self.assertAlmostEqual(
            result.tm[0], 0.91243, 6, "Parsed second TM")
        self.assertEqual(
            result.sequences[2],
            '-----AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGK-LSDGHRHDVRAPDYDDWSTPSELGH-AGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL',
            "Parsed query sequence")
        self.assertEqual(
             result.sequences[1],
             '     :::::::::::::::::::::::::::::::::::::::::::       :::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::...  ..              .::.      .:::::::::::::::::::::::::::::::::........:::::::::::::::::::::::::::::::::::::::::::::::::::: ::: ::',
             "Parsed map sequence")
        self.assertEqual(
            result.sequences[0],
            'DDGYSSYVLLQEQILTVKRSFSEALEKELNLVEVRAPILFRVGDGTQD-------AVQVPVKAIPNASFEVVHSLAKWKRRTLANYKFAPGHGLYTHMTALRVDD-VLDNIHSVVVDQWDWEMVMKDDQRNLAFLKEVVCKVYAAIRKTELAVCEKYKQKPILPETIQFVHAEHLLLAYPNLTAKEREREIAREYGAVFLIGIGA-VLS--------------SLSS-----LKGLNGDILLYNPTLDDSLEVSSMGIRVNAEALRHQISLTGDDSLLKSEWHQQLLNGEFPQTVGGGIGQSRMVMFMLRKKHIGEVQCSVWPEEIR-KKH-NL',
            "Parsed template sequence")


class TestTMScore(unittest.TestCase):
    """Test :py:class:`phyre_engine.tools.strucaln.TMScore`."""

    SAMPLE_OUTPUT = textwrap.dedent("""\
         *****************************************************************************
         *                                 TM-SCORE                                  *
         * A scoring function to assess the similarity of protein structures         *
         * Based on statistics:                                                      *
         *       0.0 < TM-score < 0.17, random structural similarity                 *
         *       0.5 < TM-score < 1.00, in about the same fold                       *
         * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
         * For comments, please email to: zhng@umich.edu                             *
         *****************************************************************************

        Structure1: 2n77_B/04-  Length=   23
        Structure2: 2n77_B/44-  Length=   24 (by which all scores are normalized)
        Number of residues in common=   23
        RMSD of  the common residues=    3.001

        TM-score    = 0.5490  (d0= 0.78)
        MaxSub-score= 0.8184  (d0= 3.50)
        GDT-TS-score= 0.8438 %(d<1)=0.7500 %(d<2)=0.8333 %(d<4)=0.8750 %(d<8)=0.9167
        GDT-HA-score= 0.7396 %(d<0.5)=0.5000 %(d<1)=0.7500 %(d<2)=0.8333 %(d<4)=0.8750

         -------- rotation matrix to rotate Chain-1 to Chain-2 ------
         i          t(i)         u(i,1)         u(i,2)         u(i,3)
         1     27.6071335522   0.3353627526  -0.9050856164   0.2614418695
         2     34.6565275227  -0.8719482403  -0.1931273944   0.4498978504
         3     83.5013902368  -0.3567044862  -0.3788427595  -0.8539555451

        Superposition in the TM-score: Length(d<5.0)= 21  RMSD=  0.98
        (":" denotes the residue pairs of distance < 5.0 Angstrom)
        -ETERAAVAIQSQFRKFQKKKAGS
         :::::::::::::::::::::··
        PETERAAVAIQSQFRKFQKKKAGS
        123456789012345678901234

    """)

    def test_parser(self):
        """Test TMScore output parser."""
        # Disable warning for broken class method detection.
        # pylint: disable=no-member
        result = sa.TMScore.Result.parse_str(self.SAMPLE_OUTPUT)
        self.assertAlmostEqual(result.tm, 0.5490, 4)
        self.assertAlmostEqual(result.maxsub, 0.8184, 4)
        self.assertAlmostEqual(result.gdt_ts, 0.8438, 4)
        self.assertAlmostEqual(result.gdt_ha, 0.7396, 4)
        self.assertEqual(result.sequences[0], "-ETERAAVAIQSQFRKFQKKKAGS")
        self.assertEqual(result.sequences[1], "PETERAAVAIQSQFRKFQKKKAGS")

    @phyre_engine.test.requireFields(["bin_dir"], ["tools", "tmscore"])
    def test_alignment(self):
        """Align a PDB file with itself using TM-score."""

        # pylint: disable=unsubscriptable-object
        config = phyre_engine.test.config["tools"]["tmscore"]
        tmscore_dir = config["bin_dir"]
        tmscore_exec = config.get("executable", "TMscore")
        pdb = str(Path(DATA_DIR, "pdb_chains/2a/12as_A.pdb"))

        with tempfile.TemporaryDirectory() as tmpdir:
            tmscore = sa.TMScore(bin_dir=tmscore_dir, executable=tmscore_exec)
            result = tmscore.align(pdb, pdb, str(Path(tmpdir, "sup")))
            self.assertAlmostEqual(result.tm, 1.0, 4)
            self.assertAlmostEqual(result.maxsub, 1.0, 4)
            self.assertAlmostEqual(result.gdt_ts, 1.0, 4)
            self.assertAlmostEqual(result.gdt_ha, 1.0, 4)
            self.assertEqual(
                result.sequences[0], (
                "AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPD"
                "AQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGD"
                "GERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRE"
                "RAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAF"
                "ELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLP"
                "HIGQVQAGVWPAAVRESVPSLL"))
            self.assertEqual(result.sequences[0], result.sequences[1])
            for sup in result.superpositions:
                self.assertTrue(Path(sup).exists())
