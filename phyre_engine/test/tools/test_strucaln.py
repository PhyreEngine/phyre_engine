import os
from pathlib import Path
import unittest
import tempfile

import phyre_engine.test
import phyre_engine.tools.strucaln as sa

DATA_DIR = os.path.join(os.path.dirname(phyre_engine.test.__file__), 'data')


class TestTMAlign(unittest.TestCase):
    """Test common methods for TMAlign tools."""

    @phyre_engine.test.requireFields(["tmalign"], ["tools"])
    def test_runner(self):
        """Try running TMalign"""
        # pylint: disable=unsubscriptable-object
        tmalign_exec = phyre_engine.test.config["tools"]["tmalign"]
        pdb = str(Path(DATA_DIR, "pdb_chains/2a/12as_A.pdb"))
        tmalign = sa.TMAlign(tmalign=tmalign_exec)

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
