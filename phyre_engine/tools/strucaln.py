'''
Contains structural alignment algorithms.
'''
import re
import subprocess
import sys
import phyre_engine.tools.external as ext
from phyre_engine.tools.util import NamedTuple

class TMAlign(object):
    '''
    Aligns two structures with TM-align.

    :ivar invert: (Default: `False`) Whether to invert the order of arguments
    when calling the TMalign executable. This only affects the class internally:
    aligning *a* to *b* will give the same results in the same order whether or
    not this variable is true. This only affects the raw output of TMalign and
    the order in which structures appear in any superpositions. It may be useful
    to set this variable when aligning lots of structures to a reference
    structure because TMalign always rotates the first structure onto the
    second, so superpositions produced in this way will all be consistent with
    the reference structure.
    :vartype invert: bool
    '''

    class Result:
        '''
        Represents an alignment between two structures produced by TMalign.

        :ivar tm: Tuple of TM-scores (normalised by chain 1 and 2)
        :ivar length: Tuple of chain lengths.
        :ivar rmsd: RMSD over aligned region.
        :ivar seqid: Sequence identity over aligned region.
        :ivar sequences: Tuple of sequences (chain 1, map, chain 2)
        '''

        tm1_regex = re.compile(r"^TM-score=\s*([0-9.]+).*Chain_1", re.MULTILINE)
        tm2_regex = re.compile(r"^TM-score=\s*([0-9.]+).*Chain_2", re.MULTILINE)
        len1_regex = re.compile(r"^Length of Chain_1:\s*(\d+)", re.MULTILINE)
        len2_regex = re.compile(r"^Length of Chain_2:\s*(\d+)", re.MULTILINE)
        rmsd_regex = re.compile(r"RMSD=\s*([0-9.]+)", re.MULTILINE)
        seqid_regex = re.compile(r"Seq_ID=n_identical/n_aligned=\s*([0-9.]+)", re.MULTILINE)
        sup_suffixes = ["", "_all", "_atm", "_all_atm", "_all_atm_lig"]


        def __init__(self, scores, sequences, superpositions=None):
            '''
            Construct a new Result.

            :param tuple scores: Dictionary of scores produced by TMAlign.
            :param tuple sequences: Tuple of sequences produced by TMAlign.
            :param tuple superpositions: Tuple of superpositions produced by TMAlign.
            '''
            self.scores = scores
            self.sequences = sequences
            self.superpositions = superpositions

            self.tm = scores["TM"]
            self.length = scores["length"]
            self.rmsd = scores["rmsd"]
            self.seqid = scores["seqid"]

        def sup(self, all_regions=False, all_atoms=False, ligands=False):
            '''
            Returns the path to an alignment with the corresponding parameters.
            :param bool all_regions: Include only the aligned regions?
            :param bool all_atoms: Include all atoms or just the C-alpha trace?
            :param bool ligands: Include ligands? Only valud if `all_regions` and
            `all_atoms` are true.
            :return: Path to file containing the superposition.
            '''
            if self.superpositions is None:
                return None

            #0: .
            #1: _all
            #2: _atm
            #3: _all_atm
            #4: _all_atm_lig
            if not all_regions and not all_atoms:
                return self.superpositions[0]
            elif not all_regions and all_atoms:
                return self.superpositions[1]
            elif all_regions and not all_atoms:
                if ligands:
                    return self.superpositions[4]
                else:
                    return self.superpositions[2]
            else:
                return self.superpositions[3]

        @classmethod
        def parse_str(cls, out, invert=False, sup_prefix=None):
            '''
            Parse TMalign output into a :py:class:`.TMAlign.Result` object.
            :param str out:  Standard output of TMalign.
            :param bool invert: Should we invert scores and sequences?
            :param str sup_prefix: Prefix of superposition files.
            :return: New alignment.
            '''
            tm1 = cls.tm1_regex.search(out).group(1)
            tm2 = cls.tm2_regex.search(out).group(1)
            len1 = cls.len1_regex.search(out).group(1)
            len2 = cls.len2_regex.search(out).group(1)
            rmsd = cls.rmsd_regex.search(out).group(1)
            seqid = cls.seqid_regex.search(out).group(1)
            scores = {
                "TM": (float(tm1), float(tm2)),
                "rmsd": float(rmsd),
                "length": (int(len1), int(len2)),
                "seqid": float(seqid)
            }
            #Sequence lines are the final three non-empty lines
            seqs = tuple([l for l in out.split("\n") if len(l) > 0][-3:])
            if invert:
                scores["TM"] = (scores["TM"][1], scores["TM"][0])
                scores["length"] = (scores["length"][1], scores["length"][0])
                seqs = (seqs[2], seqs[1], seqs[0])

            sups = None
            if sup_prefix:
                sups = tuple([sup_prefix + suffix for suffix in cls.sup_suffixes])
            return cls(scores, seqs, sups)



    def __init__(self, invert=False, bin_dir=None, executable="TMalign"):
        '''
        Initialise a new TMalign object, ready to align pairs of structures.

        :param bool invert: Whether to invert arguments to TMalign.
        :param str bin_dir: Path to TMalign executable.
        :param str executable: Name of executable.
        '''
        self.bin_dir = bin_dir
        self.executable = executable
        self.invert = invert

    def align(self, file1, file2, superpos=None):
        '''
        To a single alignment between two PDB files.

        :param str file1: Path of first file.
        :param str file2: Path of second file.
        :param str superpos: Prefix of superposition files to generate, or
        `None` if no superposition is to be written. Files will be generated
        with the standard suffixes generated by TMalign with the ``-o`` option.

        :return: Alignment object
        :rtype: phyre_engine.tools.strucaln.TMAlign.Result
        '''

        cmd_line = ext.ExternalTool()(
            (self.bin_dir, self.executable),
            positional=[file2, file1] if self.invert else [file1, file2],
            options={"o": superpos} if superpos is not None else None)
        tmalign_out = subprocess.run(
            cmd_line, stdout=subprocess.PIPE, check=True,
            universal_newlines=True)
        return self.Result.parse_str(
            tmalign_out.stdout, self.invert, superpos)

class TMScore:
    """
    Align two structures with identical sequences using
    `TMscore <https://zhanglab.ccmb.med.umich.edu/TM-score/>`_.
    """

    class Result(NamedTuple):
        """
        The results of an alignment by TMscore.

        All scores are normalised by the length of the second structure.

        :ivar tm: TM-scores (float between 0 and 1) of the alignment.

        :ivar maxsub: MaxSub score (between 0 and 1) for the alignment.

        :ivar gdt_ts: Global Distance Test (Total Score), between 0 and 1.

        :ivar gdt_ha: Global Distance Test (High Accuracy), between 0 and 1.

        :ivar sequences: 2-tuple of sequences, showing the alignment
            between each structure.
        """
        FIELDS = "tm maxsub gdt_ts gdt_ha sequences superpositions"
        _TMSCORE_REGEX = re.compile(r"^TM-score\s+=\s*([0-9.]+)",
                                    re.MULTILINE)
        _MAXSUB_REGEX = re.compile(r"^MaxSub-score=\s*([0-9.]+)",
                                   re.MULTILINE)
        _GDT_TS_REGEX = re.compile(r"^GDT-TS-score=\s*([0-9.]+)",
                                   re.MULTILINE)
        _GDT_HA_REGEX = re.compile(r"^GDT-HA-score=\s*([0-9.]+)",
                                   re.MULTILINE)
        _SEQUENCE_LINES = (-4, -2)

        _SUP_SUFFIXES = ["", "_atm"]

        @classmethod
        def parse_str(cls, str_result, sup_prefix=None):
            """
            Parse the standard output of TMscore.

            The parameter `str_result` should look like this:

            .. code-block:: none

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
                Structure2: 2n77_B/04-  Length=   23 (by which all scores are normalized)
                Number of residues in common=   23
                RMSD of  the common residues=    0.000

                TM-score    = 1.0000  (d0= 0.68)
                MaxSub-score= 1.0000  (d0= 3.50)
                GDT-TS-score= 1.0000 %(d<1)=1.0000 %(d<2)=1.0000 %(d<4)=1.0000 %(d<8)=1.0000
                GDT-HA-score= 1.0000 %(d<0.5)=1.0000 %(d<1)=1.0000 %(d<2)=1.0000 %(d<4)=1.0000

                 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
                 i          t(i)         u(i,1)         u(i,2)         u(i,3)
                 1      0.0000000000   1.0000000000  -0.0000000000   0.0000000000
                 2      0.0000000000   0.0000000000   1.0000000000  -0.0000000000
                 3     -0.0000000000  -0.0000000000   0.0000000000   1.0000000000

                Superposition in the TM-score: Length(d<5.0)= 23  RMSD=  0.00
                (":" denotes the residue pairs of distance < 5.0 Angstrom)
                ETERAAVAIQSQFRKFQKKKAGS
                :::::::::::::::::::::::
                ETERAAVAIQSQFRKFQKKKAGS
                12345678901234567890123

            :return: TMscore results
            :rtype: :py:class:`.TMScore.Result`
            """
            # I think data model magic in NamedTuple is breaking pylint.
            # pylint: disable=no-member

            # Get non-empty lines
            lines = [l for l in str_result.split("\n") if l]
            sequences = tuple([lines[i] for i in cls._SEQUENCE_LINES])

            superpositions = None
            if sup_prefix is not None:
                superpositions = [sup_prefix + suffix
                                  for suffix in cls._SUP_SUFFIXES]

            return cls(
                tm=float(cls._TMSCORE_REGEX.search(str_result).group(1)),
                maxsub=float(cls._MAXSUB_REGEX.search(str_result).group(1)),
                gdt_ts=float(cls._GDT_TS_REGEX.search(str_result).group(1)),
                gdt_ha=float(cls._GDT_HA_REGEX.search(str_result).group(1)),
                sequences=sequences,
                superpositions=superpositions
            )

    def __init__(self, bin_dir=None, executable="TMscore"):
        self.bin_dir = bin_dir
        self.executable = executable

    def align(self, file1, file2, superpos=None):
        """
        Align PDB files `file1` and `file2` using TMscore.

        :param str file1: PDB file 1.

        :param str file2: PDB file 2. Scores are normalised by the length of
            this protein.

        :param str superpos: Prefix used for superposition files. If `None`,
            no superposition files are generated.
        """
        cmd_line = ext.ExternalTool()(
            (self.bin_dir, self.executable),
            positional=[file1, file2],
            options={"o": superpos} if superpos is not None else None)

        tmalign_out = subprocess.run(cmd_line, stdout=subprocess.PIPE,
                                     check=True, universal_newlines=True)
        # pylint: disable=no-member
        return self.Result.parse_str(tmalign_out.stdout, superpos)
