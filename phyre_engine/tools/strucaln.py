'''
Contains structural alignment algorithms.
'''
import re
import subprocess
import sys

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



    def __init__(self, tmalign="TMalign", invert=False):
        '''
        Initialise a new TMalign object, ready to align pairs of structures.

        :param str tmalign: Path to TMalign executable.
        :param bool invert: Whether to invert arguments to TMalign.
        '''
        self.tmalign = tmalign
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

        cmd_line = [self.tmalign]
        #Add file to command line
        if self.invert:
            cmd_line.extend([file2, file1])
        else:
            cmd_line.extend([file1, file2])

        #Generate superposition?
        if superpos is not None:
            cmd_line.extend(["-o", superpos])

        tmalign_out = subprocess.run(
            cmd_line, stdout=subprocess.PIPE, check=True)
        stdout = tmalign_out.stdout.decode(sys.stdout.encoding)
        return self.Result.parse_str(
            stdout, self.invert, superpos)
