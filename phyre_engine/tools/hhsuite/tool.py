"""Module containing wrappers for running hh-suite tools."""
from phyre_engine.tools.external import ExternalTool
import subprocess
import os

#: Wrapper for the hhblits tool included within the hh-suite package.
#:
#: As well as the standard short-form flags of hhblits (see ``hhblits -help``),
#: the following flags are accepted:
#:
#: :param str database: HHsearch database, equivalent to the ``-d`` flag.
#: :param str input: Input alignment file, equivalent to the ``-i`` flag.
#: :param int iterations: Number of search iterations, equivalent to the
#:     ``-n`` flag.
#: :param float evalue_cutoff: e-value cutoff, equivalent to the ``-e`` flag.
#: :param str output: output file, equivalent to the ``-o`` flag.
#: :param int verbose: Verbosity from 0 (silent) to 2 (verbose) (``-v`` flag).
hhblits = ExternalTool(
    flag_map={
        "database": "d",
        "input": "i",
        "iterations": "n",
        "evalue_cutoff": "e",
        "output": "o",
        "verbose": "v"
    }, long_prefix="-")

#: Wrapper for the hhsearch tool included within the hh-suite package.
#:
#: As well as the standard short-form flags of hhsearch (see ``hhsearch
#: -help``), the following flags are accepted:
#:
#: :param str database: HHsearch database, equivalent to the ``-d`` flag.
#: :param str input: input alignment file, equivalent to the ``-i`` flag.
#: :param float evalue_cutoff: e-value cutoff, equivalent to the ``-e`` flag.
#: :param str output: output file, qquivalent to the ``-o`` flag.
#: :param int verbose: Verbosity from 0 (silent) to 2 (verbose) (``-v`` flag).
hhsearch = ExternalTool(
    flag_map = {
        "database":"d",
        "input": "i",
        "evalue_cutoff": "e",
        "output": "o",
        "verbose": "v"
    }, long_prefix = "-")

#: Wrapper for the hhmake tool included within the hh-suite package.
#:
#: All options to ``hhmake`` may be included in the constructor. The input file
#: (the ``-i`` option) is mandatory. The following aliases are accepted as
#: well as the standard options.
#:
#: :param str input: Query alignment file.
#: :param str output: HMM file to be written (``-o`` flag)
#: :param str append: HMM file to be appended (``-a`` flag)
#: :param int verbose: Increase verbosity from 0 (silent) to 2 (verbose). This
#:     is the ``-v`` flag.
hhmake = ExternalTool(
    flag_map={
        "input": "i",
        "output": "o",
        "append": "a",
        "verbose": "v",
    }, long_prefix="-")

#: Wrapper around hhsuite's cstranslate.
#:
#: The following options are available. Hyphens may be replaced with
#: underscores to make valid python identifiers.
#:
#: .. code-block:: none
#:
#:     infile <file>            Input file with alignment or sequence
#:     outfile <file>           Output file for generated abstract state sequence (def: <infile>.as)
#:     append <file>            Append generated abstract state sequence to this file
#:     informat prf|seq|fas|... Input format: prf, seq, fas, a2m, or a3m (def=auto)
#:     outformat seq|prf        Outformat: abstract state sequence or profile (def=seq)
#:     match-assign [0:100]     Make all FASTA columns with less than X% gaps match columns
#:                              (def: make columns with residue in first sequence match columns)
#:     alphabet <file>          Abstract state alphabet consisting of exactly 219 states (def=off)
#:     context-data <file>      Add context-specific pseudocounts using given context-data (def=off)
#:     pc-admix [0,1]           Pseudocount admix for context-specific pseudocounts (def=0.90)
#:     pc-ali [0,inf[           Constant in pseudocount calculation for alignments (def=12.0)
#:     weight [0,inf[           Weight of abstract state column in emission calculation (def=1000.00)
#:     binary                   Write binary instead of character sequence (def=off)
#:     ffindex                  Read from -i <ffindex>, write to -o <ffindex>; enables openmp if possible (def=off)
cstranslate = ExternalTool()

#: Wrapper around ffindex_build.
#:
#: :param str data: Output data file.
#: :param str index: Output index file.
#: :param bool append: (``-a``) Append files/indexes, also needed for sorting
#:     an already existing ffindex.
#: :param str ffdata_file: (``-d``) A second ffindex data file for
#:     inserting/appending.
#: :param str ffindex_file: (``-i``) A second ffindex index file for
#:     inserting/appending.
#: :param str file_list: (``-f``) File containing a list of file names, one per
#:     line.
#: :param str sort: (``-s``)  Sort index file, so that the index can queried.
#:     Append operations can be done without sorting, but should not be.
#:
ffindex_build = ExternalTool(
    flag_map = {
        "append": "a",
        "ffdata_file": "d",
        "ffindex_file": "i",
        "file_list": "f",
        "sort": "s",
    })


#: Wrapper around ``ffindex_modify``.
#:
#: :param bool sort: (``-s``) Sort index file.
#: :param bool unlink: (``-u``) Unlink entry from index.
#: :param bool version: (``-v``) Print version info and exit.
#: :param str file_list: (``-f``) File containing a list of names to remove,
#:     one per line.
#:
ffindex_modify = ExternalTool(
    flag_map = {
        "file_list": "f",
        "sort": "s",
        "version": "v",
        "unlink": "u",
    })


def run(command, *args, HHLIB=None, **kwargs):
    """
    Convenience function to execute a tool from hh-suite.

    This function simply passes all arguments to :py:func:`subprocess.run`, with
    the exception of the ``env`` parameter. If ``HHLIB`` is passed to this
    function, then the ``HHLIB`` environment variable is overriden.
    """
    env = None
    if HHLIB is not None:
        env = os.environ.copy()
        env["HHLIB"] = HHLIB
    subprocess.run(command, *args, env=env, **kwargs)
