"""
Various tools for working with PDB files.
"""
import gzip
import pathlib
import Bio.Seq
import Bio.PDB
import Bio.Alphabet.IUPAC
import contextlib
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import json
import re

# Used when searching for a PDB file of arbitrary type
_STRUCTURE_SUFFIXES = (".pdb", ".pdb.gz", ".cif", ".cif.gz")

@contextlib.contextmanager
def open_pdb(path):
    """
    Open a PDB or MMCIF file. This should be used as a context manager, and
    automatically handles decompressing gzipped files.
    """
    path = pathlib.Path(path)
    try:
        stream = None
        if path.suffix == ".gz":
            stream = gzip.open(str(path), "rt")
            yield stream
        else:
            stream = path.open("r")
            yield stream
    finally:
        if stream is not None:
            stream.close()

def find_pdb(pdb_name, *args, suffix_list=_STRUCTURE_SUFFIXES, **kwargs):
    """
    For each suffix, try and find a matching file using :py:func:`.pdb_path`.

    :return: Path to the file or ``None`` if no matching file could be found.
    :rtype: :py:class:`pathlib.Path` object.
    """

    for suffix in suffix_list:
        path = pdb_path(pdb_name, suffix, *args, **kwargs)
        if path.exists():
            return path
    return None

def pdb_path(pdb_name, suffix, chain_id=None, base_dir=None):
    """
    Gives the directory in which to save a PDB or MMCIF file.

    If the ``base_dir`` parameter is supplied, paths are given relative to the
    base directory. If the ``chain_id`` is not supplied, then the path of the
    PDB file ``1xyz.pdb`` will be ``xy/1xyz.pdb``. If the ``chain_id`` is
    supplied, then the PDB ID is a separate trailing subdirectory and the chain
    ID is appended before any file extensions.

    :param str pdb_name: Four-letter PDB identifier.
    :param str suffix: Suffix to use (``.pdb`` is probably useful).
    :param str chain_id: Optional chain ID. If supplied, it is assumed that
        chains are stored in separate files, in a subdirectory named after the
        PDB ID.
    :param str base_dir: Optional base directory.

    :return: The path of the PDB file.
    :rtype: Object from :py:mod:`pathlib`.

    >>> from phyre_engine.component.db.db import pdb_path
    >>> pdb_path("1ABC", ".pdb")
    PosixPath("ab/1abc.pdb")
    >>> pdb_path("1ABC", ".cif.gz", "A")
    PosixPath("ab/1abc/1abc_A.cif.gz")
    >>> pdb_path("1ABC", ".cif.gz", "A", "/cif/dir")
    PosixPath("/cif/dir/ab/1abc/1abc_A.cif.gz")
    """

    middle = pdb_name[1:3].lower()
    pdb_id = pdb_name.lower()

    root = pathlib.Path(base_dir) if base_dir is not None else pathlib.Path(".")
    root = root / middle
    if chain_id is not None:
        root = root / pdb_id / "{}_{}{}".format(pdb_id, chain_id, suffix)
    else:
        root = root / "{}{}".format(pdb_id, suffix)
    return root

def atom_seq(chain):
    """
    Return the amino acid sequence of ``chain``, along with a list of
    :py:class:`Bio.PDB.Residue.Residue` objects corresponding to the chosen
    sequence.

    This function does the following filtering:

    1. Discards all HETATMs.
    2. Discards everything but the standard 20 amino acids.
    3. Discard all residues without a CA atom defined.

    This is intended to be the canonical source of a PDB structure's sequence.

    :param Bio.PDB.Chain.Chain chain: PDB chain.

    :return: 2-Tuple containing a one-letter amino acid sequence and a list of
        :py:class:`Bio.PDB.Residue.Residue` classes.
    :rtype: tuple(str, Bio.PDB.Residue.Residue)
    """

    sequence = ""
    residues = []
    for res in chain:
        if (res.get_id()[0] == ' '
                and Bio.PDB.Polypeptide.is_aa(res, True)
                and "CA" in res):
            sequence += Bio.PDB.Polypeptide.three_to_one(res.get_resname())
            residues.append(res)
    return sequence, residues

def renumber(chain, new_id=" "):
    """
    Renumber a chain from 1, stripping insertion codes.

    :param `Bio.PDB.Chain` chain: structure to sanitise.
    :param str new_id: ID of the new chain.
    :return: A 2-tuple containing the following:

        1. The new :py:class:`Bio.PDB.Chain.Chain` object.
        2. A list of 2-tuples containing the new residue index and old
           residue ID. The new residue index is an integer and the old is
           the 3-tuple returned by :py:meth:`Bio.PDB.Chain.Chain.get_id`.
    """
    mapping = []
    sanitised_chain = Chain(new_id)

    for res_index, res in enumerate(chain):
        sanitised_res = Residue(
            (res.get_id()[0], res_index + 1, ' '),
            res.get_resname(),
            res.get_segid())
        mapping.append((res_index + 1, res.get_id()))

        for atom in res:
            sanitised_res.add(atom.copy())
        sanitised_chain.add(sanitised_res)
    return mapping, sanitised_chain

def _chunk_string(string, size):
    """Split string into equal length chunks. Used for dumping JSON data into
    the REMARK section."""
    return (string[i:i + size] for i in range(0, len(string), size))

def write_json_remark(stream, data, remark_num):
    """
    Dump JSON-formatted data to a stream. Data will be written in a REMARK field
    and wrapped to 80 characters.

    :param stream: File handle to which data will be written.
    :param data: Data that will be converted to JSON and written to the stream.
    :param int remark_num: Numeric ID of the REMARK to write.
    """
    json_str = json.dumps(data)
    for map_line in _chunk_string(json_str, 69):
        print("REMARK {:3d} ".format(remark_num) + map_line, file=stream)

def read_json_remark(stream, remark_num):
    """
    Read and parse JSON data from a numbered REMARK field of a PDB file.
    Arguments are the same as for :py:func:`.write_json_remark`.
    """
    json_str = read_remark(stream, remark_num, "")
    return json.loads(json_str)

def read_remark(stream, remark_num, join="\n"):
    """
    Read all REMARK lines with the given number from a PDB file.

    :param stream: File handle from which to read data.
    :param int remark_num: Number of the remark to read.
    :param str join: String used to join lines.
    """
    remark_lines = []
    search_regex = re.compile("^REMARK +{} (.*)$".format(remark_num))
    for line in stream:
        match = search_regex.match(line)
        if match is not None:
            remark_lines.append(match.group(1))
    return join.join(remark_lines)
