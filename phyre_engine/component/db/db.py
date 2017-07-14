from phyre_engine.component import Component
from enum import Enum
import gzip
import pathlib
import json
import urllib.request
import Bio.PDB
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import contextlib
import logging

log = lambda: logging.getLogger(__name__)

class StructureType(Enum):
    PDB = "pdb"
    MMCIF = "cif"

# Used when searching for a PDB file of arbitrary type
_STRUCTURE_SUFFIXES = (".pdb", ".pdb.gz", ".cif", ".cif.gz")

class StructureRetriever(Component):
    """
    Downloads a structure from the RCSB, saving it with the usual naming
    convention in a specified base directory.

    For example, if the PDB file ``4hhb.pdb`` is to be downloaded, it will be
    saved in a directory underneath the base directory as ``hh/4hhb.pdb``. The
    type of file to be downloaded can be set using the ``struc_type`` parameter
    of the constructor.

    The files to be downloaded will be determined by the ``PDB`` key of each
    element in the ``templates`` list in the pipeline state.

    :param .StructureType struc_type: Type of file to download.
    :param str base_dir: Base directory in which to save PDB files.
    """

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    #: URL template from which files are retrieved.
    URL = "https://files.rcsb.org/download/{PDB}.{type}.gz"

    def __init__(self, struc_type, base_dir="."):
        self.struc_type = StructureType(struc_type)
        self.base_dir = pathlib.Path(base_dir)

    def run(self, data, config=None, pipeline=None):
        """Run component."""
        templates = self.get_vals(data)
        for template in templates:
            url = self.URL.format(
                PDB=template["PDB"].lower(),
                type=self.struc_type.value)

            path = pdb_path(
                template["PDB"],
                ".{}.gz".format(self.struc_type.value),
                base_dir=self.base_dir)
            path.parent.mkdir(parents=True, exist_ok=True)

            urllib.request.urlretrieve(url, str(path))
        return data

class ChainPDBBuilder(Component):
    """For each structure, extract that chain to a PDB file file
    and extract the sequence (based on ATOM records) of that structure.

    This component will read the chain of each template from the corresponding
    MMCIF file and write it to a PDB file.

    Sometimes it is useful to preserve the arcane mappings used by the authors
    of PDB structures. For our purposes, we are often mapping between sequence
    and structure, so it is useful to treat PDB files as arrays of residues. To
    preserve the mappings, then, we write a JSON-encoded array of residue IDs
    to the header of the PDB file under ``REMARK 149`` (the decimal ascii codes
    of "P" and "E" added together).
    """

    REQUIRED = ["templates"]
    ADDS     = []
    REMOVES  = []

    def __init__(self, mmcif_dir, pdb_dir, conf_sel=None, overwrite=False):
        """Initialise new component.

        :param str mmcif_dir: Base directory of the MMCIF archive.
        :param str pdb_dir: Base directory in which to store PDB files.
        :param list[phyre_engine.tools.conformation.ConformationSelector] conf_sel:
            List of selectors applied in order to each chain to remove pick a
            single conformation.
        :param bool overwrite: If ``True``, overwrite existing PDB files.
        """
        self.mmcif_dir = pathlib.Path(mmcif_dir)
        self.pdb_dir = pathlib.Path(pdb_dir)
        self.conf_sel = conf_sel if conf_sel is not None else []
        self.overwrite = overwrite

    def run(self, data, config=None, pipeline=None):
        """Run the component."""
        templates = self.get_vals(data)

        parser = Bio.PDB.MMCIFParser()
        pdbio = Bio.PDB.PDBIO()
        # Create a new list rather than modifying the existing list in place.
        # This is so that we can handle failures gracefully, if you can call
        # squawking an exception and then ignoring the template graceful.
        new_templates = []
        for template in templates:
            pdb_id = template["PDB"]
            chain = template["chain"]

            source_file = find_pdb(pdb_id, base_dir=self.mmcif_dir)
            if source_file is None:
                log().error(
                    "Could not find MMCIF file '%s' in '%s'",
                    pdb_id, self.mmcif_dir)
                continue

            pdb_file = pdb_path(pdb_id, ".pdb", chain, self.pdb_dir)
            pdb_file.parent.mkdir(parents=True, exist_ok=True)

            if not pdb_file.exists() or self.overwrite:
                with source_file.open("r") as pdb_in:
                    structure = parser.get_structure(
                        "{}_{}".format(pdb_id, chain),
                        pdb_in)
                    chain = structure[0][chain]

                for selector in self.conf_sel:
                    chain = selector.select(chain)
                mapping, chain = self._renumber(chain, ".")

                with pdb_file.open("w") as pdb_out:
                    json_str = json.dumps(mapping)
                    for map_line in _chunk_string(json_str, 69):
                        print("REMARK 149 " + map_line, file=pdb_out)

                    pdbio.set_structure(chain)
                    pdbio.save(pdb_out)
            template["structure"] = str(pdb_file)
            new_templates.append(template)

        data["templates"] = new_templates
        return data

    def _renumber(self, chain, new_id=" "):
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
                (' ', res_index + 1, ' '),
                res.get_resname(),
                res.get_segid())
            mapping.append((res_index + 1, res.get_id()))

            for atom in res:
                sanitised_res.add(atom.copy())
            sanitised_chain.add(sanitised_res)
        return mapping, sanitised_chain

@contextlib.contextmanager
def open_pdb(path):
    """
    Open a PDB or MMCIF file. This should be used as a context manager, and
    automatically handles decompressing gzipped files.
    """
    try:
        if path.suffix == ".gz":
            stream = gzip.open(str(path), "rt")
            yield stream
        else:
            stream = path.open("r")
            yield stream
    finally:
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

def _chunk_string(string, size):
    """Split string into equal length chunks. Used for dumping JSON data into
    the REMARK section."""
    return (string[i:i + size] for i in range(0, len(string), size))
