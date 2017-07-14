from phyre_engine.component import Component
import copy
from enum import Enum
import gzip
import pathlib
import Bio.SeqUtils
import json
import sys
import urllib.request
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.PDB
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.PDBExceptions import PDBConstructionException
import contextlib

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

    This component will read the chain of each representative from the
    corresponding MMCIF file, write it to a PDB file, and store the sequence of
    the ATOM records in that PDB file.

    Sometimes it is useful to preserve the arcane mappings used by the authors
    of PDB structures. For our purposes, we are often mapping between sequence
    and structure, so it is useful to treat PDB files as arrays of residues. To
    preserve the mappings, then, we write a JSON-encoded array of residue IDs
    to a "map" file.
    """

    REQUIRED = ["templates"]
    ADDS     = []
    REMOVES  = []

    def __init__(
        self, mmcif_dir, pdb_dir, map_dir,
        mutation_selector, microhet_selector=None):
        """Initialise new component.

        :param str mmcif_dir: Base directory of the MMCIF archive.
        :param str pdb_dir: Base directory in which to store PDB files.
        :param str map_dir: Base directory in which sequence map files will be
            saved.
        :param phyre_engine.component.db.conformation.MutationSelector mutation_selector:
            Subclass of `phyre_engine.component.db.conformation.MutationSelector`
            used for selecting which mutation of a template should be used.
        :param phyre_engine.component.db.conformation.MicroConformationSelector microhet_selector:
            Subclass of
            `phyre_engine.component.db.conformation.MicroConformationSelector`
            used to pick a single conformation. All conformations will be kept
            if this is not supplied.

        """
        self.mmcif_dir = pathlib.Path(mmcif_dir)
        self.pdb_dir = pathlib.Path(pdb_dir)
        self.map_dir = pathlib.Path(map_dir)
        self.mutation_selector = mutation_selector
        self.microhet_selector = microhet_selector

    def run(self, data, config=None, pipeline=None):
        """Run the component."""
        templates = self.get_vals(data)

        # Create a new list rather than modifying the existing list in place.
        # This is so that we can handle failures gracefully, if you can call
        # squawking an exception and then ignoring the template graceful.
        new_templates = []
        for template in templates:
            pdb_id = template["PDB"].lower()
            chain = template["chain"]
            middle = pdb_id[1:3].lower()

            # Create dir in which to store pdb file
            pdb_dir = self.pdb_dir / middle
            map_dir = self.map_dir / middle
            pdb_dir.mkdir(exist_ok=True)
            map_dir.mkdir(exist_ok=True)

            mmcif_file = self.mmcif_dir / middle / "{}.cif.gz".format(pdb_id)

            # If the output files already exist, we can just read from them.
            try:
                (pdb_files, map_files) = self.find_template_files(pdb_id, chain)
                if pdb_files and map_files:
                    for pdb_file, map_file in zip(pdb_files, map_files):
                        new_template = copy.copy(template)
                        self.read_chain(
                            pdb_id, chain, new_template, pdb_file, map_file)
                        new_templates.append(new_template)
                else:
                    extracted_templates = self.extract_chain(
                        pdb_id, chain, template,
                        mmcif_file)
                    new_templates.extend(extracted_templates)
            except PDBConstructionException as e:
                # TODO: When we start to use a logging framework, log this error
                # properly.
                print(
                    "Error parsing chain {} of {}".format(chain, pdb_id),
                    file=sys.stderr)
            except Exception as e:
                raise ChainExtractionError(pdb_id, chain) from e

        data["templates"] = new_templates
        return data

    def find_template_files(self, pdb_id, chain_id):
        pdb_id = pdb_id.lower()
        middle = pdb_id[1:3]
        pdb_dir = self.pdb_dir / middle
        map_dir = self.map_dir / middle

        # First, check if the correct files exist for a single conformation.
        # Next, check if multiple conformations exist.
        # Otherwise, return None
        pdb_file = pdb_dir / "{}_{}.pdb".format(pdb_id, chain_id)
        map_file = map_dir / "{}_{}.json".format(pdb_id, chain_id)
        if pdb_file.exists() and map_file.exists():
            return ([pdb_file], [map_file])
        else:
            pdb_files = sorted(
                list(pdb_dir.glob("{}_{}-*.pdb")),
                key=lambda p: p.name)
            map_files = sorted(
                list(map_dir.glob("{}_{}-*.json")),
                key=lambda p: p.name)
            if pdb_files and map_files:
                return (pdb_files, map_files)
                #TODO: Throw error if lists don't match
            else:
                return (None, None)

    def read_chain(
            self, pdb_id, chain_id, template,
            pdb_file, map_file):
        """
        Read a previously-extracted chain.

        Sets the ``sequence``, ``structure`` and ``map`` fields of the
        ``template`` argument. Modifies ``template`` in place.

        :param str pdb_id: PDB ID.
        :param dict template: Dictionary describing template.
        :param pathlib.Path pdb_file: Path object pointing to the PDB file into
            which the chain will be saved.
        :param pathlib.Path map_file: Path object poniting to the file into
            which the residue mapping will be saved.
        """
        with pdb_file.open("r") as pdb_fh:
            structure = PDBParser().get_structure(pdb_id, pdb_fh)
            model = next(structure.get_models())
            chain = model[' ']

        # Get sequence of the residues.
        pdb_seq = "".join([
            Bio.SeqUtils.seq1(r.get_resname())
            for r in chain
        ])

        # Build a Bio.PDB.SeqRecord object containing this sequence.
        seq_name = "{}_{}".format(pdb_id, chain_id)
        bio_seq = SeqRecord(
                id=seq_name,
                description="",
                seq=Seq(pdb_seq)
        )
        template["name"] = seq_name
        template["sequence"] = bio_seq
        template["structure"] = pdb_file
        template["map"] = map_file


    def extract_chain(
            self, pdb_id, chain, template,
            mmcif_file):
        """
        Extract a chain from an MMCIF file.

        This method will extract the chain, save a PDB file containing that
        chain, and write a map file containing the map between residue index and
        author-assigned residue IDs.

        :return: A list of templates, shallow-copied from the original template,
            and with "sequence", "structure" and "map" keys added.

        :param str pdb_id: PDB ID.
        :param str chain: PDB chain.
        :param str template: Dictionary describing template.
        :param pathlib.Path mmcif_file: Path object pointing to the MMCIF from
            which the chain  will be extracted.
        """
        mmcif_parser = MMCIFParser()
        pdbio = Bio.PDB.PDBIO()

        pdb_id = pdb_id.lower()
        middle = pdb_id[1:3]
        pdb_dir = self.pdb_dir / middle
        map_dir = self.map_dir / middle

        with gzip.open(str(mmcif_file), "rt") as mmcif_fh:
            structure = mmcif_parser.get_structure(pdb_id, mmcif_fh)
        struc_model = next(structure.get_models())

        templates = []
        mutations = self.mutation_selector.select(struc_model[chain])
        for index, conformation in enumerate(mutations):
            if self.microhet_selector:
                conformation = self.microhet_selector.select(conformation)
            struc_chain, res_map = self.sanitise_chain(conformation)

            if len(mutations) > 1:
                pdb_file = pdb_dir / "{}_{}-{}.pdb".format(pdb_id, chain, index+1)
                map_file = map_dir / "{}_{}-{}.json".format(pdb_id, chain, index+1)
            else:
                pdb_file = pdb_dir / "{}_{}.pdb".format(pdb_id, chain)
                map_file = map_dir / "{}_{}.json".format(pdb_id, chain)

            # Build mapping between sequence index and residue ID. This is just
            # an array of residue IDs that can be indexed in the same way as
            # everything else.
            with map_file.open("w") as map_fh:
                json.dump(res_map, map_fh)

            # Write PDB file.
            pdbio.set_structure(struc_chain)
            pdbio.save(str(pdb_file))

            # Get sequence of the residues.
            pdb_seq = "".join([
                Bio.SeqUtils.seq1(r.get_resname())
                for r in struc_chain
            ])

            # Build a Bio.PDB.SeqRecord object containing this sequence.
            bio_seq = SeqRecord(
                    id="{}_{}".format(pdb_id, chain),
                    description="",
                    seq=Seq(pdb_seq)
            )
            new_template = copy.copy(template)
            new_template["sequence"] = bio_seq
            new_template["structure"] = pdb_file
            new_template["map"] = map_file
            templates.append(new_template)
        return templates

    def sanitise_chain(self, chain, new_id=" "):
        """Strip insertion codes and disordered atoms from a chain.

        Given a `Bio.PDB` structure (`Bio.PDB.Chain.Chain`), renumber residues
        from 1, stripping out insertion codes. At the same time, we remove
        disordered atoms and residues. For now, we just keep the first
        conformation that occurs.

        :param `Bio.PDB.Chain` chain: structure to sanitise.
        :param str new_id: ID of the new chain.

        :return: A tuple containing:

            1. A new `Bio.PDB.Chain` object consisting of sanitised residues.
            2. An array containing the author-assigned IDs of the sanitised residues.
        """
        mapping = []
        sanitised_chain = Chain(new_id)
        res_index = 1
        for res in chain:
            # Discard HETATMs
            if res.get_id()[0] != ' ':
                continue

            sanitised_res = Residue(
                (' ', res_index, ' '),
                res.get_resname(),
                res.get_segid())
            mapping.append(res.get_id())

            for atom in res:
                sanitised_res.add(atom.copy())
            sanitised_chain.add(sanitised_res)
            res_index += 1
        return sanitised_chain, mapping

class ChainExtractionError(Exception):
    """
    Exeception thrown when an error is encountered extracting a chain from an
    MMCIF file.
    """

    ERR_MSG = "Error extracting chain {} from structure {}"
    def __init__(self, structure, chain):
        super().__init__(self.ERR_MSG.format(chain, structure))

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
