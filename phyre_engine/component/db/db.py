from phyre_engine.component import Component
import phyre_engine.tools.pdb as pdb
from enum import Enum
import pathlib
import urllib.request
import Bio.PDB
import logging

log = lambda: logging.getLogger(__name__)

class StructureType(Enum):
    PDB = "pdb"
    MMCIF = "cif"

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

            path = pdb.pdb_path(
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

            source_file = pdb.find_pdb(pdb_id, base_dir=self.mmcif_dir)
            if source_file is None:
                log().error(
                    "Could not find MMCIF file '%s' in '%s'",
                    pdb_id, self.mmcif_dir)
                continue

            pdb_file = pdb.pdb_path(pdb_id, ".pdb", chain, self.pdb_dir)
            pdb_file.parent.mkdir(parents=True, exist_ok=True)

            if not pdb_file.exists() or self.overwrite:
                with source_file.open("r") as pdb_in:
                    structure = parser.get_structure(
                        "{}_{}".format(pdb_id, chain),
                        pdb_in)
                    chain = structure[0][chain]

                for selector in self.conf_sel:
                    chain = selector.select(chain)
                mapping, chain = pdb.renumber(chain, ".")

                with pdb_file.open("w") as pdb_out:
                    pdb.write_json_remark(pdb_out, mapping, 149)
                    pdbio.set_structure(chain)
                    pdbio.save(pdb_out)
            template["structure"] = str(pdb_file)
            new_templates.append(template)

        data["templates"] = new_templates
        return data

class PDBSequence(Component):
    """
    Read a sequence from a list of ATOM records for each template by parsing the
    PDB chain from each template.

    The templates must have the ``structure`` key defined, pointing to a PDB
    file from which ATOM records will be parsed. This component adds the
    ``sequence`` key to each element of the ``templates`` list in the pipeline
    state. The value of the ``sequence`` field will be a
    :py:class:`Bio.SeqIO.SeqRecord.SeqRecord` object.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = ["templates"]

    def run(self, data, config=None, pipeline=None):
        templates = self.get_vals(data)

        parser = Bio.PDB.PDBParser()
        for template in templates:
            structure = parser.get_structure("", template["structure"])
            chain = list(structure[0].get_chains())[0]
            template["sequence"], _ = pdb.atom_seq(chain)
        return data

