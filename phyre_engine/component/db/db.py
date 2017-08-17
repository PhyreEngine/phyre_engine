from phyre_engine.component import Component
import phyre_engine.tools.pdb as pdb
import phyre_engine.tools.conformation
from enum import Enum
import pathlib
import urllib.request
import Bio.PDB
import logging
import collections
import json

log = lambda: logging.getLogger(__name__)

class StructureType(Enum):
    PDB = "pdb"
    MMCIF = "cif"

class StructureRetriever(Component):
    """
    Downloads a structure from the RCSB, saving it with the usual naming
    convention in a specified base directory. The PDB code specified by the
    ``PDB`` key of the pipeline state will be

    For example, if the PDB file ``4hhb.pdb`` is to be downloaded, it will be
    saved in a directory underneath the base directory as ``hh/4hhb.pdb``. The
    type of file to be downloaded can be set using the ``struc_type`` parameter
    of the constructor.


    :param .StructureType struc_type: Type of file to download.
    :param str base_dir: Base directory in which to save PDB files.
    """

    REQUIRED = ["PDB"]
    ADDS = []
    REMOVES = []

    #: URL template from which files are retrieved.
    URL = "https://files.rcsb.org/download/{PDB}.{type}.gz"

    def __init__(self, struc_type, base_dir="."):
        self.struc_type = StructureType(struc_type)
        self.base_dir = pathlib.Path(base_dir)

    def run(self, data, config=None, pipeline=None):
        """Run component."""
        pdb_id = self.get_vals(data)
        url = self.URL.format(
            PDB=pdb_id.lower(),
            type=self.struc_type.value)

        path = pdb.pdb_path(
            pdb_id,
            ".{}.gz".format(self.struc_type.value),
            base_dir=self.base_dir)
        path.parent.mkdir(parents=True, exist_ok=True)

        urllib.request.urlretrieve(url, str(path))
        return data

class ChainPDBBuilder(Component):
    """
    Extract a chain from an MMCIF file and save it as a PDB file. The PDB ID
    and chain ID are given by the ``PDB`` and ``chain`` fields of the pipeline
    state. The ``structure`` field, containing the file path of the new PDB file
    is added.

    Sometimes it is useful to preserve the arcane mappings used by the authors
    of PDB structures. For our purposes, we are often mapping between sequence
    and structure, so it is useful to treat PDB files as arrays of residues. To
    preserve the mappings, then, we write a JSON-encoded array of residue IDs
    to the header of the PDB file under ``REMARK 149`` (the decimal ascii codes
    of "P" and "E" added together).

    If the ``conf_sel`` parameter is not supplied (i.e. the default value of
    ``None`` is used), then the selectors
    :py:class:`phyre_engine.tool.conformation.PopulationMutationSelection` and
    :py:class:`phyre_engine.tool.conformation.PopulationMicroHetSelector` are
    applied in turn. To disable filtering, pass an empty list.
    """

    REQUIRED = ["PDB", "chain"]
    ADDS = ["structure"]
    REMOVES  = []

    class MissingSourceError(RuntimeError):
        """Raised when an MMCIF file is unexpectedly missing."""

        _ERROR = "Could not find source file for PDB ID {}"

        def __init__(self, pdb_id):
            super().__init__(self._ERROR.format(pdb_id))
            self.pdb_id = pdb_id


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
        self.overwrite = overwrite
        if conf_sel is not None:
            self.conf_sel = conf_sel
        else:
            self.conf_sel = [
                phyre_engine.tools.conformation.PopulationMutationSelector(),
                phyre_engine.tools.conformation.PopulationMicroHetSelector()
            ]

    def run(self, data, config=None, pipeline=None):
        """Run the component."""
        pdb_id, chain = self.get_vals(data)

        parser = Bio.PDB.MMCIFParser()
        pdbio = Bio.PDB.PDBIO()

        source_file = pdb.find_pdb(pdb_id, base_dir=self.mmcif_dir)
        if source_file is None:
            log().error(
                "Could not find MMCIF file '%s' in '%s'",
                pdb_id, self.mmcif_dir)
            raise self.MissingSourceError(pdb_id)

        pdb_file = pdb.pdb_path(pdb_id, ".pdb", chain, self.pdb_dir)
        pdb_file.parent.mkdir(parents=True, exist_ok=True)

        if not pdb_file.exists() or self.overwrite:
            with pdb.open_pdb(source_file) as pdb_in:
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
        data["structure"] = str(pdb_file)
        return data

class PDBSequence(Component):
    """
    Read a sequence from the ATOM records of a PDB structure.

    The pipeline state must have the ``structure`` key defined, pointing to a
    PDB file from which ATOM records will be parsed. This component adds the
    ``sequence`` key to the pipeline state. The value of the ``sequence`` field
    will be a Python string consisting of single-letter amino acids.
    """

    ADDS = ["sequence"]
    REMOVES = []
    REQUIRED = ["structure"]

    def run(self, data, config=None, pipeline=None):
        structure_path = self.get_vals(data)

        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure("", structure_path)
        chain = list(structure[0].get_chains())[0]
        data["sequence"], _ = pdb.atom_seq(chain)
        return data

class Reduce(Component):
    """
    Group elements of the pipeline state with identical fields together.

    This component will examine each element of the list specified by the
    ``item_list`` parameter. Within that list, elements are grouped by the field
    or fields specified by the ``compare_field`` parameter. Identical elements
    are replaced in the ``item_list`` by a single element.

    A ``reduction`` key will be added to the pipeline state containing the
    grouped elements. For example, the ``reduction`` key might look like the
    following after grouping on the ``type`` field:

    .. code-block:: none

        [
            [
                {"type": 1, "PDB": "1abc", "chain": "B", ...},
                {"type": 1, "PDB": "1xyz", "chain": "C", ...}
                # ...
            ],
            [
                {"type": 2, "PDB": "1foo", "chain": "A", ...},
                {"type": 2, "PDB": "2bar", "chain": "X", ...},
                # ...
            ]
            # ...
        ]

    :param str item_list: Name of the field in the pipeline state containing the
        elements to be compared.
    :param tuple(str) compare_fields: Tuple of fields to be compared.
    """
    ADDS = ["reduction"]
    REMOVES = []
    REQUIRED = []

    CONFIG_SECTION = "reduce"

    def __init__(self, item_list, compare_fields):
        self.item_list = item_list
        self.compare_fields = compare_fields

    def run(self, data, config=None, pipeline=None):
        """Prune elements by identical field."""

        # Indexed by tuple of field values
        identical = collections.defaultdict(lambda: [])

        for element in data[self.item_list]:
            fields = tuple(element[f] for f in self.compare_fields)
            identical[fields].append(element)

        log().info(
            ("Reduced '%s' with %d elements to %d non-identical elements by "
             "fields %s"),
            self.item_list, len(data[self.item_list]),
            len(identical), self.compare_fields)

        # Generate the new list of representatives and the mapping.
        data[self.item_list] = [ident[0] for ident in identical.values()]
        data["reduction"] = list(identical.values())
        return data

class Expand(Component):
    """
    Expand a list of elements using the definitions in the ``reduction`` key of
    the pipeline state.

    This component is intended to complement :py:class:`.Reduce` by expanding a
    list of templates to contain all templates with an identical sequence. Any
    fields in the representative that were not in the clustered elements are
    copied directly to the cluster members.

    .. seealso::
       :py:class:`.Reduce`
           For constructor parameters.
    """
    REQUIRED = ["reduction"]
    ADDS = []
    REMOVES = ["reduction"]

    CONFIG_SECTION = "reduce"

    def __init__(self, item_list, compare_fields):
        self.item_list = item_list
        self.compare_fields = compare_fields

    @staticmethod
    def _update(orig, other):
        """
        Update the template 'other' with all the parameters from 'orig' that
        are not present in 'other'.
        """
        for key, value in orig.items():
            if key not in other:
                other[key] = value
        return other

    def _reduction_dict(self, reduction):
        """Build a dict of cluster members indexed by 'compare_field'."""
        reduction_dict = {}

        for clus_mem in reduction:
            field_index = tuple(clus_mem[0][f] for f in self.compare_fields)
            reduction_dict[field_index] = clus_mem
        return reduction_dict

    def run(self, data, config=None, pipeline=None):
        """Expand a key with grouped sequences."""
        reduction = self.get_vals(data)

        reduction_dict = self._reduction_dict(reduction)
        extra = []
        for representative in data[self.item_list]:
            field_index = tuple(representative[f] for f in self.compare_fields)
            if field_index in reduction_dict:
                clus_members = reduction_dict[field_index]
                for clus_member in clus_members[1:]:
                    expanded = self._update(representative, clus_member)
                    extra.append(expanded)
        log().info(
            "Adding %d templates to the %d already present.",
            len(extra), len(data[self.item_list]))
        data[self.item_list].extend(extra)
        return data
