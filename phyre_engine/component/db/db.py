import phyre_engine.component.component
from phyre_engine.component import Component
import phyre_engine.tools.pdb as pdb
import phyre_engine.tools.conformation
from phyre_engine.tools.template import Template
import phyre_engine.logutils
from enum import Enum
import pathlib
import urllib.request
import Bio.PDB
import logging
import collections
import json
from builtins import FileNotFoundError
import shutil
import requests

class StructureType(Enum):
    PDB = "pdb"
    MMCIF = "cif"

class PDBList(Component):
    """
    Reads (after optionally downloading) a list of PDB IDs in the format
    supplied by the RCSB's
    `web service <https://www.rcsb.org/pdb/rest/getCurrent>`_. For each PDB ID a
    dictionary is added to the ``templates`` list. Each element of ``templates``
    will contain the PDB ID in the key ``PDB``.

    :param str file: File containing PDB IDs. If this parameter is not supplied,
        it is downloaded from the RCSB website.
    """
    ADDS = ["templates"]
    REQUIRED = []
    REMOVES = []

    PDB_ID_URL = "https://www.rcsb.org/pdb/rest/getCurrent"

    def __init__(self, file=None):
        self.file = file

    def run(self, data, config=None, pipeline=None):
        """Read a list of PDB IDs."""
        if self.file is not None:
            with open(self.file, "r") as file_in:
                entries = pdb.get_current(file_in.read())
        else:
            entries = pdb.get_current()

        data["templates"] = [{"PDB": entry} for entry in entries]
        return data


class HTTPMap(phyre_engine.component.component.Map):
    """
    Operates identically to :py:class:`phyre_engine.component.component.Map`,
    but copies a `requests.Session
    <http://docs.python-requests.org/en/latest/user/advanced/#session-objects>`_
    object into the ``session`` field of the pipeline state for each list
    element. This can be used to speed up repetitive HTTP requests using HTTP
    Keep-Alive.

    The ``session`` field is automatically cleaned up when this component
    exits.
    """

    def run(self, data, config=None, pipeline=None):
        """Iterate over a list of components with a common HTTP session."""
        session = requests.Session()

        try:
            for item in data[self.field]:
                item["session"] = session
            data = super().run(data, config, pipeline)
        finally:
            for item in data[self.field]:
                del item["session"]
        return data

class StructureRetriever(Component):
    """
    Downloads a structure from the RCSB, saving it with the usual naming
    convention in a specified base directory. The PDB code specified by the
    ``PDB`` key of the pipeline state will be

    For example, if the PDB file ``4hhb.pdb`` is to be downloaded, it will be
    saved in a directory underneath the base directory as ``hh/4hhb.pdb``. The
    type of file to be downloaded can be set using the ``struc_type`` parameter
    of the constructor.

    If it is set, the ``mmcif_dir`` field from the ``foldlib`` section of the
    pipeline configuration is used as the `base_dir` parameter.

    :param phyre_engine.component.db.db.StructureType struc_type: Type of file
        to download.

    :param str base_dir: Base directory in which to save PDB files.

    :param bool overwrite: If `False`, structures that have already been
        downloaded are skipped.
    """

    REQUIRED = ["PDB"]
    ADDS = []
    REMOVES = []

    #: URL template from which files are retrieved.
    URL = "https://files.rcsb.org/download/{PDB}.{type}.gz"

    @classmethod
    def config(cls, params, config):
        """
        Sets `base_dir` to ``foldlib.mmcif_dir`` if present.

        .. csv-table:: Configuration mapping
            :header: "Section", "Field", "Parameter"

            ``foldlib``,   ``mmcif_dir``,   ``base_dir``
                       ,   ``overwrite``,   ``overwrite``
        """
        return config.extract({"foldlib": [
                ("mmcif_dir", "base_dir"),
                "overwrite",
            ]}).merge_params(params)

    def __init__(self, struc_type, base_dir=".", overwrite=False):
        self.struc_type = StructureType(struc_type)
        self.base_dir = pathlib.Path(base_dir)
        self.overwrite = overwrite

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
        if (not self.overwrite) and path.exists():
            return data

        path.parent.mkdir(parents=True, exist_ok=True)

        session = data["session"] if "session" in data else requests.Session()
        with session.get(url, stream=True) as r:
            with path.open("wb") as out_fh:
                shutil.copyfileobj(r.raw, out_fh)
        return data

class ChainPDBBuilder(Component):
    """
    Extract each chain from an MMCIF file and save it as a PDB file. The PDB ID
    must be given by the ``PDB`` field of the pipeline state. This component
    returns a *list* of elements, one for each field. Each element will contain
    the following fields:

    ``chain``
        The chain ID.

    ``structure``
        The file path of the new PDB file.

    ``sequence``
        The canonical sequence of the template, as determined by
        :py:class:`phyre_engine.tools.template.Template`.

    ``original_residues``
        The original residue IDs (as BioPython-style tuples of hetero-flag,
        residue number and insertion code) of the residues in the canonical
        sequence.

    ``canonical_indices``
        Residue IDs (indexed from 1) of the renumbered residues chosen for the
        canonical sequence. Joining ``canonical_indices`` and
        ``original_residues`` will give a mapping between the new residue IDs
        and old residues IDs for all residues in the canonical sequence.

    .. note::

        It would make sense for this component to return a ``template_obj`` key,
        as a :py:class:`phyre_engine.tools.template.Template` object is
        generated internally. However, storing (and pickling) those objects is
        relatively expensive. Since this component may return a list of
        templates, it must be the final step in a
        :py:class:`~phyre_engine.component.component.Map` pipeline, which means
        that large numbers of :py:class:`phyre_engine.tools.template.Template`
        objects would pile up. Instead, use :py:class:`.BuildTemplate` to load
        a template into the ``template_obj`` key.

    The chains written by this component are renumbered consecutively starting
    from 1 without any insertion codes. Unless an extra filter is added using
    the `conf_sel` argument, all residue and atom types are written to the PDB
    file.

    This component generates the "canonical" sequence of a PDB structure. This
    sequence is generated from all residues in the ``ATOM`` records of the PDB
    file that are standard amino acids and have a ``CA`` atom. This sequence is
    written as a single line to ``REMARK 150`` (ASCII "C" + "S" for Canonical
    Sequence).

    Finally, a ``REMARK 156`` ("S" + "I" for Sequence Index) will be written
    containing the *renumbered* residue ID corresponding to the canonical
    sequence.

    For example, we might start with the following original PDB file (with
    everything right of the residue index excised for legibility):

    .. code-block:: none

        ATOM      1 CA A  ALA A 10
        ATOM      2 CA A  GLY A 11
        HETATM    3 CA A  AMP A 11A
        ATOM      4 CB A  ALA A 12

    This PDB file will be renumbered beginning from 1, and insertion codes will
    be stripped. It will then look like this:

    .. code-block:: none

        ATOM      1 CA A  ALA A 1
        ATOM      2 CA A  GLY A 2
        HETATM    3 CA A  AMP A 3
        ATOM      4 CB A  ALA A 4

    Next, the canonical sequence and sequence mapping will be added. Residues 3
    and 4 will not be included in the sequence: residue 3 is a ``HETATM`` and
    residue 4 does not contain a ``CA`` atom. The mapping is onto the
    *renumbered* residue IDs:

    .. code-block:: none

        REMARK 150 AG
        REMARK 156 [1, 2]
        ATOM      1 CA A  ALA A 1
        ATOM      2 CA A  GLY A 2
        HETATM    3 CA A  AMP A 3
        ATOM      4 CB A  ALA A 4

    .. warning::

        This module may write JSON in any (valid) format, not necessarily as it
        is shown in these examples. At the moment, JSON data is actually written
        on a single line. I hope that this doesn't break too many programs.

    .. warning::

        This component returns a *list* of items which will almost certainly
        break any components in the pipeline after this. This component should
        be wrapped in a :py:class:`phyre_engine.component.component.Map`
        component, which will take care of merging the output of this component.

    .. note::

        If the ``conf_sel`` parameter is not supplied (i.e. the default value of
        ``None`` is used), then the selectors
        :py:class:`phyre_engine.tool.conformation.PopulationMutationSelection`
        and
        :py:class:`phyre_engine.tool.conformation.PopulationMicroHetSelector`
        are applied in turn. To disable filtering, pass an empty list.
    """

    REQUIRED = ["PDB"]
    ADDS = [
        "structure", "name", "chain", "sequence", "original_residues",
        "canonical_indices",
    ]
    REMOVES  = []

    class MissingSourceError(RuntimeError):
        """Raised when an MMCIF file is unexpectedly missing."""

        _ERROR = "Could not find source file for PDB ID {}"

        def __init__(self, pdb_id):
            super().__init__(self._ERROR.format(pdb_id))
            self.pdb_id = pdb_id


    def __init__(self, mmcif_dir, chain_dir, conf_sel=None, overwrite=False):
        """Initialise new component.

        :param str mmcif_dir: Base directory of the MMCIF archive.
        :param str chain_dir: Base directory in which to store PDB files.
        :param list[phyre_engine.tools.conformation.ConformationSelector] conf_sel:
            List of selectors applied in order to each chain to remove pick a
            single conformation.
        :param bool overwrite: If ``True``, overwrite existing PDB files.
        """
        self.mmcif_dir = pathlib.Path(mmcif_dir)
        self.chain_dir = pathlib.Path(chain_dir)
        self.overwrite = overwrite
        if conf_sel is not None:
            self.conf_sel = conf_sel
        else:
            self.conf_sel = [
                phyre_engine.tools.conformation.PopulationMutationSelector(),
                phyre_engine.tools.conformation.PopulationMicroHetSelector()
            ]

    @classmethod
    def config(cls, params, config):
        """
        Extract ``foldlib.chain_dir`` and ``foldlib.mmcif_dir``.

        .. csv-table:: Configuration mapping
            :header: "Section", "Field", "Parameter"

            ``foldlib``,   ``chain_dir``,   ``chain_dir``
                       ,   ``mmcif_dir``,   ``mmcif_dir``
        """
        return config.extract(
            {"foldlib": ["chain_dir", "mmcif_dir"]}
        ).merge_params(params)

    def run(self, data, config=None, pipeline=None):
        """Run the component."""
        pdb_id = self.get_vals(data)

        mmcif_parser = Bio.PDB.FastMMCIFParser()
        pdb_parser = Bio.PDB.PDBParser()

        source_file = pdb.find_pdb(pdb_id, base_dir=self.mmcif_dir)
        if source_file is None:
            self.logger.error(
                "Could not find MMCIF file '%s' in '%s'",
                pdb_id, self.mmcif_dir)
            raise self.MissingSourceError(pdb_id)

        results = []
        self.logger.debug(
            "Extracting chains from %s (%s)",
            pdb_id, source_file)

        with pdb.open_pdb(source_file) as pdb_in:
            structure = mmcif_parser.get_structure(pdb_id, pdb_in)
            self.logger.debug(
                "Found %d chains in %s",
                len(list(structure.get_chains())), pdb_id)

            for chain in structure[0]:
                pdb_file = pdb.pdb_path(
                    pdb_id, ".pdb", chain.id, self.chain_dir)
                pdb_file.parent.mkdir(parents=True, exist_ok=True)

                result = data.copy()
                result["chain"] = chain.id
                result["structure"] = str(pdb_file)
                result["name"] = "{}_{}".format(pdb_id.lower(), chain.id)

                if not pdb_file.exists() or self.overwrite:
                    self.logger.debug(
                        "Extracting chain %s from PDB %s",
                        chain.id, pdb_id)

                    # Store all captured log output in REMARK 999
                    general_logger = logging.getLogger("phyre_engine")
                    with phyre_engine.logutils.capture_log(general_logger) as log_buf:
                        # Select conformations.
                        for selector in self.conf_sel:
                            chain = selector.select(chain)
                        template = Template.build(pdb_id, result["chain"],
                                                  chain)

                        log_buf.seek(0)
                        template.remarks[999].extend(log_buf.readlines())

                    with pdb_file.open("w") as pdb_out:
                        template.write(pdb_out)
                else:
                    chain = pdb_parser.get_structure(
                        "", result["structure"])[0]["A"]
                    template = Template.build(pdb_id, result["chain"], chain)

                    self.logger.debug(
                        "Loaded existing chain %s of PDB %s from %s",
                        chain.id, pdb_id, pdb_file)
                result["sequence"] = template.canonical_seq
                result["original_residues"] = template.mapping
                result["canonical_indices"] = template.canonical_indices
                results.append(result)
        return results


class BuildTemplate(Component):
    """
    Build a :py:class:`phyre_engine.tools.template.Template` object from the
    keys written by :py:class:`.ChainBuilter`.

    A :py:class:`phyre_engine.tools.template.Template` object will be stored in
    the key ``template_obj`` in the pipeline state.  This object contains a
    :py:class:`Bio.PDB.Chain.Chain` object, as well as information about the
    template sequence and mapping of original residue numbers to new residue
    numbers.

    See :py:class:`.ChainPDBBuilder` for information about the fields that must
    be in the pipeline state.

    The fields ``original_residues`` and ``canonical_indices`` are removed from
    the pipeline state, as that information will be stored in the
    ``template_obj`` key.

    :param bool load_chain: If `False`, the PDB chain is not loaded from disk,
        and the `chain` attribute is set to `None`.
    """
    ADDS = ["template_obj"]
    REMOVES = ["original_residues", "canonical_indices"]
    REQUIRED = [
        "chain", "PDB", "structure", "original_residues",
        "canonical_indices", "sequence"
    ]

    def __init__(self, load_chain=True):
        self.load_chain = load_chain

    def run(self, data, config=None, pipeline=None):
        """Build Template object."""
        chain = None
        if self.load_chain:
            pdb_parser = Bio.PDB.PDBParser()
            chain = pdb_parser.get_structure("", data["structure"])[0]["A"]

        template = Template(
            data["PDB"], data["chain"], chain, data["original_residues"],
            data["sequence"], data["canonical_indices"])
        data["template_obj"] = template
        del data["original_residues"]
        del data["canonical_indices"]
        return data


class PDBSequence(Component):
    """
    Read a sequence from the ATOM records of a PDB structure.

    The pipeline state must have the ``template_obj`` key defined, pointing to
    a :py:class:`phyre_engine.tools.template.Template` object.

    This component adds the ``sequence`` key to the pipeline state. The value
    of the ``sequence`` field will be a Python string consisting of
    single-letter amino acids.
    """

    ADDS = ["sequence"]
    REMOVES = []
    REQUIRED = ["template_obj"]

    def run(self, data, config=None, pipeline=None):
        template = self.get_vals(data)
        data["sequence"] = template.canonical_seq
        return data

class AnnotateCATrace(Component):
    """
    Adds a ``ca_trace`` field indicating whether the
    :py:class:`phyre_engine.tools.template.Template` object ``template_obj`` is
    a C-α trace.

    This component will read the chain contained within the ``template_obj``,
    examining only those residues in the canonical sequence.
    """


    REQUIRED = ["template_obj"]
    ADDS = ["ca_trace"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Add a ``ca_trace`` field."""
        template = self.get_vals(data)
        # Read list of residues in canonical sequence
        residues = [template.chain[i] for i in template.canonical_indices]

        for res in residues:
            if len(res) != 1 or "CA" not in res:
                data["ca_trace"] = False
                return data
        data["ca_trace"] = True
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

    @property
    def REQUIRED(self):
        return [self.item_list]

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

        self.logger.info(
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
        self.logger.info(
            "Adding %d templates to the %d already present.",
            len(extra), len(data[self.item_list]))
        data[self.item_list].extend(extra)
        return data
