"""
This module contains components for interacting with PDB (or mmCIF) files.
"""
import enum
import Bio.PDB
import Bio.PDB.MMCIF2Dict
import phyre_engine.tools.pdb
import phyre_engine.tools.util
import phyre_engine.tools.template
from phyre_engine.component.component import Component
from phyre_engine.tools.jmespath import JMESExtensions
import jmespath
import tempfile
import os
import datetime
import urllib
import xml.etree.ElementTree
import io

class StructureType(enum.Enum):
    """Possible structure types."""
    PDB = "PDB"
    MMCIF = "MMCIF"

class UnknownStructureTypeError(RuntimeError):
    """Raised when an unknown structure type is encountered."""
    ERR_MSG = "Unknown structure type. Available types: {}".format(
        ", ".join([s.value for s in StructureType]))

    def __init__(self):
        super().__init__(self.ERR_MSG)

class ReadStructure(Component):
    """
    Read a structure (in PDB or mmCIF format) into a
    :py:class:`Bio.PDB.Structure.Structure` object. This component will read
    the structure pointed to by the ``structure`` key of the pipeline state
    and will add the ``structure_obj`` key containing the parsed structure. The
    ``structure_type`` key will also be added, indicating the type of file that
    was parsed.

    The ``structure`` key may be a file name, file handle or
    :py:class:`pathlib.Path` object.

    ..warning ::

        If a stream is passed in the ``structure`` key, it *must* be
        possible to seek back to the beginning of the file handle.
    """
    REQUIRED = ["structure"]
    ADDS = ["structure_obj", "structure_type"]
    REMOVES = []

    PARSERS = {
        StructureType.PDB: Bio.PDB.PDBParser(QUIET=True),
        StructureType.MMCIF: Bio.PDB.MMCIFParser(QUIET=True)
    }

    @staticmethod
    def _parse_structure(parser, struc_in):
        try:
            structure = parser.get_structure("", struc_in)
            if list(structure.get_models()):
                return structure
        except Exception:
            pass
        return None

    def run(self, data, config=None, pipeline=None):
        """Read a structure file."""
        structure_file = self.get_vals(data)
        with phyre_engine.tools.util.Stream(structure_file, "r") as struc_in:
            for struc_type in StructureType:
                parser = self.PARSERS[struc_type]
                structure = self._parse_structure(parser, struc_in)
                if structure is not None:
                    data["structure_obj"] = structure
                    data["structure_type"] = struc_type.value
                    return data
                struc_in.seek(0)
        raise UnknownStructureTypeError()

class PDBSeq(Component):
    """
    Read a sequence from the atom records of a protein structure into the
    ``sequence`` key of the pipeline state. This component reads the
    ``structure_obj`` key, which must be a
    :py:class:`Bio.PDB.Structure.Structure` object.
    """
    REQUIRED = ["structure_obj"]
    ADDS = ["sequence"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Reading sequence from structure."""
        structure = self.get_vals(data)
        sequence, _ = phyre_engine.tools.pdb.atom_seq(structure.get_residues())
        data["sequence"] = sequence
        return data

class TooManyChainsError(RuntimeError):
    """Raised when a monomer contains too many chains."""

    ERR_MSG = "Expected one chain, saw {}"
    def __init__(self, num_chains):
        super().__init__(self.ERR_MSG.format(num_chains))

class ConvertToMonomer(Component):
    """
    Extract a single chain from a structure. If the structure contains multiple
    chains, raise :py:exc:`phyre_engine.component.pdb.TooManyChainsError`.

    This component will modify the ``structure_obj`` key in place, converting it
    to a :py:class:`Bio.PDB.Chain.Chain` object.
    """
    REQUIRED = ["structure_obj"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Extract chain from a structure."""
        structure = self.get_vals(data)
        chains = list(structure.get_chains())
        if len(chains) > 1:
            raise TooManyChainsError(len(chains))
        data["structure_obj"] = chains[0]
        return data

class ConvertToTemplate(Component):
    """
    Convert a single chain in the ``structure_obj`` key of the pipeline state to
    a template using :py:meth:`phyre_engine.tools.template.Template.build`. This
    will write the santised structure to a new PDB file, and replace the
    ``structure`` and ``structure_obj`` keys with the template.

    This component will place the template file in the current working
    directory. By default, a unique filename is used. This can be overridden by
    setting the `file_name` parameter in order to avoid the proliferation of
    ugly file names when the pipeline is run multiple times.

    :param str file_name: Use this fixed file name instead of a unique name.
    """

    ADDS = ["structure"]
    REQUIRED = ["structure_obj"]
    REMOVES = []

    def __init__(self, file_name=None):
        self.file_name = file_name


    def _open_structure(self):
        if self.file_name is None:
            file_des, file_name = tempfile.mkstemp(
                suffix="-template.pdb", dir=os.getcwd(), text=True)
            file_handle = os.fdopen(file_des)
        else:
            file_name = self.file_name
            file_handle = open(self.file_name, "w")
        return file_handle, file_name

    def run(self, data, config=None, pipeline=None):
        """Convert ``structure_obj`` to a template."""
        structure_obj = self.get_vals(data)

        file_handle = None
        try:
            file_handle, file_name = self._open_structure()

            template = phyre_engine.tools.template.Template.build(structure_obj)
            template.write(file_handle)
            data["structure"] = file_name
            data["structure_obj"] = template.chain
        finally:
            file_handle.close()
        return data

class TemplateMapping(Component):
    """
    Read the mapping between a template and the original residue IDs.
    """

    ADDS = ["residue_mapping"]
    REQUIRED = ["structure"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Read residue mapping from template."""
        structure_file = self.get_vals(data)
        template = phyre_engine.tools.template.Template.load(structure_file)
        data["residue_mapping"] = template.mapping
        return data

class MMCIFMetadata(Component):
    """
    Parse metadata from an mmCIF file.

    This component will look up the mmCIF file specified by the ``PDB`` field
    in the pipeline state from the directory `mmcif_dir`. It will then parse
    the mmCIF file into a dictionary using
    :py:class:`~Bio.PDB.MMCIF2Dict.MMCIF2Dict`.

    Each pair of expressions in the `fields` dictionary is then evaluated: the
    value of each pair in the dictionary is a JMESPath expression evaluated
    against the metadata dictionary, and the key is the field name to be placed
    into the ``metadata`` dictionary in the pipeline state. For example:

    >>> from phyre_engine.component.pdb import MMCIFMetadata
    >>> meta = MMCIFMetadata(
    ...     mcif_dir="/path/to/mmcif/files",
    ...     fields={
    ...         "resolution": "to_number(_reflns.d_resolution_high)",
    ...         "title": "_struct.title",
    ...     })
    >>> results = meta.run({"PDB": "12AS"})
    >>> results["metadata"]
    {"resolution": 2.2, "title": "ASPARAGINE SYNTHETASE MUTANT C51A..."}

    .. note::

        Some mmCIF files can be very large, so by default mmCIF files are
        filtered to remove all lines beginning with ``ATOM`` or ``HETATM``.
        This is technically incorrect, because the mmCIF specification does not
        require the ``_atom_site`` records to begin with the ``group_PDB``
        field. However, all the mmCIF files from the PDB obey this format for
        (supposed) compatibility with some line-based PDB-format parsers.
        If this pre-filter is causing problems, it can be disabled by setting
        `prefilter=False`.

    .. note::

        Most of the keys in the mmCIF dictionary will contain a dot (``.``),
        which is parsed by JMESPath as a sub-expression unless it is escaped.
        To escape a JMESPath expression, surround it in double quotes. This
        can be slightly tricky depending on how you are passing strings into
        this class. In Python and YAML, it is probably easiest to specify the
        string with single quotes, and then use double quotes to specify a
        literal JMESPath expression: ``'"a.b"'``.

    :param str mmcif_dir: Base directory of mmCIF files, divided according
        to the middle two letters of the PDB identifier.

    :param dict[str, str] fields: Mapping of field names to JMESPath
        expressions evaluated against the parsed mmCIF dictionary.

    :param bool prefilter: Strip lines beginning with ``ATOM`` or ``HETATM``
        before parsing the mmCIF file. This will greatly speed up the parsing
        of large files, but technically violates the mmCIF specification. The
        default is `True`, because the prefilter should work with all the
        mmCIF files from the Protein Data Bank.
    """

    REQUIRED = ["PDB"]
    ADDS = ["metadata"]
    REMOVES = []

    def __init__(self, mmcif_dir, fields, prefilter=True):
        self.mmcif_dir = mmcif_dir
        self.fields = fields
        self.prefilter = prefilter

    @staticmethod
    def _prefilter(mmcif_in):
        """
        Read file handle mmcif_in, strip lines beginning with ATOM or HETATM
        and return an io.StringIO buffer containing the pre-filtered lines.
        """
        buf = io.StringIO()
        for line in mmcif_in:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                buf.write(line)
        buf.seek(0)
        return buf

    def run(self, data, config=None, pipeline=None):
        """Parse metadata from mmCIF file."""
        pdb_id = self.get_vals(data)

        mmcif_file = phyre_engine.tools.pdb.find_pdb(
            pdb_id, suffix_list=(".cif", ".cif.gz"),
            base_dir=self.mmcif_dir)
        data.setdefault("metadata", {})


        with phyre_engine.tools.pdb.open_pdb(mmcif_file) as mmcif_in:
            if self.prefilter:
                mmcif_in = self._prefilter(mmcif_in)

            mmcif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(mmcif_in)
            jmes_extensions = JMESExtensions(mmcif_dict)
            jmes_opts = jmespath.Options(custom_functions=jmes_extensions)
            for field, jmespath_expr in self.fields.items():
                value = jmespath.search(jmespath_expr, mmcif_dict, jmes_opts)
                data["metadata"][field] = value

        return data


class RCSBMetadata(Component):
    """
    Look up metadata about a PDB entry using the RCSB's web services.

    This module uses the RCSB's `RESTful web services
    <https://www.rcsb.org/pdb/software/rest.do>`_ to look up metadata about the
    PDB ID given by the ``PDB`` field of the pipeline state. If the field
    ``chain`` is present metadata for that particular chain is retrieved. For a
    list of fields, see https://www.rcsb.org/pdb/results/reportField.do.
    Metadata will be stored in the ``metadata`` dictionary.

    The RCSB's web services do not provide any metadata about type, so you will
    need to explicitly supply the type of the fields you require. This is done
    by specifying fields and types as a mapping. The following types are
    available:

    int
        Integer number.
    float
        Real number.
    str
        Arbitrary text.
    date
        A date.

    Custom types may be passed by supplying a callable rather than the name of
    a type to the constructor.

    .. warning::

        Parsing a date will cause a :py:class:`date.date` object to be stored
        in the pipeline state, which will cause errors when it is serialised.
        If you only need the date for display purposes, it is recommended that
        you parse it as text (``str``). If you do parse it as a ``date``,
        remove it from the pipeline state before serialisation.

    >>> from phyre_engine.component.pdb import RCSBMetadata
    >>> meta = RCSBMetadata({
    ...     "resolution": "float",
    ...     "releaseDate": "date",
    ... })
    >>> meta.run({"PDB": "11as"})
    >>> {"PDB": "11as", "metadata": {
    ...     "resolution": 2.5,
    ...     "releaseDate": date(1998, 12, 30),
    ... }

    :param list[str] field: List of fields to look up.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    #: Mapping of type names to conversion functions.
    VALUE_TYPES = {
        "int": int,
        "float": float,
        "str": str,
        "date": lambda date_str: datetime.datetime.strptime(
            date_str, "%Y-%m-%d").date()
    }

    #: Format string used to build the REST URL.
    REST_URL = (
        "https://www.rcsb.org/pdb/rest/customReport.xml?pdbids={pdbids}"
        "&customReportColumns={columns}"
        "&service=wsfile"
        "&format=xml"
    )

    def __init__(self, fields, chain=False):
        self.fields = fields
        self.chain = chain


    def _record_dict(self, record):
        """
        Convert an XML record returned by the web service into a python dict.
        This method does no type conversion. It returns a string for each
        value.
        """
        record_dict = {}
        for child in record:
            # Discard the "namespace" in the tag. It's not an actual XML
            # namespace: I think it indicates the report this field would
            # usually be a part of.
            tag = child.tag
            value = child.text
            if "." in tag:
                tag = tag.split(".")[-1]
            record_dict[tag] = value
        return record_dict


    def _record_key(self, record_dict):
        """
        Returns the identifier used to index templates from a record dictionary
        generated from _record_dict.
        """
        if self.chain:
            template_key = (
                record_dict["structureId"].lower(),
                record_dict["chainId"])
        else:
            template_key = record_dict["structureId"].lower()
        return template_key


    def _record_metadata(self, record_dict):
        """
        Convert an untyped metadata dictionary generated by _record_dict into
        a typed metadata dictionary.
        """
        metadata = {}
        for tag in self.fields:
            value = record_dict[tag]
            if value == "null":
                metadata[tag] = None
            else:
                if callable(self.fields[tag]):
                    value = self.fields[tag](value)
                else:
                    converter = self.VALUE_TYPES[self.fields[tag]]
                    value = converter(value)
                metadata[tag] = value
        return metadata


    def _identifiers(self, templates):
        """
        Returns the list of structure IDs and columns to retrieve. Also returns
        a mapping of template identifiers to templates so that we can look up
        each template from the returned results later.
        """
        # Always retrieve the structureId so we can look up the correct
        # template.
        columns = ["structureId"]
        if self.chain:
            # Include chain length to force the web service to list all chains.
            # Otherwise, it can default to not listing *any* chains. We include
            # the chain ID to enable reverse lookup of templates for mapping
            # the response back onto the pipeline state.
            columns += ["chainId", "chainLength"]

        columns += list(self.fields.keys())
        # Build a map of templates indexed by PDB ID and optionally chain. This
        # is necessary because we might have duplicates in the list of
        # templates, and the web service will not return duplicates so we can't
        # just zip the list of templates and web service data together.
        # We will also build the list of IDs to retrieve at the same time.
        pdb_ids = set()
        template_map = {}
        for template in templates:
            # PDB IDs are ASCII-only, so we don't need to worry about unicode
            # matching (yet…)
            pdb_id = template["PDB"].lower()
            if self.chain:
                pdb_id += "." + template["chain"]
                template_key = template["PDB"].lower(), template["chain"]
            else:
                template_key = template["PDB"].lower()
            if template_key not in template_map:
                template_map[template_key] = []
            template_map[template_key].append(template)
            pdb_ids.add(pdb_id)
        return list(pdb_ids), columns, template_map

    def run(self, data, config=None, pipeline=None):
        """Retrieve metadata from RCSB."""
        templates = self.get_vals(data)
        pdb_ids, columns, template_map = self._identifiers(templates)

        # Work in groups of 512 IDs at a time to avoid getting a 400 bad
        # request response.
        records = []
        while pdb_ids:
            pdb_ids_chunk = pdb_ids[-512:]
            pdb_ids[-512:] = []
            # Retrieve data
            url = self.REST_URL.format(
                pdbids=",".join(pdb_ids_chunk),
                columns=",".join(columns))
            self.logger.debug("Retrieving %s", url)
            with urllib.request.urlopen(url) as rest_data:
                root = xml.etree.ElementTree.parse(rest_data)
            records.extend(list(root.findall("record")))

        for record in records:
            record_dict = self._record_dict(record)
            template_key = self._record_key(record_dict)
            metadata = self._record_metadata(record_dict)

            for template in template_map[template_key]:
                if "metadata" not in template:
                    template["metadata"] = {}
                template["metadata"].update(metadata)
        return data



class FastResolutionLookup(Component):
    """
    Quickly look up the resolution of all templates from the source mmCIF
    files.

    This component loops over the ``templates`` list, looking up the best
    resolution from the source mmCIF file. Each template must have the ``PDB``
    field defined. The field ``resolution`` is added to each template. If
    the template was not resolved by X-ray crystallography, the resoltion will
    be set to `None`.

    .. note::

        The parser used here should work without any problems, but is fairly
        primitive. It simply loops over the mmCIF files line by line, looking
        for a line starting with ``_reflns.d_resolution_high``. This is fast,
        but could fail with mmCIF files that are deliberately obfuscated.

    :param str mmcif_dir: Root directory containing mmCIF files.
    """
    REQUIRED = ["templates"]
    REMOVES = []
    ADDS = []

    def __init__(self, mmcif_dir):
        self.mmcif_dir = mmcif_dir

    def run(self, data, config=None, pipeline=None):
        """Reading resolution information from mmCIF files."""
        templates = self.get_vals(data)
        cache = {}

        for template in templates:
            if template["PDB"] not in cache:
                mmcif_file = phyre_engine.tools.pdb.find_pdb(
                    template["PDB"], suffix_list=(".cif", ".cif.gz"),
                    base_dir=self.mmcif_dir)
                resolution = None
                with phyre_engine.tools.pdb.open_pdb(mmcif_file) as mmcif_in:
                    for line in mmcif_in:
                        if line.startswith("_reflns.d_resolution_high"):
                            resolution = float(line.split()[1])
                            break
                cache[template["PDB"]] = resolution

            template["resolution"] = cache[template["PDB"]]
        return data


class FindStructure(Component):
    """
    Using the ``PDB`` and ``chain`` fields, find the corresponding structure
    file on disk and add it to the pipeline state. By default, the
    ``structure`` field is set; this can be changed by setting the `field`
    parameter.


    :param str template_dir: Root directory of the template library.
    """

    REQUIRED = ["PDB", "chain"]
    ADDS = []
    REMOVES = []

    def __init__(self, template_dir, field="structure"):
        self.template_dir = template_dir
        self.field = field

    def run(self, data, config=None, pipeline=None):
        """Find template for a given PDB and chain."""
        pdb_id, chain_id = self.get_vals(data)
        pdb_path = phyre_engine.tools.pdb.find_pdb(
            pdb_id, chain_id, self.template_dir)
        if pdb_path is None:
            raise FileNotFoundError(
                "No template for PDB {} (chain {}) in {}".format(
                    pdb_id, chain_id, self.template_dir))
        data[self.field] = str(pdb_path)
        return data
