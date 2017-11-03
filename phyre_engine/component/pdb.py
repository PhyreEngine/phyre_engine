"""
This module contains components for interacting with PDB (or mmCIF) files.
"""
import enum
import Bio.PDB
import phyre_engine.tools.pdb
import phyre_engine.tools.util
import phyre_engine.tools.template
from phyre_engine.component.component import Component
import tempfile
import os
import datetime
import urllib
import xml.etree.ElementTree

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
    REQUIRED = ["PDB"]
    ADDS = ["metadata"]
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

    def __init__(self, fields):
        self.fields = fields

    def run(self, data, config=None, pipeline=None):
        """Retrieve metadata from RCSB."""

        # Always retrieve the structureId so we can look up the correct
        # template. Not used right now, but it has no overhead and might be
        # useful in the future.
        columns = ["structureId"]

        # Include chain if available
        pdb_id = data["PDB"]
        if "chain" in data:
            pdb_id += "." + data["chain"]
            # Include chain length to force the web service to list all chains.
            # Otherwise, it can default to not listing *any* chains. We include
            # the chain ID to enable reverse lookup of templates for mapping
            # the response back onto the pipeline state.
            columns += ["chainId", "chainLength"]
        columns += list(self.fields.keys())

        if "metadata" not in data:
            data["metadata"] = {}

        # Retrieve data
        url = self.REST_URL.format(
            pdbids=pdb_id,
            columns=",".join(columns))
        with urllib.request.urlopen(url) as rest_data:
            root = xml.etree.ElementTree.parse(rest_data)

        records = root.findall("record")
        for record in records:
            for child in record:
                # Discard the "namespace" in the tag. It's not an actual XML
                # namespace: I think it indicates the report this field would
                # usually be a part of.
                tag = child.tag
                value = child.text
                if "." in tag:
                    tag = tag.split(".")[-1]

                if tag in self.fields:
                    if callable(self.fields[tag]):
                        value = self.fields[tag](value)
                    else:
                        converter = self.VALUE_TYPES[self.fields[tag]]
                        value = converter(value)
                    data["metadata"][tag] = value
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
