"""
This module contains components for interacting with PDB (or mmCIF) files.
"""
import enum
import Bio.PDB
import phyre_engine.tools.pdb
import phyre_engine.tools.util
from phyre_engine.component.component import Component

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

