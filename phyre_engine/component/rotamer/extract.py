"""
Contains components for extracting information about the phi, psi and chi angles
of a residue.
"""
import math
import sys
from phyre_engine.component import Component
import Bio.PDB
from Bio.PDB import PDBParser
import phyre_engine.tools.rotamer.rotamer as rotamer
import pickle

# pylint: disable=abstract-method
class AngleExtractorBase(Component):
    """
    Extract φ, ψ and χ angles from all residues in a list of PDB files.

    This is a base class providing utility methods for extracting data. It does
    not supply a ``run()`` method directly.

    :ivar max_ca_distance: Maximum allowed Cα--Cα distance between consecutive
        residues. Residue pairs farther than this distance will be discarded.
        Default: 4.0Å
    :vartype max_ca_distance: float
    """

    max_ca_distance = 4.0

    def __init__(self, symmetric_chis=None):
        """
        :param set symmetric_chis: Set of amino acids for which the final
            rotamer should be considered symmetric. Symmetric rotamers are
            calculated modulo 180degrees. Usually, amino acids with terminating
            rings or ``O-C=O`` should be considered symmetric.

        .. seealso::

           Variable `SYMMETRIC_FINAL_CHI`
               In modules `tools.rotamer.data.dunbrack` and
               `tools.rotamer.data.molprobity`. Defines sets of amino acids with
               symmetric final χ angles according to the definitions of
               Dunbrack (1997) or MolProbity.
        """
        self.symmetric_chis = symmetric_chis if symmetric_chis else set()

    def _residue_triplets(self, residues):
        """
        Iterate over each triplet of residues.

        :return: Triplet of residues. If a residue is not present (due to a
            chain end), it is given as `None`.
        """
        prev_res = None
        current_res = next(residues, None)
        next_res = next(residues, None)

        while current_res is not None:
            yield (prev_res, current_res, next_res)
            prev_res = current_res
            current_res = next_res
            next_res = next(residues, None)

    def _require_atom(self, residue, atoms, feature, angle_index=None):
        """
        Raise an exception if a residue is missing a particular atom.

        Arguments are the same as for
        :py:class:`phyre_engine.tools.rotamer.MissingAtomError`, except that
        `atoms` must be a list (or other iterable) of atoms.

        :raise phyre_engine.tools.rotamer.MissingAtomError: Raised if `atom` is
            missing from `residue`.
        """
        for atom in atoms:
            if atom not in residue:
                raise rotamer.MissingAtomError(
                    residue,
                    feature,
                    atom,
                    angle_index)

    def _dihedrals(self, triplet):
        """
        Calculate φ and ψ angles (in degrees) for a given triplet.

        :raise phyre_engine.tools.rotamer.MissingAtomError: If a residue is
            missing a required atom.
        """

        # BioPython just raises a KeyError when we attempt to index an atom
        # that doesn't exist, which doesn't give us the residue on which the
        # lookup failed. Here, we add some tedious error checking to provide
        # more informative exceptions.
        self._require_atom(triplet[0], ("C",), rotamer.ResidueFeature.PHI)
        self._require_atom(triplet[1], ("N", "C", "CA"), rotamer.ResidueFeature.PHI)
        self._require_atom(triplet[2], ("N",), rotamer.ResidueFeature.PSI)

        phi = Bio.PDB.calc_dihedral(
            triplet[0]["C"].get_vector(),
            triplet[1]["N"].get_vector(),
            triplet[1]["CA"].get_vector(),
            triplet[1]["C"].get_vector())
        psi = Bio.PDB.calc_dihedral(
            triplet[1]["N"].get_vector(),
            triplet[1]["CA"].get_vector(),
            triplet[1]["C"].get_vector(),
            triplet[2]["N"].get_vector())
        return phi / math.pi * 180, psi / math.pi * 180


    def _triplet_angles(self, triplet):
        """
        Calculate φ, ψ, and χ angles for a triplet of residues.

        :param tuple triplet: Tuple of residues.
        :returns: `None` if the triplet is invalid, or a dictionary containing
            the keys ``residue``, ``torsion`` and ``sidechain`` as defined in
            :py:func:`.run`.

        :raise phyre_engine.tools.rotamer.MissingAtomError: If a residue is
            missing a required atom.
        :raise phyre_engine.tools.rotamer.UnknownResidueType: If a residue is of
            an unknown type. Currently, we only know about the standard 20 amino
            acids.
        """
        # Skip chain ends
        if any(r is None for r in triplet):
            return None

        self._require_atom(triplet[0], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)
        self._require_atom(triplet[1], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)
        self._require_atom(triplet[2], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)

        # Check distances between CA atoms
        if triplet[0]["CA"] - triplet[1]["CA"] > self.max_ca_distance:
            return None
        elif triplet[1]["CA"] - triplet[2]["CA"] > self.max_ca_distance:
            return None

        # Calculate angles
        dihedrals = self._dihedrals(triplet)
        sidechain = rotamer.Sidechain.calculate(
            triplet[1], self.symmetric_chis)
        return {
            "residue": triplet[1],
            "torsion": dihedrals,
            "sidechain": sidechain
        }

    def _residues(self, pdbs):

        parser = PDBParser()
        for pdb_file in pdbs:
            structure = parser.get_structure("input", pdb_file)
            for chain in structure[0]:
                chain_id = chain.get_id()
                residues = chain.get_residues()
                for triplet in self._residue_triplets(residues):
                    try:
                        residue_angles = self._triplet_angles(triplet)
                        if residue_angles is None:
                            continue
                        residue_angles["source"] = pdb_file
                        residue_angles["chain"] = chain_id
                        yield residue_angles, None
                    except rotamer.MissingAtomError as err:
                        yield None, err
                    except rotamer.UnknownResidueType as err:
                        yield None, err

class AngleExtractor(AngleExtractorBase):
    """
    Extract φ, ψ and χ angles from all residues in a list of PDB files.

    Adds a ``residues`` key to the pipeline state containing a list of
    dictionaries. See :py:meth`.run` for details on the fields of each residue
    dictionary.

    :ivar max_ca_distance: Maximum allowed Cα--Cα distance between consecutive
        residues. Residue pairs farther than this distance will be discarded.
        Default: 4.0Å
    :vartype max_ca_distance: float
    """
    REQUIRED = ["pdbs"]
    ADDS = ["residues"]
    REMOVES = ["pdbs"]

    def run(self, data, config=None, pipeline=None):
        """
        For each residue in each PDB given in the ``pdbs`` array of `data`,
        extract all residues and calculate the phi, psi and chi angles (that is,
        backbone dihedral angles and side-chain dihedral angles).

        Residues are considered in consecutive triplets: chain ends are
        discarded because phi and psi cannot both be calculated for these
        residues. If the CA--CA distance between two residues is greater than
        3.8A, the triplet is discarded.

        All chi angles are calculated for each residue.

        Residues missing required atoms are discarded with a warning.

        This method adds the ``residues`` key to `data`. This is a list of
        dictionaries, each containing the following keys:

        ``residue``
            The `Bio.PDB.Residue` object from which angles were extracted.

        ``source``
            File name of the PDB file from which the residue was extracted.

        ``chain``
            ID of the chain from which the residue was extracted.

        ``torsion``
            Tuple of (phi, psi) angles for this residue (in degrees).

        ``sidechain``
            A `phyre_engine.tools.rotamer.SideChain` object describing the chi
            angles of the residue.
        """
        pdbs = self.get_vals(data)
        residues = []
        for residue, exception in self._residues(pdbs):
            if exception:
                # FIXME: Need a better logging system
                print(exception, file=sys.stderr)
                continue

            residues.append(residue)
        del data["pdbs"]
        data["residues"] = residues
        return data

class AngleExtractorPickle(AngleExtractorBase):
    """
    Identical to :py:class`.AngleExtractor`, but residues are pickled for later
    reading.

    Each residue dictionary is appended to the pickle file, so can be read by
    unpickling the file multiple times. This class primarily exists so that
    we can extract all residues from multiple PDB files in a single pass and
    then process the residues in groups later to avoid exhausting the system
    memory.
    """
    REQUIRED = ["pdbs"]
    ADDS = []
    REMOVES = []

    class LightweightResidue():
        def __init__(self, res):
            self.res = res.get_resname()

        def get_resname(self):
            return self.res

    def __init__(self, pickle_file, *args):
        super().__init__(args)
        self.pickle_file = pickle_file

    def run(self, data, config=None, pipeline=None):
        with open(self.pickle_file, "wb") as pickle_fh:
            pdbs = self.get_vals(data)
            for residue, exception in self._residues(pdbs):
                if exception:
                    # FIXME: Need a better logging system
                    print(exception, file=sys.stderr)
                    continue

                residue["residue"] = AngleExtractorPickle.LightweightResidue(
                    residue["residue"])
                pickle.dump(residue, pickle_fh)
        return data


class AssignRotamers(Component):
    """
    For each residue in the ``residues`` list, assign a rotamer.

    Requires the ``sidechain`` attribute to be present on each element of the
    ``residues`` array.
    """
    REQUIRED = ["residues"]
    ADDS = []
    REMOVES = []

    def __init__(self, rotamers):
        """
        :param dict rotamers: Rotamer definitions.

        .. seealso::

            Module :py:data:`phyre_engine.tools.rotamer.data.dunbrack.__init__.ROTAMERS`
                Rotamer definitions according to Dunbrack's 1997 rotamer
                library.

            Module :py:data:`phyre_engine.tools.rotamer.data.molprobity.ROTAMERS`
                Rotamer definitions according to MolProbity.
        """
        self.rotamers = rotamers

    def run(self, data, config=None, pipeline=None):
        """
        Add a ``rotamer`` key to each residue.
        """
        for res in data["residues"]:
            res_rotamer = rotamer.Rotamer.find(res["sidechain"], self.rotamers)
            res["rotamer"] = res_rotamer
        return data
