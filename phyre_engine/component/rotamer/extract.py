"""
Contains components for extracting information about the phi, psi and chi angles
of a residue.
"""
import math
import sys
from phyre_engine.component import Component
import Bio.PDB
from Bio.PDB import PDBParser
import phyre_engine.tools.rotamer as rotamer
from phyre_engine.tools.rotamer import Rotamer

class AngleExtractor(Component):
    """
    Extract φ, ψ and χ angles from all residues in a list of PDB files.

    :ivar max_ca_distance: Maximum allowed Cα--Cα distance between consecutive
        residues. Residue pairs farther than this distance will be discarded.
        Default: 4.0Å
    :vartype max_ca_distance: float
    """

    REQUIRED = ["pdbs"]
    ADDS = ["residues"]
    REMOVES = ["pdbs"]

    max_ca_distance = 4.0;

    def __init__(self, symmetric_chis=set()):
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
        self.symmetric_chis = symmetric_chis

    def residue_triplets(self, residues):
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

    def require_atom(self, residue, atoms, feature, angle_index=None):
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

    def dihedrals(self, triplet):
        """
        Calculate φ and ψ angles (in degrees) for a given triplet.

        :raise phyre_engine.tools.rotamer.MissingAtomError: If a residue is
            missing a required atom.
        """

        # BioPython just raises a KeyError when we attempt to index an atom
        # that doesn't exist, which doesn't give us the residue on which the
        # lookup failed. Here, we add some tedious error checking to provide
        # more informative exceptions.
        self.require_atom(triplet[0], ("C",), rotamer.ResidueFeature.PHI)
        self.require_atom(triplet[1], ("N", "C", "CA"), rotamer.ResidueFeature.PHI)
        self.require_atom(triplet[2], ("N",), rotamer.ResidueFeature.PSI)

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


    def triplet_angles(self, triplet):
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

        self.require_atom(triplet[0], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)
        self.require_atom(triplet[1], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)
        self.require_atom(triplet[2], ("CA",), rotamer.ResidueFeature.CA_DISTANCE)

        # Check distances between CA atoms
        if triplet[0]["CA"] - triplet[1]["CA"] > self.max_ca_distance:
            return None
        elif triplet[1]["CA"] - triplet[2]["CA"] > self.max_ca_distance:
            return None

        # Calculate angles
        dihedrals = self.dihedrals(triplet)
        sidechain = rotamer.Sidechain.calculate(
            triplet[1], self.symmetric_chis)
        return {
            "residue": triplet[1],
            "torsion": dihedrals,
            "sidechain": sidechain
        }

    def run(self, data):
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
        # Output list
        res_list = []

        pdbs = self.get_vals(data)
        parser = PDBParser()
        for pdb_file in pdbs:
            structure = parser.get_structure("input", pdb_file)
            for chain in structure[0]:
                chain_id = chain.get_id()
                residues = chain.get_residues()
                for triplet in self.residue_triplets(residues):
                    try:
                        residue_angles = self.triplet_angles(triplet)
                        if residue_angles is None:
                            continue
                        residue_angles["source"] = pdb_file
                        residue_angles["chain"] = chain_id
                        res_list.append(residue_angles)
                    except rotamer.MissingAtomError as err:
                        print(err, file=sys.stderr)
                    except rotamer.UnknownResidueType as err:
                        print(err, file=sys.stderr)
        del data["pdbs"]
        data["residues"] = res_list
        return data

class AssignRotamers(Component):
    """
    For each residue in the ``residues`` list, assign a rotamer.
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

    def run(self, data):
        """
        Add a ``rotamer`` key to each residue.
        """
        for res in data["residues"]:
            rotamer = Rotamer.find(res["sidechain"], self.rotamers)
            res["rotamer"] = rotamer
        return data
