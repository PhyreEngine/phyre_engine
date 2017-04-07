"""
This module contains tools to analyse side-chain rotamers.
"""
import math

import Bio.PDB

from . import data
from Bio.PDB.Residue import Residue

class Sidechain:
    """
    Representation of a side-chain, complete with angles and rotamer.

    :ivar angles: Tuple of angles (in degrees) corresponding to each χ angle.
    :ivar res_name: Name of the residue type.

    :vartype angles: tuple of floats
    :vartype res_name: str
    """

    def __init__(self, res_name, angles):
        """
        Initialise a new side-chain.

        :param str res_name: Residue type (three letters)
        :param tuple(float) angles: Tuple of angles
        :param `Bio.PDB.Residue` residue: Residue object to parse.
        """
        self.res_name = res_name
        self.angles = angles

    @classmethod
    def calculate(cls, residue):
        """
        Calculate side-chain angles from a residue.

        :param `Bio.PDB.Residue` residue: Residue to parse.
        :raise MissingAtomError: If an atom is missing from
        :return: A `Sidechain` object populated with the calculated angles or
            `None` if the specified residue doesn't have a side-chain (i.e.
            glycine and alanine).
        """
        if data.NUM_CHI_ANGLES[residue.get_resname()] == 0:
            return None

        angles = cls.calculate_angles(residue)
        return cls(residue.get_resname(), angles)

    @staticmethod
    def calculate_angles(residue):
        """
        Calculate χ angles for a given residue.

        :param `Bio.PDB.Residue` residue: Residue to parse.
        :return: Tuple of angles (in degree)
        :rtype: tuple(float)

        :raise MissingAtomError: When a required atom is not found.
        """
        num_angles = data.NUM_CHI_ANGLES[residue.get_resname()]
        angles = [None] * num_angles

        for i in range(0, num_angles):
            chi_atoms = []
            for atom in data.CHI_ATOMS[residue.get_resname()][i]:
                if atom not in residue:
                    raise MissingAtomError(residue, i, atom)
                else:
                    chi_atoms.append(residue[atom].get_vector())

            dihedral = math.degrees(Bio.PDB.calc_dihedral(*chi_atoms))
            if dihedral < 0:
                dihedral += 360
            angles[i] = dihedral
        return tuple(angles)

    def isclose(self, other, rel_tol=1e-09, abs_tol=0.0):
        """
        Check whether two rotamers are approximately the same.

        This will return `True` if the side-chains share the same residue name
        and the angles are approximately equal. Approximate equality is
        calculated using the `math.isclose` function.

        :param `Sidechain` other: Side-chain to compare.
        :param float rel_tol: Relative tolerance (see `math.isclose`.)
        :param float abs_tol: Absolute tolerance (see `math.isclose`.)
        :return: Boolean indicating whether side-chains are similar.
        .. seealso::
        Function `math.isclose`.
        """
        if self.res_name != other.res_name:
            return False
        for a, b in zip(self.angles, other.angles):
            if not math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol):
                return False
        return True

class Rotamer:
    """
    Class representing a rotamer.

    :ivar ranges: Tuple of ranges. A range is a two-element tuple containing the
        start and end angles (in degrees, from 0--360).
    :ivar name: Name of this tuple.
    """

    def __init__(self, res_name, rot_name):
        """
        Instantiate a new `Rotamer` object.

        :param str res_name: Name of the amino acid type of this rotamer.
        :param str rot_name: Rotamer name.
        """
        self.ranges = data.ROTAMERS[res_name][rot_name]
        self.name = rot_name

    @classmethod
    def find(cls, sidechain):
        """
        Returns the rotameric state of the given sidechain.

        :param `Sidechain` sidechain: The side-chain to search for.
        :return: Rotameric state of the sidechain or `None` for outliers.
        """

        # We'll just do a linear search for the moment, because I don't expect
        # this to be part of any high-performance pipelines. It would probably
        # be possible to reduce this to O(1) because we know the start and stop
        # of each range tend to be the same values.
        for rot in ALL_ROTAMERS[sidechain.res_name].values():
            if rot.valid_chi(sidechain.angles):
                return rot
        return None

    def valid_chi(self, angles):
        """
        Check whether each χ angle falls within the allowed range.

        :param angles: Tuple of angles to check.
        :return: `True` if all angles are valid or `False` otherwise.
        """
        for i, chi in enumerate(angles):
            if not (self.ranges[i][0] <= chi < self.ranges[i][1]):
                return False
        return True

class MissingAtomError(ValueError):
    """Thrown when a χ angle cannot be calculated due to missing angles.

    :ivar residue: The residue that is missing an atom.
    :ivar chi: The index of the χ angle (starting from 0).
    :ivar atom: The name of the missing atom.

    :vartype residue: `Bio.PDB.Residue`
    :vartype chi: int
    :vartype atom: str
    """

    _ERR_MSG = "Atom {atom} missing when calculating χ angle {chi} of " \
        "residue {residue}"

    def __init__(self, residue, chi, atom):
        """
        Initialise a new exception.

        :param `Bio.PDB.Residue` residue: The residue that is missing atoms.
        :param int chi: Index of the χ angle (starting from 0).
        :param str atom: Name of the missing atom.
        """
        super().__init__(
            self._ERR_MSG.format(atom=atom, residue=residue, chi=chi))
        self.residue = residue
        self.chi = chi
        self.atom = atom


# Pre-load list of rotamers
ALL_ROTAMERS = {}
for aa, rots in data.ROTAMERS.items():
    ALL_ROTAMERS[aa] = {}
    for rot_name in rots.keys():
        ALL_ROTAMERS[aa][rot_name] = Rotamer(aa, rot_name)
