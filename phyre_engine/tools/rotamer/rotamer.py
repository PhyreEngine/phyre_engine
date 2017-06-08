"""
This module contains tools to analyse side-chain rotamers.
"""
import math
import enum

import Bio.PDB

from .data.generic import AMINO_ACIDS, NUM_CHI_ANGLES, CHI_ATOMS

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
    def calculate(cls, residue, final_chi_ranges=None):
        """
        Calculate side-chain angles from a residue.

        :param `Bio.PDB.Residue` residue: Residue to parse.
        :param dict final_chi_ranges: Definitions of the angular ranges allowed
            for the final χ atom.
        :raise MissingAtomError: If an atom is missing from
        :return: A `Sidechain` object populated with the calculated angles or
            `None` if the specified residue doesn't have a side-chain (i.e.
            glycine and alanine).

        .. seealso::

        Section :ref:`description-of-rotamer-variables`
            For a description of the :py:data:`ROTAMERS` and :py:data:`FINAL_CHI_RANGE`
            variables.
        """
        if residue.get_resname() not in AMINO_ACIDS:
            raise UnknownResidueType(residue.get_resname())
        if NUM_CHI_ANGLES[residue.get_resname()] == 0:
            return None

        angles = cls.calculate_angles(residue, final_chi_ranges)
        return cls(residue.get_resname(), angles)

    @staticmethod
    def calculate_angles(residue, final_chi_ranges=None):
        """
        Calculate χ angles for a given residue.

        Parameters are defined the same as for :py:meth:`calculate`.

        :return: Tuple of angles (in degree)
        :rtype: tuple(float)

        :raise MissingAtomError: When a required atom is not found.
        :raise UnknownResidueType: When an unknown residue type is passed in.
        """
        if residue.get_resname() not in AMINO_ACIDS:
            raise UnknownResidueType(residue.get_resname())

        num_angles = NUM_CHI_ANGLES[residue.get_resname()]
        angles = [None] * num_angles

        for i in range(0, num_angles):
            chi_atoms = []
            for atom in CHI_ATOMS[residue.get_resname()][i]:
                if atom not in residue:
                    raise MissingAtomError(residue, atom, ResidueFeature.CHI, i)
                else:
                    chi_atoms.append(residue[atom].get_vector())

            dihedral = math.degrees(Bio.PDB.calc_dihedral(
                chi_atoms[0], chi_atoms[1],
                chi_atoms[2], chi_atoms[3]))

            if dihedral < 0:
                dihedral += 360
            angles[i] = dihedral

        # Correct final chi angle of residues with symmetric final units.
        if (final_chi_ranges is not None
            and residue.get_resname() in final_chi_ranges):

            # If the angle is not in the allowed range, flip it by 180 degrees.
            if angles[-1] not in final_chi_ranges[residue.get_resname()]:
                angles[-1] = (angles[-1] + 180) % 360
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

    def __init__(self, res_name, rot_name, ranges):
        """
        Instantiate a new `Rotamer` object.

        :param str res_name: Name of the amino acid type of this rotamer.
        :param str rot_name: Rotamer name.
        """
        self.res_name = res_name
        self.ranges = ranges
        self.name = rot_name

    @classmethod
    def find(cls, sidechain, rotamers):
        """
        Returns the rotameric state of the given sidechain.

        :param `Sidechain` sidechain: The side-chain to search for.
        :return: Rotameric state of the sidechain or `None` for outliers.
        """

        # We'll just do a linear search for the moment, because I don't expect
        # this to be part of any high-performance pipelines. It would probably
        # be possible to reduce this to O(1) because we know the start and stop
        # of each range tend to be the same values.
        for rot_name, rot_angles in rotamers[sidechain.res_name].items():
            rot = Rotamer(sidechain.res_name, rot_name, rot_angles)
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
            if chi not in self.ranges[i]:
                return False
        return True

class ResidueFeature(enum.Enum):
    """
    Different types of angles that we may want to measure.
    """
    # : Phi (φ) backbone angle
    PHI = "φ"
    # : Psi (ψ) backbone angle
    PSI = "ψ"
    # : Chi (χ) side-chain angles
    CHI = "χ"
    # : Cα-Cα distance
    CA_DISTANCE = "Cα-Cα "

class MissingAtomError(ValueError):
    """Thrown when a χ angle cannot be calculated due to missing angles.

    :ivar residue: The residue that is missing an atom.
    :ivar feature: The residue feature being calculated (e.g. φ, ψ, χ angles).
    :ivar angle_index: The index of the χ angle (starting from 0), or `None` if
        not relevant.
    :ivar atom: The name of the missing atom.

    :vartype residue: `Bio.PDB.Residue`
    :vartype feature: :py:class:`.ResidueFeature`
    :vartype angle_index: `None` or int
    :vartype atom: str
    """

    _ERR_MSG_NOINDEX = "Atom {atom} missing when calculating {feature} of " \
        "residue {residue}"
    _ERR_MSG_INDEX = "Atom {atom} missing when calculating " \
        "{angle_index} of {feature} for residue {residue} "

    def __init__(self, residue, atom, feature, angle_index=None):
        """
        Initialise a new exception.

        :param `Bio.PDB.Residue` residue: The residue that is missing atoms.
        :param int chi: Index of the χ angle (starting from 0).
        :param str atom: Name of the missing atom.
        """
        if angle_index is None:
            err_template = self._ERR_MSG_NOINDEX
        else:
            err_template = self._ERR_MSG_INDEX
        err_msg = err_template.format(
            residue=residue,
            feature=feature,
            angle_index=angle_index,
            atom=atom
        )
        super().__init__(err_msg)
        self.residue = residue
        self.feature = feature
        self.angle_index = angle_index
        self.atom = atom

class UnknownResidueType(Exception):
    """
    Thrown when a residue is of an unknown type.

    :ivar type: Residue type that caused the failure.
    """
    _ERR_MSG = "Can't calculate side-chain angles for residue type {}"

    def __init__(self, res_type):
        super().__init__(self._ERR_MSG.format(res_type))
        self.type = res_type
