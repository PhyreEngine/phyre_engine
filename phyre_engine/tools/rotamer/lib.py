from phyre_engine.tools.rotamer.data.generic import NUM_CHI_ANGLES
class RotamerBin:
    """
    Represents a point in the rotamer library.

    :ivar probability: Probability of this rotamer occurring.
    :ivar mean_chi: Tuple of mean χ angles adopted by this rotamer.
    """
    probability = None
    mean_chi = None

    def __init__(self, num_chi):
        """Intialize an empty bin with ``num_chi`` angles."""
        self.mean_chi = [None] * num_chi

    def __repr__(self):
        return "<RotamerBin(probability={}, mean_chi={})>".format(
            self.probability,
            self.mean_chi)

class BackboneDependentRotamerLibrary:
    """
    Rotamer library giving rotamer probabilities and angles for different
    backbone conformations.

    This class represents a rotamer library. Users of this class may look up the
    probabilities and mean χ angles of each rotamer at each (ϕ, ψ) bin for each
    amino acid. Think of this class as containing a series of grids. A grid is
    defined for each amino acid type, and contains bins spanning (ϕ, ψ) space
    for that amino acid. Each bin contains several rotamers, which have a
    probability and a set of mean χ angles.

    Rotamers are first indexed by amino acid, then a (ϕ, ψ) tuple, then the
    rotamer name. The final element is a :py:class:`.RotamerBin`.
    """

    def __init__(self, rotamers, bin_width=10):
        """
        Create an empty rotamer library.

        :param dict rotamers: Dictionary of rotamer definitions. See, for
            example,
            :py:data:`phyre_engine.tools.rotamer.data.dunbrack.ROTAMERS`.
        :param int bin_width: Width of each bin.
        """
        self.lib = {}
        for aa, rotamer_dict in rotamers.items():
            self.lib[aa] = {}
            for phi in range(-180, 190, bin_width):
                for psi in range(-180, 190, bin_width):
                    self.lib[aa][(phi, psi)] = {}
                    for rotamer in rotamer_dict.keys():
                        rot_bin = RotamerBin(NUM_CHI_ANGLES[aa])
                        self.lib[aa][(phi, psi)][rotamer] = rot_bin

    def __getitem__(self, aa):
        """
        Index rotamer library by amino acid.

        >>> from phyre_engine.tools.rotamer.data.dunbrack import ROTAMERS
        >>> lib = BackboneDependentRotamerLibrary(ROTAMERS)
        >>> lib["SER"][(10, 10)][(1,)]
        """
        return self.lib[aa]

    def items(self):
        """
        Call to iterate over the rotamer library.

        Acts identically to the ``items()`` method of a dictionary. Keys will be
        amino acid names and values will be a dictionary indexed by a (ϕ, ψ)
        tuple.
        """
        return self.lib.items()
