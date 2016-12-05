"""
Module containing classes for selecting which conformation of a structure to
use.
"""

from abc import ABC, abstractmethod
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

class MutationSelector(ABC):
    """
    Interface defining residue-level conformation selectors.

    Structures submitted to the PDB may contain point mutations. For example,
    the first residue of structure ``3jqh`` appears as both a proline and a
    serine. The mutations are chosen according the "alternate location" field of
    a PDB file (the ``_atom_site.label_alt_id`` field of an mmCIF file).

    This abstract base class defines the interface of "mutation selector"
    classes, which are are responsible for selecting the sequence and structure
    of structures with point mutations.
    """


    @abstractmethod
    def select(self, chain):
        """
        Method called to select a single residue-level conformation.

        This method must return a tuple of ``Bio.PDB.Chain.Chain`` objects. It
        is considered acceptable to return *multiple* objects, corresponding to
        each sequence.
        """
        pass

    def clean_residue(self, residue):
        """
        Remove the "disordered" property of all atoms belonging to a residue.

        Once a concrete residue has been chosen from a DisorderedResidue object,
        the atoms within that residue will still be DisorderedAtom objects.
        However, because only one alternate location indicator field is
        available, those atoms cannot actually be disordered. This method will
        convert those atoms into Atom objects.

        Args:
            residue: A Bio.PDB.Residue.Residue object.

        Returns: A Bio.PDB.Residue.Residue object.
        """
        sanitised_res = Residue(
            residue.get_id(),
            residue.get_resname(),
            residue.get_segid())

        # Fix altloc ID of atoms. These will be DisorderedAtoms because
        # they have an altloc, but they can only contain one Atom
        # because the residue subsumes the altloc.
        for disordered_atom in residue:
            atom = disordered_atom.disordered_get_list()[0]
            atom.disordered_flag = 0
            atom.set_altloc(' ')
            sanitised_res.add(atom)
        return sanitised_res


class ArbitraryMutationSelector(MutationSelector):
    """Given a structure with point mutations, select an arbitrary sequence."""

    def select(self, chain):
        """Select the first conformation encountered."""
        new_chain = Chain(chain.get_id())

        for residue in chain:
            # is_disordered() == 2 for point mutations.
            # is_disordered() == 1 if residue contains disordered atoms
            if residue.is_disordered() == 2:
                first_res = residue.disordered_get_list()[0]
                new_chain.add(self.clean_residue(first_res))
            else:
                new_chain.add(residue)
        return (new_chain, )

class AllMutationSelector(MutationSelector):
    """Given a structure with point mutations, split it into separate
    sequences."""

    def select(self, chain):
        """
        Split a structure on point mutations.

        If a structure contains a point mutation, this selector will split it
        into two structures, one for each sequences.

        Raises:
            TooManyMutationsError: When more than two mutations are present.
        """

        # First, count the number of mutations
        mutation_conformations = set()
        for residue in chain:
            if residue.is_disordered() == 2:
                for conformation in residue.disordered_get_id_list():
                    mutation_conformations.add(conformation)
        if len(mutation_conformations) > 2:
            raise TooManyMutationsError(len(mutation_conformations))

        # Build a chain for each conformation:
        confs = {c: Chain(chain.get_id()) for c in mutation_conformations}
        for residue in chain:
            if residue.is_disordered() == 2:
                for conformation in residue.disordered_get_id_list():
                    conf_res = residue.disordered_get(conformation)
                    conf_res.disordered = 0
                    confs[conformation].add(self.clean_residue(conf_res))
            else:
                for chain in confs.values():
                    chain.add(residue)
        return tuple(confs.values())

class TooManyMutationsError(Exception):
    """Error thrown when a structure contains too many mutations."""

    def __init__(self, num_mutations):
        """Initialise a new error.

        Args:
            new_mutations: Number of mutations found.
        """
        msg = "Found too many ({}) mutations".format(num_mutations)
        super().__init__(msg)

class MicroConformationSelector(ABC):
    """
    Interface defining selectors for structures with microheterogeneity.

    Structures in the PDB can exhibit microheterogeneity; that is, individual
    atoms within each residue may be found at different locations. Heterogenous
    atoms are assigned different "alternate location" identifiers, and the
    "occupancy" field gives the probability of observing either conformation.

    This abstract base class defines the interface of "microconformation selector"
    classes, which are are responsible for selecting which conformations are
    chosen as templates.
    """


    @abstractmethod
    def select(self, chain):
        """
        Method called to select a single conformation.

        This method must return a single ``Bio.PDB.Chain.Chain`` object.
        """
        pass

class ArbitraryConformationSelector(MicroConformationSelector):
    """Select a single, arbitrary, conformation."""

    def select(self, chain):
        """Pick an arbitrary conformation."""

        sanitised_chain = Chain(chain.get_id())
        for res in chain:
            sanitised_res = Residue(
                res.get_id(),
                res.get_resname(),
                res.get_segid())

            for atom in res:
                if not atom.is_disordered():
                    # Keep atoms without alternate locations
                    sanitised_res.add(atom)
                else:
                    # Keep the first conformation
                    conformation = atom.disordered_get_list()[0]
                    # Ugly hack here. We need to flip the disordered flag, or
                    # PDBIO will complain when we try to write this atom.
                    conformation.disordered_flag = 0
                    conformation.set_altloc(' ')
                    sanitised_res.add(conformation)
            sanitised_chain.add(sanitised_res)
        return sanitised_chain

class PopulationConformationSelector(MicroConformationSelector):
    """Select the single most populated conformation.

    Each conformation is assigned a score according to the atoms it contains.
    The score is given by a weighted sum of "important" atoms. CA atoms are
    assigned a score of 4, CB atoms a score of 3, and N, C and O atoms a score
    of 2. All other atoms are assigned a score of 1.

    If the score assigned to each conformation is equal, ties are broken by
    choosing the conformation with the highest cumulative occupancy over
    backbone and CB atoms. Any remaining ties are broken by choosing the
    conformation with the lowest cumulative B factor over the same atoms.
    Any remaining ties are broken by picking an arbitrary conformation.
    """

    SCORES = {
        "CA": 4,
        "CB": 3,
        "N": 2,
        "C": 2,
        "O": 2,
    }

    def score_conformations(self, residue):
        """Assign scores to conformations, and return the ID of the
        top-ranked conformation."""

        # Build a list of all conformations.
        conformations = {}
        for atom in residue:
            if atom.is_disordered():
                weight = self.SCORES.get(atom.get_id(), 1)
                for conf in atom.disordered_get_id_list():
                    # Set initial value for each conformation
                    if conf not in conformations:
                        conformations[conf] = {
                            "conformation": conf,
                            "population": 0,
                            "occupancy": 0,
                            "b_factors": 0,
                        }
                    conf_atom = atom.disordered_get(conf)
                    conformations[conf]["population"] += weight
                    conformations[conf]["occupancy"] += conf_atom.get_occupancy()
                    conformations[conf]["b_factors"] += conf_atom.get_bfactor()

        # Take advantage of python's stable sorting to do a three-layer sort.
        # First, by b_factors, then occupancy then population. We want an
        # ascending sort on b_factors (lower == better), descending on occupancy
        # (higher == better) and descending on population.
        sorted_conf = sorted(
            conformations.values(),
            key=lambda c: c["b_factors"])
        sorted_conf = sorted(
            sorted_conf,
            key=lambda c: c["occupancy"],
            reverse=True)
        sorted_conf = sorted(
            sorted_conf,
            key=lambda c: c["population"],
            reverse=True)

        chosen = sorted_conf[0]["conformation"]
        return chosen

    def select(self, chain):
        """Pick the most-populated conformation."""

        sanitised_chain = Chain(chain.get_id())
        for res in chain:
            if res.is_disordered() == 1:
                chosen = self.score_conformations(res)
                sanitised_res = Residue(
                    res.get_id(),
                    res.get_resname(),
                    res.get_segid())

                for atom in res:
                    if not atom.is_disordered():
                        # Keep atoms without alternate locations
                        sanitised_res.add(atom)
                    else:
                        # Keep the first conformation
                        conformation = atom.disordered_get(chosen)
                        # Ugly hack here. We need to flip the disordered flag, or
                        # PDBIO will complain when we try to write this atom.
                        conformation.disordered_flag = 0
                        conformation.set_altloc(' ')
                        sanitised_res.add(conformation)
                sanitised_chain.add(sanitised_res)
            else:
                sanitised_chain.add(res)
        return sanitised_chain
