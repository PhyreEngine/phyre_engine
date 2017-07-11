"""
Module containing classes for selecting which conformation of a structure to
use.
"""

from abc import ABC, abstractmethod
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import functools

class ConformationSelector(ABC):
    """
    Interface defining conformation selectors.

    Structures submitted to the PDB may contain point mutations. For example,
    the first residue of structure ``3jqh`` appears as both a proline and a
    serine. The mutations are chosen according the "alternate location" field of
    a PDB file (the ``_atom_site.label_alt_id`` field of an mmCIF file).

    Structures in the PDB can also exhibit microheterogeneity; that is,
    individual atoms within each residue may be found at different locations.
    Heterogenous atoms are assigned different "alternate location" identifiers,
    and the "occupancy" field gives the probability of observing either
    conformation.

    This abstract base class defines the interface of "conformation selector"
    classes, which are are responsible for selecting the sequence and structure
    of structures with point mutations and/or microheterogeneity. Subclasses
    must override the :py:meth:`.select` method.
    """


    @abstractmethod
    def select(self, chain):
        """
        Method called to select a single residue-level conformation.

        This method must return a single :py:class:`Bio.PDB.Chain.Chain` object.
        The returned chain may include multiple conformations, in which case it
        should be further processed.
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

        :param Bio.PDB.Residue.Residue residue: Residue object.
        :return: A cleaned residue object
        :rtype: :py:class:`Bio.PDB.Residue.Residue`
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

    def clean_atom(self, atom):
        """
        Convert a disordered atom to a non-disordered atom.

        This method sets the alternate location field of the atom to a space,
        and sets the ``disordered_flag`` property to zero.

        :param Bio.PDB.Atom.Atom atom: Atom to clean.
        :return: A cleaned atom object.
        :rtype: :py:class`Bio.PDB.Atom.Atom`
        """
        atom.disordered_flag = 0
        atom.set_altloc(' ')
        return atom

class PopulationConformationSelector(ConformationSelector):
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

    @functools.total_ordering
    class ConformationScore:
        """
        Describes the score of a given conformation, including weighted atomic
        population, cumulative occupancy and cumulative b-factor.

        Objects of this class are comparable and sortable. An object is greater
        than another if it has a higher overall score. Scores are compared first
        by population, then ties are broken by occupancy and then b-factor.
        """

        __slots__ = ["population", "occupancy", "b_factor"]

        def __init__(self, population=0, occupancy=0, b_factor=0):
            self.population = population
            self.occupancy = occupancy
            self.b_factor = b_factor

        def __add__(self, other):
            return self.__class__(
                self.population + other.population,
                self.occupancy + other.occupancy,
                self.b_factor + other.b_factor)

        def __gt__(self, other):
            if self.population > other.population:
                return True
            elif self.population == other.population:
                if self.occupancy > other.occupancy:
                    return True
                elif self.occupancy == other.occupancy:
                    if self.b_factor < other.b_factor:
                        return True
            return False

        def __eq__(self, other):
            return (self.population == other.population
                    and self.occupancy == other.occupancy
                    and self.b_factor == other.b_factor)

        def __repr__(self):
            return ("<{name}("
                    "{population}, {occupancy:.1f}, {b_factor:.1f})>").format(
                name=type(self).__name__,
                population=self.population,
                occupancy=self.occupancy,
                b_factor=self.b_factor)

class PopulationMicroHetSelector(PopulationConformationSelector):
    """
    For each residue, extract the most populated conformation, breaking ties
    on occupancy and then b-factor.

    .. seealso::

        :py:class:`phyre_engine.tools.conformation.PopulationConformationSelector`

            For a description of the selection algorithm.
    """

    def score_conformations(self, residue):
        """Assign scores to conformations, and return a tuple containing the ID
        and score of the top-ranked conformation.

        :param Bio.PDB.Residue residue: Residue from which conformations are
            scored.
        """

        # Build a list of all conformations.
        conformations = {}
        for atom in residue:
            if atom.is_disordered() and atom.get_id() in self.SCORES:
                weight = self.SCORES[atom.get_id()]
                for conf in atom.disordered_get_id_list():
                    # Set zero for this conformation
                    if conf not in conformations:
                        conformations[conf] = self.ConformationScore()

                    conf_atom = atom.disordered_get(conf)
                    conformations[conf] += self.ConformationScore(
                        weight,
                        conf_atom.get_occupancy(),
                        conf_atom.get_bfactor())

        sorted_confs = sorted(conformations.items(), key=lambda c: c[1])
        return sorted_confs[-1]

    def select(self, chain):
        """Pick the most-populated conformation."""

        sanitised_chain = Chain(chain.get_id())
        for res in chain:
            if res.is_disordered() == 1:
                chosen, score = self.score_conformations(res)

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
                        sanitised_res.add(self.clean_atom(conformation))
                sanitised_chain.add(sanitised_res)
            else:
                sanitised_chain.add(res)
        return sanitised_chain

class PopulationMutationSelector(PopulationConformationSelector):
    """
    For each point mutation in a PDB structure that is represented by an
    alternative location field, pick the mutation with the highest weighted
    population, breaking ties on occupancy and then b-factor.

    .. seealso::

        :py:class:`phyre_engine.tools.conformation.PopulationConformationSelector`

            For a description of the selection algorithm.
    """


    def score_conformations(self, residue):
        """
        Returns a sorted list of tuples containing conformation IDs and scores.
        """
        conformations = {}
        for conf in residue.disordered_get_id_list():
            if conf not in conformations:
                conformations[conf] = self.ConformationScore()

            conf_res = residue.disordered_get(conf)
            for atom in conf_res:
                if atom.get_id() in self.SCORES:
                    conformations[conf] += self.ConformationScore(
                        self.SCORES[atom.get_id()],
                        atom.get_occupancy(),
                        atom.get_bfactor())

        sorted_confs = sorted(conformations.items(), key=lambda c: c[1])
        return sorted_confs[-1]


    def select(self, chain):
        """Pick the most-populated conformation."""

        sanitised_chain = Chain(chain.get_id())
        for res in chain:
            if res.is_disordered() == 2:
                chosen, score = self.score_conformations(res)
                conf_res = res.disordered_get(chosen)
                sanitised_chain.add(self.clean_residue(conf_res))
            else:
                sanitised_chain.add(res)  #
        return sanitised_chain
