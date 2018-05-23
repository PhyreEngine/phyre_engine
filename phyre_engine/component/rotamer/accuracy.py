r"""
===============================
Calculating side-chain accuracy
===============================

Components for measuring the accuracy of side-chain construction.

There are several different definitions of side-chain accuracy. A particular χ
angle is usually regarded as "correct" if it is within 40° of the corresponding
angle in the native structure. These modules default to this definition. This
module provides classes for calculating both conditional and absolute side-chain
accuracy, as well as classes for calculating per-residue and all-atom RMSD (root
mean square deviation).

Conditional and absolute accuracy
---------------------------------

The *conditional* accuracy of a side-chain at a particular angle :math:`\chi_i`
is defined as the fraction of side-chains with a correct :math:`\chi_i` angle
given that the preceding angles :math:`\chi_{i-1}, \chi_{i-2}, \dots` are
correct. The *absolute* accuracy of an angle :math:`\chi_i` s defined as the
fraction of side-chains with all angles :math:`\chi_i, \chi_{i-1}, \dots`
correct.

Per-residue and all-atom RMSD
-----------------------------

For a set of residues of a specific type, there are two ways of calculating the
RMSD. Per-residue RMSD is calculated by finding the RMSD of each side-chain and
averaging the RMSDs over the number of residues. All-atom RSMD is instead
calculated by taking the sum of the squared deviations for each atom of each
residue, taking the square root and dividing by the number of atoms. Formally,
for :math:`N_{\text{res}}` residues, the side-chains of which have
:math:`N_{\text{sc}}` atoms, the per-residue RMSD is given by

.. math::

    \text{per-residue RMSD} =
    \frac{1}{N_\text{res}} \sum^{N_\text{res}}_i \left[
        \frac{1}{N_{\text{sc}}} \sqrt{
            \sum_j^{N_{\text{sc}}} (r_j - r'_j)^2
        } \, \right],

where :math:`r_j` and :math:`r'_j` are the positions of atom :math:`j` in the
model and native structures respectively.

The all-atom RMSD is calculated over the same set of residues:

.. math::

    \text{all-atom RMSD} = \frac{1}{N_\text{sc} N_\text{res}}
        \sqrt{ \sum_i^{N_\text{sc} N_\text{res}} (r_i - r'_i)^2 }.

The per-residue RMSD is used most commonly, but both definitions are given here
for the sake of clarity.

Required data
-------------

Each of the components in this module require a ``models`` array in the pipeline
state. Each element of this array should be a dictionary containing the fields
``native`` and ``model``. These can either be strings, in which case they are
interpreted as file paths, or :py:class:`~Bio.PDB.Structure.Structure` objects.
If a file path is supplied, it is assumed to be a PDB file. Each component will,
if necessary, create a ``sidechain_accuracy`` field in the model dictionary and
add a value to that dictionary. After a component is run, each element in the
models arrary will look like ``{..., "sidechain_accuracy": {...}}``.

"""
import Bio.PDB

from phyre_engine.component import Component
from phyre_engine.tools.rotamer.rotamer import Sidechain, MissingAtomError, \
    UnknownResidueType
from phyre_engine.tools.rotamer.data.generic import NUM_CHI_ANGLES
from math import sqrt
import sys



# pylint: disable=abstract-method
class SidechainMetric(Component):
    """
    Abstract base class implementing common methods for measuring side-chain
    accuracy.
    """
    REQUIRED = ["models"]
    ADDS = []
    REMOVES = []

    def _get_structure(self, model, field):
        """
        Get the structure of ``field`` in ``model``.

        Treats strings as file paths and parses them with PDBParser. Otherwise
        assumes that the field is of the correct type.

        e.g.
        >>> self._get_structure(model, "native")
        """
        if isinstance(model[field], Bio.PDB.Structure.Structure):
            return model[field]
        else:
            return Bio.PDB.PDBParser(QUIET=True).get_structure(
                model[field], model[field])

    def _models(self, models):
        """
        Used for iterating over each model in ````. Automatically adds the
        ``sidechain_accuracy`` field to each model if it doesn't already exist.
        The third parameter is the output dictionary.

        e.g.
        >>> for native, model, accuracy in self._models(models):
        >>>    # Some calculation...
        >>>    accuracy["conditional"] = 0.5
        """
        for model_dict in models:
            if "sidechain_accuracy" not in model_dict:
                model_dict["sidechain_accuracy"] = {}

            # Load files into structures if necessary
            native = self._get_structure(model_dict, "native")
            model = self._get_structure(model_dict, "model")
            yield native, model, model_dict["sidechain_accuracy"]

    def _residues(self, native, model):
        """
        Iterate over each pair of residues in the native and model. Raises an
        error if a residue pair is not equivalent.

        Yields the amino acid type, native side-chain and model side-chain:
        """
        native_residues = native.get_residues()
        model_residues = {r.id: r for r in model.get_residues()}

        for n_res in native_residues:
            if n_res.id not in model_residues:
                # FIXME: We need a proper way of handling warnings
                print(
                    "{} not present in model--skipping".format(n_res.id),
                    file=sys.stderr)
            else:
                aa = n_res.get_resname()
                yield aa, n_res, model_residues[n_res.id]

    def _sidechains(self, native, model):
        """
        Iterate over each pair of sidechains in the native and model.

        Yields the amino acid type, native side-chain and model side-chain:
        """
        for aa, n_res, m_res in self._residues(native, model):
            try:
                n_sc = Sidechain.calculate(n_res)
                m_sc = Sidechain.calculate(m_res)
                if n_sc is not None and m_sc is not None:
                    yield aa, n_sc, m_sc
            except MissingAtomError as err:
                # FIXME: We need a proper warning system
                print(err, file=sys.stderr)
            except UnknownResidueType as err:
                # FIXME: We need a proper warning system
                print(err, file=sys.stderr)

    def _chi_angles(self, native_sc, model_sc):
        """
        Iterate over each chi angle in the native and model residue.
        """
        yield from zip(native_sc.angles, model_sc.angles)

class SidechainAccuracy(SidechainMetric):
    """
    Base class for metrics calculating side-chain accuracies.
    """

    def __init__(self, correct_within=40):
        """
        :param float correct_within: χ angles within this number of degrees are
            judged to be "correct".
        """
        self.correct_within = correct_within

    def _count_dicts(self):
        """
        Generate dictionaries for counting the number of correct and incorrect
        chi angles.
        """
        # Indexed first by amino acid name and then chi indices
        num_correct = {}
        num_incorrect = {}
        for aa, n in NUM_CHI_ANGLES.items():
            if n > 0:
                num_correct[aa] = [0] * n
                num_incorrect[aa] = [0] * n
        return num_correct, num_incorrect

    def _accuracy(self, num_correct, num_incorrect):
        """Calculate accuracy from count of correct and incorrect chi angles."""
        accuracies = {}
        for aa in num_correct.keys():
            accuracies[aa] = []
            for good, bad in zip(num_correct[aa], num_incorrect[aa]):
                if good + bad == 0:
                    accuracies[aa].append(None)
                else:
                    accuracies[aa].append(good / (good + bad))
        return accuracies

class ConditionalAccuracy(SidechainAccuracy):
    """
    Calculate conditional accuracy for each χ angle of each residue.

    This component will add the ``conditional`` field to the
    ``sidechain_accuracy`` dictionary of each element of the ``models`` array.
    """

    def run(self, data, config=None, pipeline=None):
        models = self.get_vals(data)
        for native, model, result in self._models(models):
            num_correct, num_incorrect = self._count_dicts()

            for aa, native_sc, model_sc in self._sidechains(native, model):
                for i, (n_chi, m_chi) in enumerate(
                    self._chi_angles(native_sc, model_sc)):

                    ang_diff = abs(n_chi - m_chi) % 360
                    if ang_diff < self.correct_within:
                        num_correct[aa][i] += 1
                    else:
                        num_incorrect[aa][i] += 1
                        break

            accuracies = self._accuracy(num_correct, num_incorrect)
            result["conditional"] = accuracies
        return data

class AbsoluteAccuracy(SidechainAccuracy):
    """
    Calculate absolute accuracy for each χ angle of each residue.

    This component will add the ``absolute`` field to the ``sidechain_accuracy``
    dictionary of each element of the ``models`` array.
    """

    def run(self, data, config=None, pipeline=None):
        models = self.get_vals(data)
        for native, model, result in self._models(models):
            num_correct, num_incorrect = self._count_dicts()

            for aa, native_sc, model_sc in self._sidechains(native, model):
                # set to True when a chi angle is bad
                residue_is_bad = False
                for i, (n_chi, m_chi) in enumerate(
                    self._chi_angles(native_sc, model_sc)):

                    ang_diff = abs(n_chi - m_chi) % 360
                    if ang_diff < self.correct_within and not residue_is_bad:
                        num_correct[aa][i] += 1
                    else:
                        # Remaining chi angles are bad because this one was bad
                        residue_is_bad = True
                        num_incorrect[aa][i] += 1

            accuracies = self._accuracy(num_correct, num_incorrect)
            result["absolute"] = accuracies
        return data

class PerResidueRMSD(SidechainMetric):
    """
    Calculate average per-residue RMSD for each residue.
    """
    def _atoms(self, native_res, model_res):
        bb_atoms = set(("N", "C", "CA", "O"))
        native_atoms = sorted(
            [atom for atom in native_res if atom.element != "H"],
            key=lambda x: x.get_name())
        model_atoms = sorted(
            [atom for atom in model_res if atom.element != "H"],
            key=lambda x: x.get_name())

        if len(native_atoms) != len(model_atoms):
            # FIXME: We need a proper warning system
            print(AtomMismatchError(native_res, model_res), file=sys.stderr)

        for n_atm, m_atm in  zip(native_atoms, model_atoms):
            if n_atm.get_name() != m_atm.get_name():
                # FIXME: We need a proper warning system
                print(AtomMismatchError(native_res, model_res), file=sys.stderr)
            if n_atm.get_name() in bb_atoms:
                continue

            yield n_atm, m_atm


    def run(self, data, config=None, pipeline=None):
        models = self.get_vals(data)
        for native, model, result in self._models(models):
            rmsds = {aa: []
                     for aa, n in NUM_CHI_ANGLES.items()
                     if n > 0}

            for aa, native_res, model_res in self._residues(native, model):
                if aa not in rmsds:
                    continue

                square_distance = 0
                n = 0
                for n_atm, m_atm in self._atoms(native_res, model_res):
                    vec = n_atm.get_vector() - m_atm.get_vector()
                    square_distance += vec.normsq()
                    n += 1
                if n > 0:
                    rmsd = sqrt(square_distance) / n
                    rmsds[aa].append(rmsd)
                else:
                    # FIXME: Need a proper warning system
                    print("Error computing RMSD between {} and {} of structures {} and {}".format(
                        native_res, model_res, native, model), file=sys.stderr)

            average_rmsds = {}
            for aa, rmsd_list in rmsds.items():
                if len(rmsd_list) == 0:
                    average_rmsds[aa] = None
                    continue
                average_rmsds[aa] = sum(rmsd_list) / len(rmsd_list)

            result["rmsd"] = average_rmsds
        return data


class ResidueMismatchError(Exception):
    """
    Raised when two structures to be compared have mismatched residues.

    :ivar residues: Tuple of length 2 containing the mismatched residues.
    """

    ERR_MSG = "Mismatched residues '{}' and '{}'."

    def __init__(self, res1, res2):
        self.residues = (res1, res2)
        super().__init__(self.ERR_MSG.format(res1, res2))

class AtomMismatchError(Exception):
    """
    Raised when two side-chains have different sets of atoms.

    :ivar residues: Tuple of length2 containing the mismatched residues.
    """

    ERR_MSG = "Different sets of atoms in residues {} and {}: {} and {}."

    def __init__(self, res1, res2):
        self.residues = (res1, res2)
        super().__init__(self.ERR_MSG.format(
            res1, res2,
            [a.get_name() for a in res1],
            [a.get_name() for a in res2]
        ))
