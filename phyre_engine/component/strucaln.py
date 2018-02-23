"""
This module contains components for aligning two protein structures.

See the documentation for each component for details. In general, each
component will align the structure pointed to by the ``model`` key of the
pipeline state onto the structured pointed to by the ``native`` key. You can
use the components in the :py:mod:`phyre_engine.component.alter` module to
move keys around the pipeline state.
"""
from phyre_engine.component.component import Component
import phyre_engine.tools.strucaln

class StructuralAlignment(Component):
    """
    Base class for components that do structural alignments.

    This class does not supply the ``ADDS`` property, because the scores added
    by each component may differ. Subclasses should explicitly define the
    scores that they add.
    """
    # pylint: disable=abstract-method
    REQUIRED = ["model", "native"]
    REMOVES = []

class TMAlign(StructuralAlignment):
    """
    Align two structures using `TMalign
    <https://zhanglab.ccmb.med.umich.edu/TM-align/>`_.

    This component will do a sequence-*independent* alignment between two
    structures. It will add the following scores to the pipeline state:

    ``TM``
        The TM-score ranges between 0 and 1, with a score above approximately
        0.5 indicating that two structures likely share the same fold. This
        will be stored as a 2-tuple of floats: the first score will indicate
        the TM-score normalised by the length of the ``model`` structure and
        the second by the length of the ``native``.

    ``rmsd``
        Root-mean-square deviation in ångström between the two structures.

    ``seqid``
        Sequence identity over the aligned region, between 0 and 1.

    ``structural_alignment``
        2-tuple containing the aligned sequences.

    ``superposition``
        A dictionary containing the file names of the superpositions generated
        by TMalign. This will not be generated unless the parameter
        `superposition` evaluates to a true value.
        The naming convention follows that of TMalign, with the prefix chosen
        by the `superposition` parameter. If `superposition` is ``sup``, the
        result will look like this:

        .. code-block:: python

            {"trace": "sup", "all": "sup_all", "atm": "sup_atm",
             "all_atm": "sup_all_atm", "all_atm_lig": "sup_all_atm_lig"}

    :param str superposition: Prefix to be used when generating superposition
        files. By default, superpositions are not generated.

    :param str bin_dir: Path to the ``TMalign`` executable. By default, it
        is assumed to be in your ``$PATH``.

    :param str executable: Alternate name for the executable.
    """
    ADDS = ["TM", "rmsd", "seqid", "structural_alignment"]

    def __init__(self, superposition=None, bin_dir=None, executable="TMalign"):
        self.superposition = superposition
        self.bin_dir = bin_dir
        self.executable = executable

    def run(self, data, config=None, pipeline=None):
        """Run TMalign to align two structures."""
        model, native = self.get_vals(data)
        aligner = phyre_engine.tools.strucaln.TMAlign(
            bin_dir=self.bin_dir,
            executable=self.executable)
        alignment = aligner.align(model, native, superpos=self.superposition)

        data["TM"] = alignment.tm
        data["rmsd"] = alignment.rmsd
        data["seqid"] = alignment.seqid
        data["structural_alignment"] = (
            alignment.sequences[0],
            alignment.sequences[2])
        if self.superposition:
            data["superposition"] = {
                "trace": alignment.superpositions[0],
                "all": alignment.superpositions[1],
                "atm": alignment.superpositions[2],
                "all_atm": alignment.superpositions[3],
                "all_atm_lig": alignment.superpositions[4],
            }
        return data

class TMScore(StructuralAlignment):
    """
    Align two structures with identical sequences using
    `TMscore <https://zhanglab.ccmb.med.umich.edu/TM-score/>`_.

    This is a sequence-*dependent* alignment: TMscore will fail if the two
    structures do not have the same sequence (or one is a subset of the other).
    This is useful, for example, for evaluating the quality of homology models.

    The following fields will be added to the pipeline state.

    ``TM``
        A 2-tuple of TM-scores, the same as for :py:class:`.TMAlign`.

    ``maxsub``
        MaxSub score, between 0 and 1.

    ``GDT_TS``
        Global Distance Test (Total Score), between 0 and 1.

    ``GDT_HA``
        Global Distance Test (High Accuracy), between 0 and 1.

    ``structural_alignment``
        Pair of aligned sequences.

    ``superposition``
        Similar to the ``superposition`` for :py:class:`.TMAlign`, but with
        only the keys ``trace`` and ``atm``.

    :param str superposition: Prefix to be used when generating superpositions.

    :param str bin_dir: Path to the TMscore executable.

    :param str executable: Name of the executable.
    """
    ADDS = ["TM", "maxsub", "GDT_TS", "GDT_HA", "sequence_alignment"]

    def __init__(self, superposition=None, bin_dir=None, executable="TMscore"):
        self.superposition = superposition
        self.bin_dir = bin_dir
        self.executable = executable

    def run(self, data, config=None, pipeline=None):
        """Run TMscore to calculate structural similarity."""
        model, native = self.get_vals(data)
        aligner = phyre_engine.tools.strucaln.TMScore(
            bin_dir=self.bin_dir,
            executable=self.executable)
        alignment = aligner.align(model, native, superpos=self.superposition)

        data["TM"] = alignment.tm
        data["maxsub"] = alignment.maxsub
        data["GDT_TS"] = alignment.gdt_ts
        data["GDT_HA"] = alignment.gdt_ha
        data["structural_alignment"] = alignment.sequences
        if self.superposition:
            data["superposition"] = {
                "trace": alignment.superpositions[0],
                "atm": alignment.superpositions[1],
            }
        return data
