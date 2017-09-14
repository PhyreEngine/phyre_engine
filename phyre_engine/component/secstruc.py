"""
This module contains components for calculating the secondary structure of a
protein according to various methods.

Each of these components will add a list to the ``secondary_structure`` mapping
in the pipeline state, indexed by the name of the tool. The list that is added
must contain a list of (possibly named) tuples of the same length as the query
sequence. The exact contents of each tuple is component-specific, but the first
element must contain the secondary structure state.

For example, after running predictor ``foo`` the pipeline state might look like
this:

.. code-block::

    {
        # Various keys required to run the predictor,
        "secondary_structure": {
            ``foo``: [
                ("C", 0.4, ...), # Coil, some value of 0.4, any extra values
                ("C", 0.1, ...),
                ("E", 0.2, ...),
                # ...
            ]
        }

"""

from collections import namedtuple
import subprocess
import sys
from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

SECONDARY_STRUCTURE_KEY = "secondary_structure"

DSSPResidue = namedtuple("DSSPResidue", "ss")
class DSSP(Component):
    """
    Calculate secondary structure state using
    `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_.

    This component requires the ``structure`` field to be set, otherwise it has
    no source of tertiary structure from which to calculate the secondary
    structure. If the ``sequence`` field is set, then an exception will be
    raised if the output from DSSP does not contain exactly the same number of
    residues as the sequence.

    :param str bin_dir: Directory containing the ``mkdssp`` executable, if it is
        not in the system ``$PATH``.
    """
    REQUIRED = ["structure"]
    ADDS = [SECONDARY_STRUCTURE_KEY]
    REMOVES = []

    TOOL_NAME = "dssp"
    MKDSSP = ExternalTool()

    def __init__(self, bin_dir=None):
        self.bin_dir = bin_dir

    def run(self, data, config=None, pipeline=None):
        """Calculate ``secondary_structure`` key."""
        structure = self.get_vals(data)

        # Create secondary_structure key if it is not present
        if SECONDARY_STRUCTURE_KEY not in data:
            data[SECONDARY_STRUCTURE_KEY] = {}

        # Run DSSP on the structure file and read the output
        mkdssp_cmd_line = self.MKDSSP(
            (self.bin_dir, "mkdssp"),
            options={"input": structure})
        dssp_proc = subprocess.run(
            mkdssp_cmd_line, universal_newlines=True,
            check=True, stdout=subprocess.PIPE)
        dssp_mapping = self.parse_dssp(dssp_proc.stdout.split("\n"))

        if "sequence" in data:
            if len(data["sequence"]) != len(dssp_mapping):
                raise LengthMismatchError(dssp_mapping, data["sequence"])

        data[SECONDARY_STRUCTURE_KEY][self.TOOL_NAME] = dssp_mapping
        return data

    @staticmethod
    def parse_dssp(dssp_lines):
        """
        Parse lines of output from DSSP.

        :return: List of tuples containing the residue ID and secondary
            structure state.
        :rtype: list[tuple(int, str)]
        """
        residue_section = False
        dssp_mapping = []
        for line in dssp_lines:
            if line.startswith("  #  RESIDUE AA"):
                residue_section = True
            elif residue_section and len(line) > 17:
                aa_type = line[13]
                sec_struc = line[16]

                # Ignore missing residues.
                if aa_type == '!':
                    continue
                # Use "C" for coils.
                if sec_struc == ' ':
                    sec_struc = 'C'

                dssp_mapping.append(DSSPResidue(sec_struc))
        return dssp_mapping

class LengthMismatchError(ValueError):
    """
    Raised when the number of residues in a secondary structure assignment does
    not match the length of the query sequence.

    :param list ss: Secondary structure assignment.
    :param str sequence: Sequence for which secondary structure was assigned.
    """
    ERR_MSG = (
        "Secondary structure assignment length did not match sequence length "
        "({} and {} residues).")

    def __init__(self, ss, sequence):
        super().__init__(self.ERR_MSG.format(len(ss), len(sequence)))
