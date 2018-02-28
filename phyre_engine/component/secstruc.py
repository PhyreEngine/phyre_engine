"""
This module contains components for calculating the secondary structure of a
protein according to various methods.

Each of these components will add a list to the ``secondary_structure`` mapping
in the pipeline state, indexed by the name of the tool. The list that is added
must be the same length as the query sequence.

Each element of the list must be a dictionary, which must contain an
``assigned`` key giving the assigned secondary structure state. Where it is
available, tools should also include a ``confidence`` key, listing each possible
state and the confidence in that state.

For example, after running predictor ``foo`` the pipeline state might look like
this:

.. code-block:: python

    {
        # Various keys required to run the predictor,
        "secondary_structure": {
            "foo": [
                {"assigned": "C", "confidence": {"C": 0.6, "H": 0.3, "E": 0.1}},
                {"assigned": "C", "confidence": {"C": 0.7, "H": 0.2, "E": 0.1}},
                {"assigned": "C", "confidence": {"C": 0.8, "H": 0.1, "E": 0.1}},
                {"assigned": "C", "confidence": {"C": 0.6, "H": 0.3, "E": 0.1}},
                # ...
            ]
        }

"""

import enum
import subprocess
from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

SECONDARY_STRUCTURE_KEY = "secondary_structure"

class ThreeStateSS(enum.Enum):
    """Possible three-state secondary structure elements."""
    COIL = "C"
    HELIX = "H"
    STRAND = "E"

class EightStateSS(enum.Enum):
    """Possible eight-state secondary structure elements."""
    COIL = "C"
    ALPHA = "H"
    BETA_BRIDGE = "B"
    STRAND = "E"
    HELIX_3 = "G"
    HELIX_5 = "I"
    TURN = "T"
    BEND = "S"


class DSSP(Component):
    """
    Calculate secondary structure state using
    `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_.

    This component requires the ``structure`` field to be set, otherwise it has
    no source of tertiary structure from which to calculate the secondary
    structure.

    :param str bin_dir: Directory containing the ``mkdssp`` executable, if it is
        not in the system ``$PATH``.
    """
    CONFIG_SECTION = "dssp"

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
                res_id = line[5:10]
                aa_type = line[13]
                sec_struc = line[16]

                # Ignore missing residues.
                if aa_type == '!':
                    continue
                # Use "C" for coils.
                if sec_struc == ' ':
                    sec_struc = 'C'

                residue_ss = {
                    "assigned": sec_struc,
                    "confidence": {},
                    "res_id": int(res_id)}
                for state in EightStateSS:
                    confidence = 1.0 if state.value == sec_struc else 0.0
                    residue_ss["confidence"][state.value] = confidence
                dssp_mapping.append(residue_ss)
        return dssp_mapping
