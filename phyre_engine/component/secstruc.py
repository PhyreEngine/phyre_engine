"""
This module contains components for calculating the secondary structure of a
protein according to various methods.

Each of these components will add a ``secondary_structure`` key to the pipeline
state, which will contain a list of named tuples. The exact contents of each
tuple is component-specific, but the first two elements must be named
``residue_id`` and ``sec_struc``. These correspond to the ID of the residue with
secondary structure state ``sec_struc``.
"""

from collections import namedtuple
import subprocess
import sys
from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

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
    REQUIRED = ["structure"]
    ADDS = ["secondary_structure"]
    REMOVES = []

    MKDSSP = ExternalTool()
    ResidueSecStruc = namedtuple("ResidueSecStruc", "residue_id sec_struc")

    def __init__(self, bin_dir=None):
        self.bin_dir = bin_dir

    def run(self, data, config=None, pipeline=None):
        """Calculate ``secondary_structure`` key."""
        structure = self.get_vals(data)

        # Run DSSP on the structure file and read the output
        mkdssp_cmd_line = self.MKDSSP(
            (self.bin_dir, "mkdssp"),
            options={"input": structure})
        dssp_proc = subprocess.run(
            mkdssp_cmd_line,
            check=True, stdout=subprocess.PIPE)
        dssp_lines = dssp_proc.stdout.decode(sys.stdout.encoding).split("\n")
        dssp_mapping = self.parse_dssp(dssp_lines)
        data["secondary_structure"] = dssp_mapping
        return data

    @classmethod
    def parse_dssp(cls, dssp_lines):
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

                # Bail if we didn't read a residue ID.
                if not res_id.strip():
                    break
                # Ignore missing residues.
                if aa_type == '!':
                    continue
                # Use "C" for coils.
                if sec_struc == ' ':
                    sec_struc = 'C'

                dssp_mapping.append(cls.ResidueSecStruc(int(res_id), sec_struc))
        return dssp_mapping
