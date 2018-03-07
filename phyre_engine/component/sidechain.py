"""Components for reconstructing side-chains."""

import re
from pathlib import Path
import subprocess

from phyre_engine.component import Component
from phyre_engine.tools.external import ExternalTool

class Scwrl4(Component):
    """
    Run `SCWRL4 <http://dunbrack.fccc.edu/scwrl4/>`_ to reconstruct
    the side-chains of a model.

    :param str bin_dir: Directory containing the ``scwrl4`` executable.
    :param bool overwrite: If `True`, always run ``scwrl4``; otherwise,
        an existing output file will be used as-is.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["model"]

    SCWRL4 = ExternalTool(flag_map={
        "input": "i",
        "output": "o",
        "sequence": "s",
        "parameters": "p",
        "frame": "f",
        "graph": "g",
        "workspace": "w",
        "symmetry": "%",
        "crystal": "#",
        "disable_subrotamers": "v",
        "omit_hydrogens": "h",
        "disable_terminal_capping": "t",
    })

    CONFIG_SECTION = "scwrl4"

    def __init__(self, bin_dir=None, overwrite=False):
        self.bin_dir = bin_dir
        self.overwrite = overwrite

    def run(self, data, config=None, pipeline=None):
        """Run SCWRL4 to reconstruct side-chains."""
        model = self.get_vals(data)
        outfile = Path(model).with_suffix(".scwrl4.pdb")

        if not outfile.exists() or self.overwrite:
            command_line = self.SCWRL4(
                (self.bin_dir, "Scwrl4"),
                options={"input": model, "output": outfile})
            program = subprocess.run(
                command_line, check=True, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Scwrl doesn't always set a sensible exit value, so if we see
            # "^Err$" on standard output, treat it as an error.
            if re.search("^Err$", program.stdout, re.MULTILINE):
                self.logger.error(
                    ("Error running SCWRL4. Command line: %s\n"
                    "Standard output: %s\n"
                    "Standard error: %s\n"),
                    command_line, program.stdout, program.stderr)
                raise subprocess.CalledProcessError(
                    program.returncode, command_line,
                program.stdout, program.stderr)
        data["model"] = str(outfile)
        return data
