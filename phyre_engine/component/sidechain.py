"""Components for reconstructing side-chains."""

from pathlib import Path
import subprocess

from phyre_engine.component import Component
from phyre_engine.tools.external import ExternalTool

class Scwrl4(Component):
    """
    Run `SCWRL4 <http://dunbrack.fccc.edu/scwrl4/>`_ to reconstruct
    the side-chains of a model.

    :param str bin_dir: Directory containing the ``scwrl4`` executable.
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


    def __init__(self, bin_dir=None):
        self.bin_dir = bin_dir

    def run(self, data, config=None, pipeline=None):
        """Run SCWRL4 to reconstruct side-chains."""
        model = self.get_vals(data)
        outfile = Path(model).with_suffix(".scwrl4.pdb")

        command_line = self.SCWRL4(
            (self.bin_dir, "scwrl4"),
            options={"input": model, "output": outfile})
        subprocess.run(command_line, check=True)
        data["model"] = str(outfile)
        return data