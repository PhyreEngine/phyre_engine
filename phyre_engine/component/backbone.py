"""
This module contains tools for reconstructing the backbone of a protein from
a C-α trace.
"""
import shutil
import subprocess
import tempfile
from phyre_engine.component.component import Component
import phyre_engine.tools.external

class PD2CA2main(Component):
    """
    Use `PD2_ca2main
    <http://onlinelibrary.wiley.com/doi/10.1002/jcc.23330/abstract>`_ to
    reconstruct the backbone of a C-α trace.

    :param str database: Database location.
    :param list[str] flags: Flags to pass to ``pd2_ca2main``.
    :param dict[str] options: Options to pass to ``pd2_ca2main``.
    :param str bin_dir: Directory containing the ``pd2_ca2main`` executable if
        not in the system ``PATH``.

    .. seealso::

        :py:mod:`phyre_engine.tools.external`
            For detailed definitions of "flags" and "options".
    """
    REQUIRED = ["structure"]
    ADDS = []
    REMOVES = []

    PD2_CA2MAIN = phyre_engine.tools.external.ExternalTool({
        "input": "pose:io:pdb:i",
        "output": "pose:io:pdb:o",
    })

    def __init__(self, database, flags=None, options=None, bin_dir=None):
        self.database = database
        self.flags = flags if flags is not None else []
        self.options = options if options is not None else {}
        self.bin_dir = bin_dir

    def run(self, data, config=None, pipeline=None):
        """Run pd2_ca2main to construct backbone."""
        structure_path = self.get_vals(data)

        with tempfile.NamedTemporaryFile("rb") as tmp_pdb:
            options = self.options.copy()
            options["database"] = str(self.database)
            options["input"] = str(structure_path)
            options["output"] = tmp_pdb.name

            cmd_line = self.PD2_CA2MAIN(
                (self.bin_dir, "pd2_ca2main"),
                positional=None,
                flags=self.flags,
                options=options)
            subprocess.run(cmd_line, check=True)
            shutil.copy(tmp_pdb.name, str(structure_path))
        return data
