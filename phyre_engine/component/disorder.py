"""
This module contains components for predicting the regions of a protein that are
intrinsically disordered.

Each tool in this module will a list to the ``disorder`` mapping in the pipeline
state, indexed by the name of the tool. The list that is added must be a list of
(possibly named) tuples of the same length as the protein sequence (in the
``sequence`` key). The first element should contain a boolean, indicating
whether the amino acid at that position is disordered (`True`) or not (`False`).
The remaining elements of each tuple are optional, and may include tool-specific
predictions such as confidence.

For example, if the disorder predictor ``foo`` is run, the pipeline state
afterwards might look like this:

.. code-block::

    {
        "sequence": "AG....", # Required for disorder predictors
        "disorder": {
            "foo": [
                (True, 0.5, ...), # Disordered, with some extra values
                (True, 0.1, ...),
                ...
            ]
        }
    }

"""
import collections
import json
import subprocess
import tempfile
from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

DISORDER_KEY = "disorder"

MobiDBResidue = collections.namedtuple("MobiDBResidue", "disordered prob")
class MobiDBLite(Component):
    """
    `MobiDB lite <http://protein.bio.unipd.it/mobidblite/>`_ is a meta-predictor
    of disorder. It combines nine fast predictors to quickly produce a consensus
    prediction.

    :param str bin_dir: Directory containing the ``mobidb-lite.py`` script.
    :param str supporting_bin_dir: Root directory of the binaries used by MobiDB
        lite (i.e. the directory passed via the ``--binDirectory`` option).
    """
    REQUIRED = ["sequence"]
    ADDS = [DISORDER_KEY]
    REMOVES = []


    TOOL_NAME = "mobidb-lite"
    EXECUTABLE_NAME = "mobidb-lite.py"

    MOBIDB_TOOL = ExternalTool()
    MOBDIB_DEFAULT_ARGS = {"threads": 1}
    MOBIDB_DEFAULT_FLAGS = {"longOutput"}

    def __init__(self, bin_dir, supporting_bin_dir):
        self.bin_dir = bin_dir
        self.supporting_bin_dir = supporting_bin_dir


    @staticmethod
    def parse_results(mobidb_output):
        mdb_results = json.loads(mobidb_output)
        disorder = []
        for state, prob in zip(mdb_results["consensus"], mdb_results["p"]):
            disorder.append(MobiDBResidue(state == "D", prob))
        return disorder

    def run(self, data, config=None, pipeline=None):
        """Run MobiDB lite to predict disorder."""
        sequence = self.get_vals(data)

        # Add disorder key if it's not present
        if DISORDER_KEY not in data:
            data[DISORDER_KEY] = {}

        # Write sequence to temp file and run mobidb-lite
        with tempfile.NamedTemporaryFile("w") as seq_file:
            print(">query", file=seq_file)
            print(sequence, file=seq_file)
            seq_file.flush()

            mobdib_options = self.MOBDIB_DEFAULT_ARGS.copy()
            mobdib_options["binDirectory"] = self.supporting_bin_dir
            command_line = self.MOBIDB_TOOL(
                (self.bin_dir, self.EXECUTABLE_NAME),
                flags=self.MOBIDB_DEFAULT_FLAGS,
                options=mobdib_options,
                positional=[seq_file.name])

            self.logger.info("Running '%s'", command_line)
            mobidb_proc = subprocess.run(
                command_line,
                stdout=subprocess.PIPE, check=True, universal_newlines=True)
            disorder = self.parse_results(mobidb_proc.stdout)
            data[DISORDER_KEY][self.TOOL_NAME] = disorder
        return data

