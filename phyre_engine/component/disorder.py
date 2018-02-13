"""
This module contains components for predicting the regions of a protein that are
intrinsically disordered.

Each tool in this module will a list to the ``disorder`` mapping in the pipeline
state, indexed by the name of the tool. The list that is added must be the same
length as the protein sequence (in the ``sequence`` key). Each element of the
list should be a dictionary, which must include at least the ``assigned`` key,
indicating the assigned state.

Each tool should also, if possible, add a ``confidence`` key enumerating each
possible state and the confidence of that state.

For example, if the disorder predictor ``foo`` is run, the pipeline state
afterwards might look like this:

.. code-block::

    {
        "sequence": "AG...", # Required for disorder predictors
        "disorder": {
            "foo": [
                {"assigned": "D", "confidence": {"O": 0.1, "D": 0.9}},
                {"assigned": "D", "confidence": {"O": 0.2, "D": 0.8}},
                # ...
            ]
        }
    }

"""
import enum
import json
import os
import pathlib
import subprocess
import tempfile
from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

DISORDER_KEY = "disorder"

class DisorderStates(enum.Enum):
    """Possible disorder states."""
    STRUCTURED = "S"
    DISORDERED = "D"

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
        """Parse long MobiDB lite result string."""
        mdb_results = json.loads(mobidb_output)
        disorder = []
        for state, prob in zip(mdb_results["consensus"], mdb_results["p"]):
            disorder.append({
                "assigned": state,
                "confidence": {
                    DisorderStates.DISORDERED.value: prob,
                    DisorderStates.STRUCTURED.value: 1 - prob
                }
            })
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

            # Mobidb will rudely return no output when no disordered regions are
            # found. In those cases, we just don't add the mobidb-lite key to
            # the disordered predictor.
            if mobidb_proc.stdout.strip():
                disorder = self.parse_results(mobidb_proc.stdout)
                data[DISORDER_KEY][self.TOOL_NAME] = disorder
        return data


class Disopred(Component):
    """
    Run `DISOPRED <http://bioinf.cs.ucl.ac.uk/psipred/?disopred=1>`_ to predict
    disordered regions of protein structure.

    This component takes a shortcut compared to the stock ``run_disopred.pl``
    script supplied with DISOPRED: it will use an existing ``mtx`` file
    generated, for example, by :py:class:`phyre_engine.component.hhsuite.PSSM`.
    This component also does not predict binding binding sites using
    ProtBind.

    :param str data_dir: Directory containing DISOPRED data files.
    :param str dso_lib_dir: Directory containing DISOPRED library files.
    :param str bin_dir: Directory containing the DISOPRED executables.
    """

    REQUIRED = ["pssm"]
    ADDS = ["disorder"]
    REMOVES = []

    #: Adjustable DISOPRED2 false positive rate, from 1-10.
    DISOPRED2_FPR = 5

    DISOPRED2 = ExternalTool()
    DISO_NEU_NET = ExternalTool()
    DISO_NEIGHB = ExternalTool()
    COMBINE = ExternalTool()

    def __init__(self, data_dir, dso_lib_dir, bin_dir=None):
        self.data_dir = data_dir
        self.dso_lib_dir = dso_lib_dir
        self.bin_dir = bin_dir

    @staticmethod
    def parse_results(diso_in):
        """
        Parse disopred output file into the format described in
        :py:mod:`.disorder`.

        The disopred format looks like this:

        .. code-block:: none

            #         ----- DISOPRED version 3.1 -----
            # Disordered residues are marked with asterisks (*)
            #    Ordered residues are marked with dots (.)
                1 M * 0.78
                2 K * 0.62
                3 T . 0.45
                4 A . 0.37
                5 Y . 0.20

        We parse this into the following list:

        .. code-block:: none

            [
                {"assigned": "D", "confidence": {"S": 0.22, "D": 0.78}},
                {"assigned": "D", "confidence": {"S": 0.38, "D": 0.62}},
                {"assigned": "S", "confidence": {"S": 0.55, "D": 0.45}},
                {"assigned": "S", "confidence": {"S": 0.63, "D": 0.37}},
                {"assigned": "S", "confidence": {"S": 0.80, "D": 0.20}},
            ]

        :param file diso_in: File handle pointing to DISOPRED output.
        """
        disorder = []
        for line in diso_in:
            line = line.strip()
            if line.startswith("#"):
                continue
            _index, _aa, state, score = line.split()

            if state == "*":
                state = DisorderStates.DISORDERED
            else:
                state = DisorderStates.STRUCTURED
            score = float(score)

            disorder.append({
                "assigned": state.value,
                "confidence": {
                    DisorderStates.DISORDERED.value: score,
                    DisorderStates.STRUCTURED.value: 1 -score
                }
            })
        return disorder

    def run(self, data, config=None, pipeline=None):
        """Run DISOPRED to predict disorder."""
        pssms = self.get_vals(data)
        mtx_file = pssms["mtx"]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)

            # $args = join ' ', "$EXE_DIR/disopred2", join('/', $out_dir, $base), $mtx_fn, "$DATA_DIR/", $DISO2_FPR, "\n";
            # system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk\n";
            disopred2_results = tmpdir / "disopred2"
            self.logger.info("Predicting disorder with DISOPRED2.")
            disopred2_cmd = self.DISOPRED2(
                (self.bin_dir, "disopred2"),
                positional=(
                    disopred2_results,
                    mtx_file,
                    self.data_dir + "/",
                    self.DISOPRED2_FPR))
            self.logger.debug("Running %s", disopred2_cmd)
            subprocess.run(disopred2_cmd, check=True)



            # $args = join ' ', "$EXE_DIR/diso_neu_net", "$DATA_DIR/weights.dat.nmr_nonpdb", $mtx_fn, ">", $nndiso_fn, "\n";
            # system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk";
            diso_neu_net_results = tmpdir / "diso_neu_net"
            with diso_neu_net_results.open("wb") as neu_net_out:
                self.logger.info("Running neural network classifier.")
                diso_neu_net_cmd = self.DISO_NEU_NET(
                    (self.bin_dir, "diso_neu_net"),
                    positional=(
                        pathlib.Path(self.data_dir, "weights.dat.nmr_nonpdb"),
                        mtx_file))
                self.logger.debug("Running %s > %s",
                                  diso_neu_net_cmd, diso_neu_net_results)
                subprocess.run(diso_neu_net_cmd, stdout=neu_net_out,
                               check=True)

            # $args = join ' ', "$EXE_DIR/diso_neighb", $mtx_fn, "$DATA_DIR/dso.lst", ">", $dnb_fn, "\n";
            # system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk\n";
            diso_neighb_results = tmpdir / "diso_neighb"
            with diso_neighb_results.open("wb") as neighb_out:
                self.logger.info("Running nearest neighbour classifier.")
                diso_neighb_cmd = self.DISO_NEIGHB(
                    (self.bin_dir, "diso_neighb"),
                    positional=(
                        mtx_file,
                        pathlib.Path(self.data_dir, "dso.lst")))
                self.logger.debug("Running %s > %s",
                                  diso_neighb_cmd, diso_neighb_results)
                environment = dict(os.environ)
                environment["DSO_LIB_PATH"] = str(self.dso_lib_dir) + "/"
                subprocess.run(diso_neighb_cmd, stdout=neighb_out, check=True,
                               env=environment)


            # $args = join ' ', "$EXE_DIR/combine", "$DATA_DIR/weights_comb.dat", $diso2_fn, $nndiso_fn, $dnb_fn, ">", $diso3_fn, "\n";
            # system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk";
            disopred_results = pathlib.Path("disorder.diso")
            with disopred_results.open("wb") as diso_out:
                self.logger.info("Combining disordered residue predictions.")
                combine_cmd = self.COMBINE(
                    (self.bin_dir, "combine"),
                    positional=(
                        pathlib.Path(self.data_dir, "weights_comb.dat"),
                        # Add ".diso" suffix to disopred2 output
                        str(disopred2_results) + ".diso",
                        diso_neu_net_results,
                        diso_neighb_results))
                self.logger.debug("Running %s > %s",
                                  combine_cmd, disopred_results)
                subprocess.run(combine_cmd, stdout=diso_out, check=True)

        with disopred_results.open("r") as diso_in:
            disorder = self.parse_results(diso_in)
        if "disorder" not in data:
            data["disorder"] = {}
        data["disorder"]["disopred"] = disorder
        return data
