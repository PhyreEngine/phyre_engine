"""
Tools for predicting residue-residue contacts from sequence information.

These tools all create a ``contacts`` list in the pipeline state. This list
contains a list of tuples (optionally, named tuples created with
:py:class:`phyre_engine.tools.util.NamedTuple`). The first two elements of each
tuple must be the residue IDs of the contacting residues, and the next element
must be the confidence of that contact. The remaining fields are unspecified.
For example:

.. code-block:: python

    {"contacts": [
        <Contact i=1, j=10, confidence=0.99>,
        <Contact i=2, j=9, confidence=0.95>,
        # ...
    ]}
"""
import pathlib
import subprocess

import phyre_engine.component as component
import phyre_engine.tools.util as util
import phyre_engine.tools.external as external

class MetaPsicov(component.Component):
    """
    Use `MetaPsicov <http://bioinf.cs.ucl.ac.uk/MetaPSICOV/>`_ to predict
    contacts.

    This component is fairly stupid. All it does is call the modified version
    of the ``run_metapsicov`` script packaged with the conda version of
    MetaPsicov. A more robust version of this component would replicate the
    logic of ``run_metapsicov`` and call the individual components. Pull
    requests will be accepted.

    :param str uniref90: Path to UniRef90 database. This must have been
        preprocessed with ``formatdb`` into a (legacy) BLAST database.

    :param str uniref100: Path to UniRef100 database. This should be
        preprocessed with ``esl-sfetch`` for use with HMMER.

    :param str hhblitsdb: Path to an hhsuite database such as ``uniboost``
        for generating deep MSAs.

    :param str pdb70: Path to the PDB70 hhsuite database, used for domain
        parsing.

    :param str bin_dir: Directory containing the ``run_metapsicov`` script
        and its supporting binaries.

    :param int cpu: Number of CPUs to use for components of MetaPsicov that
        support parallelisation. Some elements are parallelised using MPI, so
        you may need to set up a hostfile to use large numbers of processors.

    :param str work_dir: Directory in which to store intermediate files.

    :param bool overwrite: If `True`, remove the working directory and run
        MetaPsicov from scratch. Otherwise, if the ``stage3`` output file
        is present, it is parsed directly and MetaPsicov is not run. If the
        ``stage3`` output file does not exist, then MetaPsicov will resume
        where it left off.
    """
    REQUIRED = ["sequence"]
    ADDS = ["contacts"]
    REMOVES = []

    CONFIG_SECTION = "metapsicov"

    class MetaPsicovContact(util.NamedTuple):
        """
        Named tuple containing the following fields:

        ``i``, ``j``:
            Residue indices of the predicted contact.

        ``confidence``
            Confidence as assigned by MetaPsicov.
        """
        FIELDS = "i j confidence"

    METAPSICOV = external.ExternalTool()

    def __init__(self, uniref90, uniref100, hhblitsdb, pdb70,
                 bin_dir=None, cpu=1, work_dir="metapsicov", overwrite=False):
        self.uniref90 = uniref90
        self.uniref100 = uniref100
        self.hhblitsdb = hhblitsdb
        self.pdb70 = pdb70
        self.bin_dir = bin_dir
        self.cpu = cpu
        self.work_dir = work_dir
        self.overwrite = overwrite

    @classmethod
    def parse_results(cls, results_file):
        """Parse MetaPsicov state3 result."""
        # Records look like this:
        # i j 0 8 confidence

        results = []
        with util.Stream(results_file, "r") as results_in:
            for line in results_in:
                records = line.strip().split()
                contact = cls.MetaPsicovContact(
                    int(records[0]), int(records[1]),
                    float(records[4]))
                results.append(contact)
        return results

    def _run_metapsicov(self, work_dir, query_seq):
        """Run metapsicov in the given working dir."""
        # Write query sequence to a FASTA file
        query_seq_file = work_dir / "query.fasta"
        with query_seq_file.open("w") as query_out:
            query_out.write(">Query\n{}\n".format(query_seq))

        cmd_line = self.METAPSICOV(
            (self.bin_dir, "run_metapsicov"),
            positional=[
                query_seq_file,
                self.uniref90,
                self.uniref100,
                self.hhblitsdb,
                self.pdb70,
            ],
            flags=["keep-tempfiles"],
            options={
                "job-id": "metapsicov",
                "cpu": self.cpu,
                "work-dir": work_dir,
            })
        self.logger.debug("Running %s", cmd_line)
        subprocess.run(cmd_line, check=True)

    def run(self, data, config=None, pipeline=None):
        """Run MetaPsicov to predict contacts."""
        query_seq = self.get_vals(data)

        work_dir = pathlib.Path(self.work_dir)
        output_file = work_dir / "metapsicov.metapsicov.stage3"

        if self.overwrite or not output_file.exists():
            work_dir.mkdir(exist_ok=True)
            self._run_metapsicov(work_dir, query_seq)

        contacts = self.parse_results(output_file)
        data["contacts"] = contacts
        return data
