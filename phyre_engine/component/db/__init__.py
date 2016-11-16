from phyre_engine.component import Component
from phyre_engine.tools.mmcif import SimplifySelector
import gzip
import pathlib
import urllib.request
import Bio.SeqUtils
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB import PDBIO


class RCSBClusterDownload(Component):
    """Download a cluster file from the RCSB.

    The files retrieved by this module are described at
    `<http://www.rcsb.org/pdb/statistics/clusterStatistics.do>`, and are listed
    at `<http://www.rcsb.org/pdb/static.do?p=download/ftp/resources.jsp>`.

    This class adds the following keys to the pipeline data:

    ``cluster_file``:
        File containing clusters.
    """
    REQUIRED = []
    ADDS     = ['cluster_file']
    REMOVES  = []

    BASE_URL = "ftp://resources.rcsb.org/sequence/clusters/bc-{}.out"
    VALID_THRESHOLDS = {30, 40, 50, 70, 90, 95, 100}

    def __init__(self, threshold, filename="clusters"):
        """Initialise downloader at a given threshold.

        Args:
            ``threshold``: Sequence identity of clusters. Only the values
                listed at
                `<http://www.rcsb.org/pdb/static.do?p=download/ftp/resources.jsp>`
                are valid. Other values will cause an exception to be raised when
                `run` is called.
            ``filename``: Optional argument specifying the the filename at
                which to save the cluster file.

        Raises:
            ValueError: Invalid threshold value supplied.
        """
        err_msg = "Invalid threshold {}. Valid values: {}"
        try:
            self.threshold = int(threshold)
        except ValueError as e:
            raise ValueError(
                err_msg.format(threshold, self.VALID_THRESHOLDS)) from e

        if self.threshold not in self.VALID_THRESHOLDS:
            raise ValueError(err_msg.format(threshold, self.VALID_THRESHOLDS))

        self.filename = filename

    def run(self, data):
        """Download and parse the cluster file.

        Args:
            ``data``: Data carried through the pipeline.

        Raises:
            URLError: Error downloading the cluster file, most likely a result
                of specifying an invalid threshold in the constructor.
        """

        clus_file, headers = urllib.request.urlretrieve(
                self.BASE_URL.format(self.threshold),
                self.filename)
        data["cluster_file"] = clus_file
        return data

class ClusterParser(Component):
    """Parse a cluster file.

    The following keys are required when running this component.
    ``cluster_file``:
        File containing clusters. Each line of this file should represent one
        cluster, with each structure in that cluster separated by whitespace.


    The following keys are added when running this component:
    ``clusters``:
        Array of arrays containing PDB identifiers. Each sub-array is a
        cluster, with the ordering preserved from the original file. PDB IDs
        are of the form ``1XYZ_A``; that is, the PDB ID and chain ID separated
        by an underscore.
    """
    REQUIRED = ['cluster_file']
    ADDS     = ['clusters']
    REMOVES  = []

    def run(self, data):
        """Download and parse the cluster file.

        Args:
            ``data``: Data carried through the pipeline.
        """
        clus_file = self.get_vals(data)

        clusters = []
        with open(clus_file) as clus_fh:
            for clus_ln in clus_fh:
                clusters.append(clus_ln.rstrip().split())
        data["clusters"] = clusters
        return data

class SimpleRepresentativePicker(Component):
    """Simply pick the first element to represent each cluster.

    The following keys are required when running this component.
    ``clusters``: List of clusters. See `ClusterParser` for details on the data
        structure.

    The following keys are added when running this component:
    ``representatives``: Array of PDB IDs. Each ID consists of the PDB ID and
        the PDB chain, separated by an underscore.
    """
    REQUIRED = ['clusters']
    ADDS     = ['representatives']
    REMOVES  = []

    def run(self, data):
        """Extract representatives."""

        clusters = self.get_vals(data)
        representatives = [clus[0] for clus in clusters]
        data['representatives'] = representatives
        return data

class ChainPDBBuilder(Component):
    """For each representative structure, extract that chain to a PDB file file
    and extract the sequence (based on ATOM records) of that structure.

    This component will read the chain of each representative from the
    corresponding MMCIF file, write it to a PDB file, and store the sequence of
    the ATOM records in that PDB file.
    """

    REQUIRED = ["representatives"]
    ADDS     = ["templates"]
    REMOVES  = []

    def __init__(self, mmcif_dir, pdb_dir):
        """Initialise new component.

        Arguments:
            ``mmcif_dir``: Base directory of the MMCIF archive.
            ``pdb_dir``: Base directory in which to store PDB files.
        """
        self.mmcif_dir = pathlib.Path(mmcif_dir)
        self.pdb_dir = pathlib.Path(pdb_dir)

    def run(self, data):
        """Run the component."""
        representatives = self.get_vals(data)
        if "templates" not in data:
            data["templates"] = []

        mmcif_parser = MMCIFParser()
        pdb_parser = PDBParser()
        pdbio = Bio.PDB.PDBIO()
        simple_selector = SimplifySelector()

        for rep in representatives:
            id, chain = rep.split("_")
            id = id.lower()
            middle = id[1:3].lower()

            # Create dir in which to store pdb file
            (self.pdb_dir / middle).mkdir(exist_ok=True)
            mmcif_file = self.mmcif_dir / middle / "{}.cif.gz".format(id)
            pdb_file = self.pdb_dir / middle / "{}_{}.pdb".format(id, chain)

            with gzip.open(str(mmcif_file), "rt") as mmcif_fh:
                structure = mmcif_parser.get_structure(rep, mmcif_fh)
            struc_model = next(structure.get_models())
            struc_chain = struc_model[chain]

            # Set chain ID to " " while writing
            struc_chain.id = " "
            pdbio.set_structure(struc_chain)
            pdbio.save(str(pdb_file), simple_selector)

            # Get sequence of the residues. To be absolutely sure that our
            # sequence matches with the ATOM records, we will parse the PDB file
            # that we just wrote.
            pdb_struc = pdb_parser.get_structure(rep, str(pdb_file))
            pdb_res = pdb_struc.get_residues()
            pdb_seq = "".join([Bio.SeqUtils.seq1(r.get_resname()) for r in pdb_res])

            # Build a Bio.PDB.SeqRecord object containing this sequence.
            seq_description = {
                "PDB": id,
                "chain": chain,
            }
            bio_seq = SeqRecord(
                    id=rep,
                    description=json.dumps(seq_description),
                    seq=Seq(pdb_seq)
            )
            data["templates"].append(bio_seq)

        return data
