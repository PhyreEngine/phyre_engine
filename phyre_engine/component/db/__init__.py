import phyre_engine.tools.hhsuite as hh
from phyre_engine.component import Component
from phyre_engine.tools.mmcif import SimplifySelector
import gzip
import pathlib
import urllib.request
import Bio.SeqUtils
import Bio.SeqIO
import json
import os
import subprocess
import tempfile
import fileinput
import shutil
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
    ``templates``: Array of templates. Each template is a dictionary. This
        component sets the keys ``PDB`` and ``chain`` for each template,
        corresponding to the PDB ID and PDB chain of the cluster
        representatives.
    """
    REQUIRED = ['clusters']
    ADDS     = ['templates']
    REMOVES  = []

    def run(self, data):
        """Extract representatives."""

        clusters = self.get_vals(data)
        representatives = [clus[0] for clus in clusters]
        if "templates" not in data:
            data["templates"] = []

        for rep in representatives:
            (pdb, chain) = rep.split("_")
            data["templates"].append({"PDB":pdb, "chain":chain})
        return data

class ChainPDBBuilder(Component):
    """For each structure, extract that chain to a PDB file file
    and extract the sequence (based on ATOM records) of that structure.

    This component will read the chain of each representative from the
    corresponding MMCIF file, write it to a PDB file, and store the sequence of
    the ATOM records in that PDB file.
    """

    REQUIRED = ["templates"]
    ADDS     = []
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
        templates = self.get_vals(data)

        mmcif_parser = MMCIFParser()
        pdb_parser = PDBParser()
        pdbio = Bio.PDB.PDBIO()
        simple_selector = SimplifySelector()

        for template in templates:
            id    = template["PDB"].lower()
            chain = template["chain"]
            middle = id[1:3].lower()

            # Create dir in which to store pdb file
            (self.pdb_dir / middle).mkdir(exist_ok=True)
            mmcif_file = self.mmcif_dir / middle / "{}.cif.gz".format(id)
            pdb_file = self.pdb_dir / middle / "{}_{}.pdb".format(id, chain)

            with gzip.open(str(mmcif_file), "rt") as mmcif_fh:
                structure = mmcif_parser.get_structure(id, mmcif_fh)
            struc_model = next(structure.get_models())
            struc_chain = struc_model[chain]

            # Set chain ID to " " while writing
            struc_chain.id = " "
            pdbio.set_structure(struc_chain)
            pdbio.save(str(pdb_file), simple_selector)

            # Get sequence of the residues. To be absolutely sure that our
            # sequence matches with the ATOM records, we will parse the PDB file
            # that we just wrote.
            pdb_struc = pdb_parser.get_structure(id, str(pdb_file))
            pdb_res = pdb_struc.get_residues()
            pdb_seq = "".join([Bio.SeqUtils.seq1(r.get_resname()) for r in pdb_res])

            # Build a Bio.PDB.SeqRecord object containing this sequence.
            bio_seq = SeqRecord(
                    id="{}_{}".format(id, chain),
                    description="",
                    seq=Seq(pdb_seq)
            )
            template["sequence"] = bio_seq
            template["structure"] = pdb_file

        return data

class MSABuilder(Component):
    """Build MSA for each template sequence using hhblits."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, hhblits_db, overwrite=False, **hhblits_args):
        """Initialise a new MSA builder; this is essentially a parallel hhblits
        component.

        This differs slightly from ``phyre_engine.component.hhsuite.HHBlits``:
        that class will build a profile for a single sequence, while this will
        operate on an array of sequences.

        Files will be placed in the current working directory. Subdirectories
        named ``a3m`` and ``hhr`` will be created containing the MSAs and
        hhblits reports, respectively.

        This component reads sequences (Bio.SeqRecord objects) from the
        ``sequences`` key of the pipeline data, and adds the ``msas`` and
        ``reports`` keys.
        """
        self.hhblits_db   = hhblits_db
        self.hhblits_args = hhblits_args
        self.overwrite = overwrite

    def run(self, data):
        templates = self.get_vals(data)

        msa_path = pathlib.Path("a3m")
        hhr_path = pathlib.Path("hhr")

        msa_path.mkdir(exist_ok=True)
        hhr_path.mkdir(exist_ok=True)

        for template in templates:
            sequence = template["sequence"]
            seq_name = "{}_{}".format(template["PDB"], template["chain"])

            with tempfile.NamedTemporaryFile(suffix=".fasta") as query_file:
                msa_name    = "{}.a3m".format(seq_name)
                report_name = "{}.hhr".format(seq_name)
                msa_file = msa_path / msa_name
                hhr_file = hhr_path / report_name

                # No need to recreate
                if (not msa_file.exists()) or self.overwrite:
                    # Build a machine-readable (JSON) description for the sequence
                    desc_dict = {k: template[k] for k in ("PDB", "chain")}
                    sequence.description = json.dumps(desc_dict)

                    Bio.SeqIO.write(sequence, query_file.name, "fasta")
                    hhblits = hh.HHBlits(
                            database=self.hhblits_db,
                            input=query_file.name,
                            output=str(hhr_file),
                            oa3m=str(msa_file),
                            **self.hhblits_args)
                    hhblits.run()

                template["a3m"] = str(msa_file)
                template["hhm"] = str(hhr_file)
                template["name"] = seq_name
        return data

class AddSecondaryStructure(Component):
    """Add predicted and actual secondary structure to MSAs."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def run(self, data):
        templates = self.get_vals(data)

        for template in templates:
            # Add secondary structure
            hhlib = os.environ["HHLIB"]
            addss_pl = pathlib.Path(hhlib, "scripts/addss.pl")
            subprocess.run([str(addss_pl), "-i", template["a3m"]])
        return data

class HMMBuilder(Component):
    """Use hhmake to build an hhm file for each template."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, overwrite=False):
        self.overwrite = overwrite

    def run(self, data):
        templates = self.get_vals(data)
        hhm_path = pathlib.Path("hhm")
        hhm_path.mkdir(exist_ok=True)

        for template in templates:
            # Generate hhm file
            hhm_name = "{}.hhm".format(template["name"])
            hhm_file = hhm_path / hhm_name

            if (not hhm_file.exists()) or self.overwrite:
                subprocess.run(["hhmake",
                    "-i", template["a3m"],
                    "-o", str(hhm_file)
                ])
            template["hhm"] = str(hhm_file)
        return data

class CS219Builder(Component):
    """Use cstranslate to build coolumn state sequences for each sequence."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, overwrite=False):
        self.overwrite = overwrite

    def run(self, data):
        hhlib = os.environ["HHLIB"]

        templates = self.get_vals(data)
        cs219_path = pathlib.Path("cs219")
        cs219_path.mkdir(exist_ok=True)

        for template in templates:
            # Generate cs219 file
            cs219_name = "{}.cs219".format(template["name"])
            cs219_file = cs219_path / cs219_name

            if (not cs219_file.exists()) or self.overwrite:
                subprocess.run([
                    "cstranslate",
                    "-A", str(pathlib.Path(hhlib, "data/cs219.lib")),
                    "-D", str(pathlib.Path(hhlib, "data/context_data.lib")),
                    "-x", str(0.3), "-c", str(4), "-b",
                    "-I", "a3m",
                    "-i", template["a3m"],
                    "-o", str(cs219_file)
                ])
            template["cs219"] = str(cs219_file)
        return data

class DatabaseBuilder(Component):
    """Build ffindex/ffdata files for an hhsuite database."""

    REQUIRED = ["templates"]
    ADDS = ["database"]
    REMOVES = ["templates"]

    def __init__(self, db_prefix, overwrite=False):
        self.db_prefix = db_prefix
        self.overwrite = overwrite

    def run(self, data):
        """Collect and index the files that form an hhsuite database."""
        templates = self.get_vals(data)

        # Collect a3m/hhm/cs219 files into ffindex/ffdata databases
        to_collect = ["a3m", "hhm", "cs219"]
        ff_dbs = {}
        for type in to_collect:
            db_name = "{}_{}".format(self.db_prefix, type)
            if self.overwrite:
                ffindex = pathlib.Path("{}.ffindex".format(db_name))
                ffdata  = pathlib.Path("{}.ffdata".format(db_name))
                if ffindex.exists():
                    ffindex.unlink()
                if ffdata.exists():
                    ffdata.unlink()

            with tempfile.NamedTemporaryFile("w") as index:
                # Write all files of type `type` to a temp file
                for file in [template[type] for template in templates]:
                    print(file, file=index)
                index.flush()

                # Run ffindex_build using the the temp file as the list of files
                # to incude in the DB.
                subprocess.run(["ffindex_build",
                    "{}.ffdata".format(db_name),
                    "{}.ffindex".format(db_name),
                    "-f", index.name])
                ff_dbs[type] = db_name

        # Sort the cs219 ffindex
        subprocess.run([
            "ffindex_build", "-as",
            "{}.ffdata".format(ff_dbs["cs219"]),
            "{}.ffindex".format(ff_dbs["cs219"])])

        # Cut useless information from the indices of each file.
        for type, ff_db in ff_dbs.items():
            self._trim_index_names(templates, type, ff_db)

        # Get column state lengths and sort the a3m and hhm databases according
        # to that order.
        with open("{}.ffindex".format(ff_dbs["cs219"]), "r") as idx:
            lines = [line.split() for line in idx.readlines()]
            # Sort by length (index 2), and return the name (index 0)
            names = [ln[0] for ln in sorted(lines, key=lambda ln: int(ln[2]))]

            with tempfile.NamedTemporaryFile("w") as ordered_tmp:
                for n in names:
                    print(n, file=ordered_tmp)
                ordered_tmp.flush()
                self._ffindex_order(ff_dbs["a3m"], ordered_tmp.name)
                self._ffindex_order(ff_dbs["hhm"], ordered_tmp.name)

        del data["templates"]
        data["database"] = self.db_prefix


    def _trim_index_names(self, templates, type, db):
        """Replace names of components in an ffindex with their stem.

        When we generate an index using ``ffindex_build``, the names of each
        element in that database are simply taken from the filenames used as
        input. Since we are trying to build a valid database for hhblits and
        hhsearch, each name in the index must be consistent across database. To
        do that, we use the name of each template as the identifier.
        """

        # Build a dictionary of file -> name
        file_to_name = {t[type]: t["name"] for t in templates}

        index = "{}.ffindex".format(db)
        with fileinput.input(index, inplace=True) as fh:
            for line in fh:
                fields = line.split("\t")
                fields[0] = file_to_name[str(pathlib.Path(type, fields[0]))]
                print("\t".join(fields), end="")

    def _ffindex_order(self, database, by):
        """Run ffindex_order to order an ffindex/ffdata database according with
        the same order as the list given in the file "by". This function will
        overwrite the unsorted database.
        """
        subprocess.run([
            "ffindex_order", by,
            "{}.ffdata".format(database),
            "{}.ffindex".format(database),
            "{}_ordered.ffdata".format(database),
            "{}_ordered.ffindex".format(database)
        ])
        shutil.move("{}_ordered.ffdata".format(database),
                "{}.ffdata".format(database))
        shutil.move("{}_ordered.ffindex".format(database),
                "{}.ffindex".format(database))
