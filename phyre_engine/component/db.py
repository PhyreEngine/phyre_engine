import phyre_engine.tools.hhsuite.tool as hh
from phyre_engine.component import Component
import copy
import gzip
import pathlib
import urllib.request
import Bio.SeqUtils
import Bio.SeqIO
import json
import os
import sys
import subprocess
import tempfile
import fileinput
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.PDB
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.PDBExceptions import PDBConstructionException

class NameTemplate(Component):
    """Add a ``name`` attribute to each template."""
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, name_fn):
        """
        :param function name_fn: Function accepting a template dictinoary and
            returning a name.
        """
        self.name_fn = name_fn

    def run(self, data):
        templates = self.get_vals(data)
        for template in templates:
            template["name"] = self.name_fn(template)
        return data

class DescribeTemplate(Component):
    """
    Alter the description of the template sequence.

    This component alters the ``.description`` field of the
    :py:class:`Bio.PDB.SeqRecord` object pointed to by the ``sequence`` key of
    each template.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, description_fn):
        """
        :param function name_fn: Function accepting a template dictinoary and
            returning a description.
        """
        self.description_fn = description_fn

    def run(self, data):
        templates = self.get_vals(data)
        for template in templates:
            template["sequence"].description = self.description_fn(template)
        return data

class RCSBClusterDownload(Component):
    """Download a cluster file from the RCSB.

    The files retrieved by this module are described at
    http://www.rcsb.org/pdb/statistics/clusterStatistics.do, and are listed
    at http://www.rcsb.org/pdb/static.do?p=download/ftp/resources.jsp.

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

        :param int threshold: Sequence identity of clusters. Only the values
            listed at
            `<http://www.rcsb.org/pdb/static.do?p=download/ftp/resources.jsp>`
            are valid. Other values will cause an exception to be raised when
            `run` is called.
        :param str filename: Optional argument specifying the the filename at
            which to save the cluster file.

        :raises ValueError: Invalid threshold value supplied.
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

        :param data: Data carried through the pipeline.

        :raises URLError: Error downloading the cluster file, most likely a result
            of specifying an invalid threshold in the constructor.
        """

        clus_file, _ = urllib.request.urlretrieve(
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

        :param data: Data carried through the pipeline.
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

    Sometimes it is useful to preserve the arcane mappings used by the authors
    of PDB structures. For our purposes, we are often mapping between sequence
    and structure, so it is useful to treat PDB files as arrays of residues. To
    preserve the mappings, then, we write a JSON-encoded array of residue IDs
    to a "map" file.
    """

    REQUIRED = ["templates"]
    ADDS     = []
    REMOVES  = []

    def __init__(
        self, mmcif_dir, pdb_dir, map_dir,
        mutation_selector, microhet_selector=None):
        """Initialise new component.

        :param str mmcif_dir: Base directory of the MMCIF archive.
        :param str pdb_dir: Base directory in which to store PDB files.
        :param str map_dir: Base directory in which sequence map files will be
            saved.
        :param phyre_engine.component.db.conformation.MutationSelector mutation_selector:
            Subclass of `phyre_engine.component.db.conformation.MutationSelector`
            used for selecting which mutation of a template should be used.
        :param phyre_engine.component.db.conformation.MicroConformationSelector microhet_selector:
            Subclass of
            `phyre_engine.component.db.conformation.MicroConformationSelector`
            used to pick a single conformation. All conformations will be kept
            if this is not supplied.

        """
        self.mmcif_dir = pathlib.Path(mmcif_dir)
        self.pdb_dir = pathlib.Path(pdb_dir)
        self.map_dir = pathlib.Path(map_dir)
        self.mutation_selector = mutation_selector
        self.microhet_selector = microhet_selector

    def run(self, data):
        """Run the component."""
        templates = self.get_vals(data)

        # Create a new list rather than modifying the existing list in place.
        # This is so that we can handle failures gracefully, if you can call
        # squawking an exception and then ignoring the template graceful.
        new_templates = []
        for template in templates:
            pdb_id = template["PDB"].lower()
            chain = template["chain"]
            middle = pdb_id[1:3].lower()

            # Create dir in which to store pdb file
            pdb_dir = self.pdb_dir / middle
            map_dir = self.map_dir / middle
            pdb_dir.mkdir(exist_ok=True)
            map_dir.mkdir(exist_ok=True)

            mmcif_file = self.mmcif_dir / middle / "{}.cif.gz".format(pdb_id)

            # If the output files already exist, we can just read from them.
            try:
                (pdb_files, map_files) = self.find_template_files(pdb_id, chain)
                if pdb_files and map_files:
                    for pdb_file, map_file in zip(pdb_files, map_files):
                        new_template = copy.copy(template)
                        self.read_chain(
                            pdb_id, chain, new_template, pdb_file, map_file)
                        new_templates.append(new_template)
                else:
                    extracted_templates = self.extract_chain(
                        pdb_id, chain, template,
                        mmcif_file)
                    new_templates.extend(extracted_templates)
            except PDBConstructionException as e:
                # TODO: When we start to use a logging framework, log this error
                # properly.
                print(
                    "Error parsing chain {} of {}".format(chain, pdb_id),
                    file=sys.stderr)
            except Exception as e:
                raise ChainExtractionError(pdb_id, chain) from e

        data["templates"] = new_templates
        return data

    def find_template_files(self, pdb_id, chain_id):
        pdb_id = pdb_id.lower()
        middle = pdb_id[1:3]
        pdb_dir = self.pdb_dir / middle
        map_dir = self.map_dir / middle

        # First, check if the correct files exist for a single conformation.
        # Next, check if multiple conformations exist.
        # Otherwise, return None
        pdb_file = pdb_dir / "{}_{}.pdb".format(pdb_id, chain_id)
        map_file = map_dir / "{}_{}.json".format(pdb_id, chain_id)
        if pdb_file.exists() and map_file.exists():
            return ([pdb_file], [map_file])
        else:
            pdb_files = sorted(
                list(pdb_dir.glob("{}_{}-*.pdb")),
                key=lambda p: p.name)
            map_files = sorted(
                list(map_dir.glob("{}_{}-*.json")),
                key=lambda p: p.name)
            if pdb_files and map_files:
                return (pdb_files, map_files)
                #TODO: Throw error if lists don't match
            else:
                return (None, None)

    def read_chain(
            self, pdb_id, chain_id, template,
            pdb_file, map_file):
        """
        Read a previously-extracted chain.

        Sets the ``sequence``, ``structure`` and ``map`` fields of the
        ``template`` argument. Modifies ``template`` in place.

        :param str pdb_id: PDB ID.
        :param dict template: Dictionary describing template.
        :param pathlib.Path pdb_file: Path object pointing to the PDB file into
            which the chain will be saved.
        :param pathlib.Path map_file: Path object poniting to the file into
            which the residue mapping will be saved.
        """
        with pdb_file.open("r") as pdb_fh:
            structure = PDBParser().get_structure(pdb_id, pdb_fh)
            model = next(structure.get_models())
            chain = model[' ']

        # Get sequence of the residues.
        pdb_seq = "".join([
            Bio.SeqUtils.seq1(r.get_resname())
            for r in chain
        ])

        # Build a Bio.PDB.SeqRecord object containing this sequence.
        seq_name = "{}_{}".format(pdb_id, chain_id)
        bio_seq = SeqRecord(
                id=seq_name,
                description="",
                seq=Seq(pdb_seq)
        )
        template["name"] = seq_name
        template["sequence"] = bio_seq
        template["structure"] = pdb_file
        template["map"] = map_file


    def extract_chain(
            self, pdb_id, chain, template,
            mmcif_file):
        """
        Extract a chain from an MMCIF file.

        This method will extract the chain, save a PDB file containing that
        chain, and write a map file containing the map between residue index and
        author-assigned residue IDs.

        :return: A list of templates, shallow-copied from the original template,
            and with "sequence", "structure" and "map" keys added.

        :param str pdb_id: PDB ID.
        :param str chain: PDB chain.
        :param str template: Dictionary describing template.
        :param pathlib.Path mmcif_file: Path object pointing to the MMCIF from
            which the chain  will be extracted.
        """
        mmcif_parser = MMCIFParser()
        pdbio = Bio.PDB.PDBIO()

        pdb_id = pdb_id.lower()
        middle = pdb_id[1:3]
        pdb_dir = self.pdb_dir / middle
        map_dir = self.map_dir / middle

        with gzip.open(str(mmcif_file), "rt") as mmcif_fh:
            structure = mmcif_parser.get_structure(pdb_id, mmcif_fh)
        struc_model = next(structure.get_models())

        templates = []
        mutations = self.mutation_selector.select(struc_model[chain])
        for index, conformation in enumerate(mutations):
            if self.microhet_selector:
                conformation = self.microhet_selector.select(conformation)
            struc_chain, res_map = self.sanitise_chain(conformation)

            if len(mutations) > 1:
                pdb_file = pdb_dir / "{}_{}-{}.pdb".format(pdb_id, chain, index+1)
                map_file = map_dir / "{}_{}-{}.json".format(pdb_id, chain, index+1)
            else:
                pdb_file = pdb_dir / "{}_{}.pdb".format(pdb_id, chain)
                map_file = map_dir / "{}_{}.json".format(pdb_id, chain)

            # Build mapping between sequence index and residue ID. This is just
            # an array of residue IDs that can be indexed in the same way as
            # everything else.
            with map_file.open("w") as map_fh:
                json.dump(res_map, map_fh)

            # Write PDB file.
            pdbio.set_structure(struc_chain)
            pdbio.save(str(pdb_file))

            # Get sequence of the residues.
            pdb_seq = "".join([
                Bio.SeqUtils.seq1(r.get_resname())
                for r in struc_chain
            ])

            # Build a Bio.PDB.SeqRecord object containing this sequence.
            bio_seq = SeqRecord(
                    pdb_id="{}_{}".format(pdb_id, chain),
                    description="",
                    seq=Seq(pdb_seq)
            )
            new_template = copy.copy(template)
            new_template["sequence"] = bio_seq
            new_template["structure"] = pdb_file
            new_template["map"] = map_file
            templates.append(new_template)
        return templates

    def sanitise_chain(self, chain, new_id=" "):
        """Strip insertion codes and disordered atoms from a chain.

        Given a `Bio.PDB` structure (`Bio.PDB.Chain.Chain`), renumber residues
        from 1, stripping out insertion codes. At the same time, we remove
        disordered atoms and residues. For now, we just keep the first
        conformation that occurs.

        :param `Bio.PDB.Chain` chain: structure to sanitise.
        :param str new_id: ID of the new chain.

        :return: A tuple containing:

            1. A new `Bio.PDB.Chain` object consisting of sanitised residues.
            2. An array containing the author-assigned IDs of the sanitised residues.
        """
        mapping = []
        sanitised_chain = Chain(new_id)
        res_index = 1
        for res in chain:
            # Discard HETATMs
            if res.get_id()[0] != ' ':
                continue

            sanitised_res = Residue(
                (' ', res_index, ' '),
                res.get_resname(),
                res.get_segid())
            mapping.append(res.get_id())

            for atom in res:
                sanitised_res.add(atom.copy())
            sanitised_chain.add(sanitised_res)
            res_index += 1
        return sanitised_chain, mapping

class ChainExtractionError(Exception):
    """
    Exeception thrown when an error is encountered extracting a chain from an
    MMCIF file.
    """

    ERR_MSG = "Error extracting chain {} from structure {}"
    def __init__(self, structure, chain):
        super().__init__(self.ERR_MSG.format(chain, structure))


class MSABuilder(Component):
    """Build MSA for each template sequence using hhblits."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(
            self, hhblits_db, overwrite=False,
            basedir=".", **hhblits_args):
        """Initialise a new MSA builder; this is essentially a parallel hhblits
        component.

        This differs slightly from ``phyre_engine.component.hhsuite.HHBlits``:
        that class will build a profile for a single sequence, while this will
        operate on an array of sequences.

        Files will be placed in the directory `basedir`. Subdirectories under
        `basedir` named ``a3m`` and ``hhr`` will be created containing the MSAs
        and hhblits reports, respectively.

        This component reads sequences (Bio.SeqRecord objects) from the
        ``sequences`` key of the pipeline data, and adds the ``msas`` and
        ``reports`` keys.
        """
        self.hhblits_db   = hhblits_db
        self.hhblits_args = hhblits_args
        self.overwrite = overwrite
        self.basedir = pathlib.Path(basedir)

    def run(self, data):
        templates = self.get_vals(data)

        msa_path = self.basedir / "a3m"
        hhr_path = self.basedir / "hhr"

        msa_path.mkdir(exist_ok=True)
        hhr_path.mkdir(exist_ok=True)

        for template in templates:
            sequence = template["sequence"]

            with tempfile.NamedTemporaryFile(suffix=".fasta") as query_file:
                msa_name    = "{}.a3m".format(template["name"])
                report_name = "{}.hhr".format(template["name"])
                msa_file = msa_path / msa_name
                hhr_file = hhr_path / report_name

                # No need to recreate
                if (not msa_file.exists()) or self.overwrite:
                    Bio.SeqIO.write(sequence, query_file.name, "fasta")
                    hhblits = hh.HHBlits(
                            database=self.hhblits_db,
                            input=query_file.name,
                            output=hhr_file,
                            oa3m=msa_file,
                            **self.hhblits_args)
                    hhblits.run()

                template["a3m"] = str(msa_file)
                template["hhr"] = str(hhr_file)
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

    def __init__(self, overwrite=False, basedir=".", **hhmake_args):
        """Initialise a new HMMBuilder

        :param bool overwrite: If true, force an overwrite of existing files.
            Otherwise, files that already exist will not be touched.
        :param str basedir: Base directory in which to store HMM files. Files
            will be stored in the subdirectory ``hhm`` below this.
        :param **hhmake_args: Extra arguments to pass to
            :py:class:`phyre_engine.tools.hhsuite.HHMake`.
        """
        self.overwrite = overwrite
        self.basedir = pathlib.Path(basedir)
        self.hhmake_args = hhmake_args

    def run(self, data):
        templates = self.get_vals(data)
        hhm_path = self.basedir / "hhm"
        hhm_path.mkdir(exist_ok=True)

        for template in templates:
            # Generate hhm file
            hhm_name = "{}.hhm".format(template["name"])
            hhm_file = hhm_path / hhm_name

            if (not hhm_file.exists()) or self.overwrite:
                hhmake = hh.HHMake(
                    template["a3m"],
                    output=hhm_file,
                    **self.hhmake_args
                )
                hhmake.run()
            template["hhm"] = str(hhm_file)
        return data

class CS219Builder(Component):
    """Use cstranslate to build column state sequences for each sequence."""

    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, overwrite=False, basedir=".", **cstranslate_args):
        """Initialise a new CS219Builder.

        :param bool overwrite: Overwrite any existing cs219 files.
        :param str basedir: Base directory in which to store files. Files will
            be stored in the subdirectory ``cs219`` below this directory.
        :param **cstranslate_args: Extra arguments to pass to
            :py:class:`phyre_engine.tools.hhsuite.CSTranslate`.
        """
        self.overwrite = overwrite
        self.basedir = pathlib.Path(basedir)
        self.cstranslate_args = cstranslate_args

    def run(self, data):
        hhlib = os.environ["HHLIB"]

        templates = self.get_vals(data)
        cs219_path = self.basedir / "cs219"
        cs219_path.mkdir(exist_ok=True)

        for template in templates:
            # Generate cs219 file
            cs219_name = "{}.cs219".format(template["name"])
            cs219_file = cs219_path / cs219_name

            if (not cs219_file.exists()) or self.overwrite:
                cs_args = {
                    "-A": pathlib.Path(hhlib, "data/cs219.lib"),
                    "-D": pathlib.Path(hhlib, "data/context_data.lib"),
                    "-x": 0.3, "-c": 4, "-b": True,
                    "-I": "a3m",
                    "-i": template["a3m"],
                    "-o": cs219_file
                }
                cs_args.update(self.cstranslate_args)
                cstranslate = hh.CSTranslate(**cs_args)
                cstranslate.run()
            template["cs219"] = str(cs219_file)
        return data

class DatabaseBuilder(Component):
    """Build ffindex/ffdata files for an hhsuite database."""

    REQUIRED = ["templates"]
    ADDS = ["database"]
    REMOVES = ["templates"]

    def __init__(
        self,
        db_prefix, overwrite=False, basedir=".", **ffindex_args):
        """
        Initialise a new DatabaseBuilder component.

        :param str db_prefix: Prefix for database. The databases used by hhblits
            consist of multiple files, named like
            ``<prefix>_{a3m,hhm,cs219}.ff{index,data}``.
        :param bool overwrite: If ``True``, delete existing database files.
            Otherwise, ``ffindex_build`` may be called on existing files.
        :param str basedir: Base directory in which to save the database files.
        :param **ffindex_args: Extra arguments to pass to
            :py:class:`phyre_engine.tools.hhsuite.FFIndexBuild`.
        """
        self.db_prefix = db_prefix
        self.overwrite = overwrite
        self.basedir = basedir
        self.ffindex_args = ffindex_args

    def run(self, data):
        """Collect and index the files that form an hhsuite database."""
        templates = self.get_vals(data)

        # Make basedir if it doesn't exist.
        pathlib.Path(self.basedir).mkdir(parents=True, exist_ok=True)

        # Collect a3m/hhm/cs219 files into ffindex/ffdata databases
        to_collect = ["a3m", "hhm", "cs219"]
        ff_dbs = {}
        db_prefix = pathlib.Path(self.basedir, self.db_prefix)

        for file_type in to_collect:
            db_name = pathlib.Path("{}_{}".format(str(db_prefix), file_type))
            ffindex = pathlib.Path("{}.ffindex".format(str(db_name)))
            ffdata = pathlib.Path("{}.ffdata".format(str(db_name)))

            if self.overwrite:
                if ffindex.exists():
                    ffindex.unlink()
                if ffdata.exists():
                    ffdata.unlink()

            with tempfile.NamedTemporaryFile("w") as index:
                # Write all files of file_type `file_type` to a temp file
                for file in [template[file_type] for template in templates]:
                    print(file, file=index)
                index.flush()

                # Run ffindex_build using the the temp file as the list of files
                # to incude in the DB.
                ffindex_builder = hh.FFIndexBuild(
                    ffdata, ffindex,
                    file_list=index.name,
                    **self.ffindex_args)
                ffindex_builder.run()
                ff_dbs[file_type] = db_name

            # Sort the indices
            ffindex_sorter = hh.FFIndexBuild(
                "{}.ffdata".format(ff_dbs[file_type]),
                "{}.ffindex".format(ff_dbs[file_type]),
                append=True, sort=True,
                **self.ffindex_args)
            ffindex_sorter.run()

        # Cut useless information from the indices of each file.
        for ff_db in ff_dbs.values():
            self._trim_index_names(ff_db)

        del data["templates"]
        data["database"] = str(db_prefix)
        return data


    def _trim_index_names(self, db):
        """Replace names of components in an ffindex with their stem.

        When we generate an index using ``ffindex_build``, the names of each
        element in that database are simply taken from the filenames used as
        input. Since we are trying to build a valid database for hhblits and
        hhsearch, each name in the index must be consistent across database. To
        do that, we use the name of each template as the identifier.
        """

        index = "{}.ffindex".format(db)
        with fileinput.input(index, inplace=True) as fh:
            for line in fh:
                fields = line.split("\t")
                fields[0] = pathlib.PurePath(fields[0]).stem
                print("\t".join(fields), end="")
