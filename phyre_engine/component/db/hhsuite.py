import phyre_engine.tools.hhsuite.tool as hh
from phyre_engine.component import Component
import Bio.SeqIO
import os
import tempfile
import fileinput
import pathlib
import logging

log = lambda: logging.getLogger(__name__)

class MSABuilder(Component):
    """
    Build an MSA using hhblits. The sequence is taken from the ``sequence``
    field of the pipeline state. Files are named according to the ``name`` field
    of the pipeline state.
    """

    REQUIRED = ["sequence", "name"]
    ADDS = ["a3m", "hhr"]
    REMOVES = []
    CONFIG_SECTION = "hhsuite"

    def __init__(
            self, database, *args, bin_dir=None, HHLIB=None, basedir=".",
            overwrite=False, **kwargs):

        self.bin_dir = bin_dir
        self.HHLIB = HHLIB
        self.basedir = pathlib.Path(basedir)
        self.overwrite = overwrite

        self.flags = args
        self.options = kwargs
        self.options["database"] = database

    def hhblits(self, template):
        msa_path = self.basedir / "a3m"
        hhr_path = self.basedir / "hhr"

        with tempfile.NamedTemporaryFile(
            prefix="query-", suffix=".fasta") as query_file:

            msa_name = "{}.a3m".format(template["name"])
            report_name = "{}.hhr".format(template["name"])
            msa_file = msa_path / msa_name
            hhr_file = hhr_path / report_name

            # No need to recreate
            if (not msa_file.exists()) or self.overwrite:
                print(
                    ">{name}\n{sequence}\n".format(**template),
                    file=query_file)
                query_file.flush()

                options = self.options.copy()
                options["input"] = query_file.name
                options["oa3m"] = msa_file
                options["output"] = hhr_file

                cmd = hh.hhblits(
                    (self.bin_dir, "hhblits"),
                    positional=None,
                    flags=self.flags,
                    options=options)
                hh.run(cmd, check=True, HHLIB=self.HHLIB)

            template["a3m"] = str(msa_file)
            template["hhr"] = str(hhr_file)

    def run(self, data, config=None, pipeline=None):
        templates = self.get_vals(data)
        (self.basedir / "a3m").mkdir(exist_ok=True)
        (self.basedir / "hhr").mkdir(exist_ok=True)

        self.hhblits(data)
        return data

class AddSecondaryStructure(Component):
    """Add predicted and actual secondary structure to MSAs."""

    REQUIRED = ["a3m"]
    ADDS = []
    REMOVES = []
    CONFIG_SECTION = "hhsuite"

    def __init__(self, HHLIB=None, **_kwargs):
        # We need to have _kwargs because we define CONFIG_SECTION and so might
        # be given extra parameters.

        if HHLIB is None and "HHLIB" not in os.environ:
            raise ValueError(
                "HHLIB not set as parameter or environment variable.")
        self.HHLIB = HHLIB

        hhlib = HHLIB if HHLIB is not None else os.environ["HHLIB"]
        self.addss = str(pathlib.Path(hhlib, "scripts/addss.pl"))

    def run(self, data, config=None, pipeline=None):
        a3m = self.get_vals(data)

        cmd_line = [self.addss, "-i", a3m]
        log().debug("Running %s", cmd_line)
        hh.run(cmd_line, check=True, HHLIB=self.HHLIB)

        return data



class HMMBuilder(Component):
    """Use hhmake to build an hhm file from an a3m file."""

    REQUIRED = ["a3m", "name"]
    ADDS = ["hhm"]
    REMOVES = []
    CONFIG_SECTION = "hhsuite"

    def __init__(
            self, *args, bin_dir=None, HHLIB=None, basedir=".", overwrite=False,
            **kwargs):

        self.bin_dir = bin_dir
        self.HHLIB = HHLIB
        self.basedir = pathlib.Path(basedir)
        self.overwrite = overwrite

        self.flags = args
        self.options = kwargs

    def run(self, data, config=None, pipeline=None):
        a3m, name = self.get_vals(data)
        hhm_path = self.basedir / "hhm"
        hhm_path.mkdir(exist_ok=True)

        # Generate hhm file
        hhm_name = "{}.hhm".format(name)
        hhm_file = hhm_path / hhm_name

        if (not hhm_file.exists()) or self.overwrite:
            options = self.options.copy()
            options["input"] = a3m
            options["output"] = hhm_file

            cmd_line = hh.hhmake(
                (self.bin_dir, "hhmake"),
                positional=None,
                flags=self.flags,
                options=options)
            hh.run(cmd_line, HHLIB=self.HHLIB, check=True)

        data["hhm"] = str(hhm_file)
        return data

class CS219Builder(Component):
    """Use cstranslate to build column state sequences for an a3m."""

    REQUIRED = ["a3m", "name"]
    ADDS = []
    REMOVES = []
    CONFIG_SECTION = "hhsuite"

    def __init__(
            self, *args, bin_dir=None, HHLIB=None, basedir=".", overwrite=False,
            **kwargs):

        if HHLIB is None and "HHLIB" not in os.environ:
            raise ValueError(
                "HHLIB not set as parameter or environment variable.")

        self.bin_dir = bin_dir
        self.HHLIB = HHLIB
        self.basedir = pathlib.Path(basedir)
        self.overwrite = overwrite

        # Set some default options
        hhlib = HHLIB if HHLIB is not None else os.environ["HHLIB"]
        data_dir = pathlib.Path(hhlib, "data")
        options = {
            "context-data": str(data_dir / "context_data.lib"),
            "alphabet": str(data_dir / "cs219.lib"),
            "pc-admix": 0.3,
            "pc-ali": 4,
            "informat": "a3m"
        }
        # Overwrite default options with explicit options if set
        options.update(kwargs)
        self.options = options

        self.flags = list(args)
        if "binary" not in args and "b" not in args:
            self.flags.append("binary")

    def run(self, data, config=None, pipeline=None):
        a3m, name = self.get_vals(data)
        cs219_path = self.basedir / "cs219"
        cs219_path.mkdir(exist_ok=True)

        # Generate cs219 file
        cs219_name = "{}.cs219".format(name)
        cs219_file = cs219_path / cs219_name

        if (not cs219_file.exists()) or self.overwrite:
            options = self.options.copy()
            options["outfile"] = cs219_file
            options["infile"] = a3m
            cmd_line = hh.cstranslate(
                (self.bin_dir, "cstranslate"),
                positional=None,
                flags=self.flags,
                options=options)
            hh.run(cmd_line, check=True, HHLIB=self.HHLIB)
        data["cs219"] = str(cs219_file)
        return data

class DatabaseBuilder(Component):
    """Build ffindex/ffdata files for an hhsuite database."""

    REQUIRED = ["templates"]
    ADDS = ["database"]
    REMOVES = ["templates"]
    CONFIG_SECTION = "hhsuite"

    def __init__(self, db_prefix, bin_dir=None, overwrite=False, basedir=".",
                 **_kwargs):
        """
        Initialise a new DatabaseBuilder component.

        :param str db_prefix: Prefix for database. The databases used by hhblits
            consist of multiple files, named like
            ``<prefix>_{a3m,hhm,cs219}.ff{index,data}``.
        :param str bin_dir: Optional directory containing the ``ffindex_build``
            executable.
        :param bool overwrite: If ``True``, delete existing database files.
            Otherwise, ``ffindex_build`` may be called on existing files.
        :param str basedir: Base directory in which to save the database files.
        """
        # _kwargs is only here because we are in the hhsuite CONFIG_SECTION and
        # might need to eat some arguments that we don't care about.

        self.db_prefix = db_prefix
        self.bin_dir = bin_dir
        self.overwrite = overwrite
        self.basedir = pathlib.Path(basedir)

    def run(self, data, config=None, pipeline=None):
        """Collect and index the files that form an hhsuite database."""
        templates = self.get_vals(data)

        # First, sort templates by sequence length.
        templates.sort(key=lambda t: len(t["sequence"]))

        # Make basedir if it doesn't exist.
        self.basedir.mkdir(parents=True, exist_ok=True)

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
                for template in templates:
                    print(template[file_type], file=index)
                index.flush()

                # Run ffindex_build using the the temp file as the list of files
                # to include in the DB.
                cmd_line = hh.ffindex_build(
                    (self.bin_dir, "ffindex_build"),
                    positional=[ffdata, ffindex],
                    flags=["sort"],
                    options={"file_list": index.name})
                log().debug("Running command %s", cmd_line)
                hh.run(cmd_line, check=True)
                ff_dbs[file_type] = db_name

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
