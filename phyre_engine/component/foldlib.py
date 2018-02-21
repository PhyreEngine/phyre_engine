"""
This module contains the components that are specific to the process of
building the fold library.

"""
import collections
import collections.abc
import copy
import datetime
from pathlib import Path
import urllib.parse
import urllib.request
import xml.etree.ElementTree

import Bio.PDB.PDBParser

import phyre_engine.component.component
from phyre_engine.component import Component
import phyre_engine.component.hhsuite
import phyre_engine.component.pdb
import phyre_engine.component.secstruc
import phyre_engine.pipeline
from phyre_engine.tools.template import Template, TemplateDatabase
import phyre_engine.tools.pdb

class Open(Component):
    """
    Open the fold library for reading or writing.

    This component will open the fold library database as a persistent
    "connection". The :py:class:`phyre_engine.tools.template.TemplateDatabase`
    object will be stored in the ``template_db`` key in the pipeline state.  To
    close the connection, use either the :py:class:`~.Commit` or
    :py:class:`~.Rollback` components.

    If the `write` parameter is `True`, then the database is opened
    exclusively, meaning that other processes cannot read from it or write to
    it.

    :param str template_db: Location of the template database SQLite file.
    :param str chain_dir: Root directory containing template files.
    """
    ADDS = ["template_db"]
    REMOVES = []
    REQUIRED = []

    @classmethod
    def config(cls, params, config):
        return config.extract(
            {"foldlib": ["template_db", "chain_dir"]}
        ).merge_params(params)

    def __init__(self, template_db, chain_dir, write=False, trace=False):
        self.template_db = template_db
        self.chain_dir = chain_dir
        self.write = write
        self.trace = trace

    def run(self, data, config=None, pipeline=None):
        """Open connection to fold library database."""
        trace_callback = self.logger.debug if self.trace else None
        template_db = TemplateDatabase(self.template_db, self.chain_dir,
                                       self.write, trace_callback)
        data["template_db"] = template_db
        return data

class Commit(Component):
    """
    Commit all changes to the fold library.

    This will remove the ``template_db`` key from the pipeline state.
    """
    ADDS = []
    REMOVES = ["template_db"]
    REQUIRED = ["template_db"]

    def run(self, data, config=None, pipeline=None):
        """Commit changes to the template database."""
        template_db = self.get_vals(data)
        template_db.commit()
        del data["template_db"]
        return data

class Rollback(Component):
    """
    Roll back all changes to the template database.

    This will remove the ``template_db`` key from the pipeline state.
    """
    ADDS = []
    REMOVES = ["template_db"]
    REQUIRED = ["template_db"]

    def run(self, data, config=None, pipeline=None):
        """Commit changes to the template database."""
        template_db = self.get_vals(data)
        template_db.rollback()
        del data["template_db"]
        return data


class UpdateMetadata(Component):
    """Set the latest update date of the template library."""
    REQUIRED = ["template_db"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Set update date of template library."""
        template_db = self.get_vals(data)
        template_db.updated = datetime.date.today()
        return data


class Metadata(Component):
    """
    Look up template metadata in the template library based on the ``PDB``
    field, and set the ``metadata`` field.

    The items in the ``metadata`` dictionary are described in
    :py:meth:`phyre_engine.tools.template.TemplateDatabase.add_pdb`.
    """

    REQUIRED = ["template_db", "PDB"]
    ADDS = ["metadata"]
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Retrieve metadata from the fold library."""
        template_db, pdb_id = self.get_vals(data)
        data["metadata"] = template_db.get_pdb(pdb_id)
        return data


class Map(phyre_engine.component.component.Map):
    """
    Subclass of :py:class:`phyre_engine.component.component.Map` that
    automatically iterates over the ``templates`` list and transfers the
    ``template_db`` component into each template.
    """
    REQUIRED = ["templates", "template_db"]

    def __init__(self, **kwargs):
        super().__init__(field="templates", **kwargs)

    def run(self, data, config=None, pipeline=None):
        """Iterating over each template while building fold library."""
        for template in data["templates"]:
            template["template_db"] = data["template_db"]

        data = super().run(data, config, pipeline)

        for template in data["templates"]:
            del template["template_db"]
        return data


class RetrieveNewPDBs(Component):
    """
    Retrieve a list of PDB entries that have either been released or updated
    since the last revision of the fold library. If the fold library has never
    been updated, all PDB IDs will be returned.

    Results are stored in the ``templates`` array as dictionaries with a
    ``PDB`` field.

    :param str override_date: Override the date in the fold library. All PDBs
        since this date (in ``YY-MM-DD`` format) will be retrieved.
    """
    ADDS = ["templates"]
    REQUIRED = ["template_db"]
    REMOVES = []


    SEARCH_URL = "https://www.rcsb.org/pdb/rest/search"
    CURRENT_ENTRIES_URL = "https://www.rcsb.org/pdb/rest/getCurrent"

    XML_SEARCH = """\
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ReleaseDateQuery</queryType>
        <description>Released on 2018-01-10 and later</description>
        <pdbx_audit_revision_history.revision_date.comparator>between</pdbx_audit_revision_history.revision_date.comparator>
        <pdbx_audit_revision_history.revision_date.min>{date}</pdbx_audit_revision_history.revision_date.min>
        <pdbx_audit_revision_history.ordinal.comparator>=</pdbx_audit_revision_history.ordinal.comparator>
        <pdbx_audit_revision_history.ordinal.value>1</pdbx_audit_revision_history.ordinal.value>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>or</conjunctionType>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ReviseDateQuery</queryType>
        <description>Revised on 2018-01-10 and later</description>
        <pdbx_audit_revision_history.revision_date.comparator>between</pdbx_audit_revision_history.revision_date.comparator>
        <pdbx_audit_revision_history.revision_date.min>{date}</pdbx_audit_revision_history.revision_date.min>
        <pdbx_audit_revision_history.ordinal.comparator><![CDATA[>]]></pdbx_audit_revision_history.ordinal.comparator>
        <pdbx_audit_revision_history.ordinal.value>1</pdbx_audit_revision_history.ordinal.value>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>
    """

    def __init__(self, override_date=None):
        self.override_date = datetime.datetime.strptime(
            override_date, "%Y-%m-%d").date()

    @classmethod
    def config(cls, params, config):
        return config.extract(
            {"foldlib": ["override_date"]}
        ).merge_params(params)

    @classmethod
    def search(cls, date):
        """Search the RCSB for PDB entries released or revised after `date`."""
        query = cls.XML_SEARCH.format(date=date.strftime("%Y-%m-%d"))
        result = urllib.request.urlopen(
            cls.SEARCH_URL,
            data=urllib.parse.quote_plus(query).encode("UTF-8"))
        return result.read().decode("ASCII").split()

    def run(self, data, config=None, pipeline=None):
        """Retrieve new or updated PDB IDs from the RCSB."""
        template_db = data["template_db"]
        # Use configured date first, then fall back to fold library update
        # date, and then to the full list of PDBs.
        if self.override_date is not None:
            update_date = self.override_date
        else:
            update_date = template_db.updated

        if update_date is None:
            self.logger.warning(
                "Fold library has never been updated: getting all PDBs")
            pdb_ids = phyre_engine.tools.pdb.get_current()
        else:
            pdb_ids = self.search(update_date)

        data["templates"] = [{"PDB": entry} for entry in pdb_ids]
        self.logger.info("Retrieved %d PDBs updated since %s",
                         len(data["templates"]), update_date)
        return data

class CompressTemplate(Component):
    """
    Compress ``template_obj`` to take much less space when pickled.

    This class can be used to "compress" the ``template_obj`` key into a simple
    dictionary containing the metadata. This is useful for pickling purposes:
    there is no need to pickle the bulky `chain` attribute of the template when
    a PDB file containing the template has already been written.

    This component will replace the ``template_obj`` key with the
    ``template_metadata`` key. The ``template_metadata`` key is a dictionary
    containing the fields passed to the constructor of
    :py:class:`phyre_engine.tools.template.Template` except for the `chain`
    parameter.
    """
    ADDS = ["template_metadata"]
    REMOVES = ["template_obj"]
    REQUIRED = ["template_obj"]

    def run(self, data, config=None, pipeline=None):
        """Compress template object."""
        template = self.get_vals(data)
        data["template_metadata"] = {
            "mapping": template.mapping,
            "canonical_seq": template.canonical_seq,
            "canonical_indices": template.canonical_indices,
        }
        del data["template_obj"]
        return data

class UncompressTemplate(Component):
    """
    Uncompress the ``template_metadata`` field into the ``template_obj`` field.

    This component will convert the ``template_metadata`` field into a
    :py:class:`phyre_engine.tools.template.Template` object in the
    ``template_obj`` field.

    :param str chain_dir: Root directory of the template library. If this is
        `None`, the template will have its `chain` attribute set to `None`.
        This will speed up loading significantly, but cause any components
        that attempt to use the template chain to fail.
    """
    ADDS = ["template_obj"]
    REMOVES = ["template_metadata"]
    REQUIRED = ["template_metadata", "PDB", "chain"]

    def __init__(self, chain_dir):
        self.chain_dir = chain_dir

    @classmethod
    def config(cls, params, config):
        return config.extract({"foldlib": ["chain_dir"]}).merge_params(params)

    def run(self, data, config=None, pipeline=None):
        """Uncompress template object."""
        metadata, pdb_id, chain_id = self.get_vals(data)
        metadata["pdb_id"] = pdb_id
        metadata["chain_id"] = chain_id

        if self.chain_dir is not None:
            pdb_parser = Bio.PDB.PDBParser(QUIET=True)
            pdb_file = phyre_engine.tools.pdb.find_pdb(pdb_id, chain_id,
                                                       self.chain_dir)
            with phyre_engine.tools.pdb.open_pdb(pdb_file) as pdb_in:
                structure = pdb_parser.get_structure("", pdb_in)
                metadata["chain"] = structure[0]["A"]
        else:
            metadata["chain"] = None

        data["template_obj"] = Template(**metadata)
        del data["template_metadata"]
        return data


class DropIndices(Component):
    """
    Drop all indices in the fold library SQL database, in preparation for bulk
    inserts.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["template_db"]

    def run(self, data, config=None, pipeline=None):
        """Drop fold library indices."""
        data["template_db"].drop_indices()


class CreateIndices(Component):
    """
    Create all indices in the fold library SQL database, in preparation for
    fast searching.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["template_db"]

    def run(self, data, config=None, pipeline=None):
        """Create fold library indices."""
        data["template_db"].create_indices()


class FoldLibMetadata(phyre_engine.component.pdb.MMCIFMetadata):
    """
    Convenience sub-class of
    :py:class:`phyre_engine.component.pdb.MMCIFMetadata` for extracting
    the metadata required by the fold library.

    This component will also do some extra cleanup of the metadata:

    * For structures such as ``1iob`` that are determined by a combination of
      multiple methods, the first method listed is chosen to represent the
      structure.
    """

    FIELDS = {
        "deposition_date": '''date("_pdbx_database_status.recvd_initial_deposition_date", '%Y-%m-%d')''',
        "last_update_date": '''date(to_array("_pdbx_audit_revision_history.revision_date")[-1], '%Y-%m-%d')''',
        "release_date": '''date(to_array("_pdbx_audit_revision_history.revision_date")[0], '%Y-%m-%d')''',
        "method": '''to_array("_exptl.method")[0]''',
        "resolution": '''to_number("_reflns.d_resolution_high")''',
        "organism_name": '''to_array("_entity_src_gen.pdbx_host_org_scientific_name" || "_pdbx_entity_src_syn.organism_scientific")[0]''',
        "organism_id": '''to_number(to_array("_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id" || "_pdbx_entity_src_syn.ncbi_taxonomy_id")[0])''',
        "title": '''"_struct.title"''',
        "descriptor": '''"_struct.pdbx_descriptor"''',
    }

    @classmethod
    def config(cls, params, config):
        return config.extract({"foldlib": ["mmcif_dir"]}).merge_params(params)

    def __init__(self, mmcif_dir):
        super().__init__(mmcif_dir, self.FIELDS)

    def run(self, data, config=None, pipeline=None):
        """Retrieve structure metadata required by the fold library."""
        data = super().run(data, config, pipeline)
        return data


class AddPDB(Component):
    """
    Insert a PDB entry into the template database.

    New PDB entries will be added to the database, and old PDB entries will
    be updated.

    The ``metadata`` field in the pipeline state must contain the fields
    required by
    :py:meth:`phyre_engine.tools.template.TemplateDatabase.add_pdb`.

    .. note::

        This component does *not* commit any data to the database. Call
        :py:class:`.CommitChanges` to write data to the database.
    """
    REQUIRED = ["PDB", "metadata", "template_db"]
    ADDS = []
    REMOVES = []


    def run(self, data, config=None, pipeline=None):
        """Add PDB entry to template database."""
        pdb_id, metadata, template_db = self.get_vals(data)
        try:
            template_db.get_pdb(pdb_id)
            template_db.update_pdb(pdb_id, metadata)
        except template_db.PdbNotFoundException:
            template_db.add_pdb(pdb_id, metadata)
        return data


class AddTemplate(Component):
    """
    Insert a template into the template database.

    This component will first delete *all* data for each template to be added.
    This is in an effort to avoid any obsolete data remaining in the database.

    This class requires the keys ``PDB``, ``chain`` and ``template_obj`` to be
    present in the pipeline state.

    .. note::

        This component does *not* commit any data to the database. Call
        :py:class:`.CommitChanges` to write data to the database.
    """

    REQUIRED = ["template_obj", "template_db"]
    ADDS = []
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Add template to template database."""
        template, template_db = self.get_vals(data)
        template_db.del_template(template.pdb_id, template.chain_id)
        template_db.add_template(template)
        return data


class UpdateSequenceRepresentatives(Component):
    """
    Update list of sequence representatives in the fold library.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = ["template_db"]


    def run(self, data, config=None, pipeline=None):
        """Update sequence representatives."""
        template_db = self.get_vals(data)
        template_db.update_seq_reps()
        return data


class SequenceRepresentatives(Component):
    """
    Retain only those entries in the ``templates`` list that are sequence
    representatives.

    This component will get the IDs of all sequence representatives in the
    database using
    :py:meth:`phyre_engine.tools.template.TemplateDatabase.sequence_reps`.  It
    will then filter the ``templates`` list to contain only the intersection of
    the sequence representatives and current templates.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = ["template_db", "templates"]

    def run(self, data, config=None, pipeline=None):
        """Reduce templates list to sequence representatives."""
        template_db, templates = self.get_vals(data)

        # (PDB, chain) IDs for sequence representatives.
        sequence_reps = template_db.sequence_reps()

        # Mapping of (PDB, chain) to template
        template_map = {(i["PDB"].lower(), i["chain"]): i for i in templates}

        # Intersection of the two. That is, all sequences that need to have
        # HMMs rebuilt.
        rebuild = set(template_map.keys()).intersection(set(sequence_reps))
        templates = [template_map[i] for i in rebuild]

        self.logger.info("Reducing %d templates to %d representatives",
                         len(data["templates"]), len(templates))

        data["templates"] = templates
        return data


class ExpandSequenceRepresentatives(Component):
    """

    Return a list of templates with the same canonical sequence as the current
    template.

    This component will look up all chains with identical sequences to the chain
    given by the ``PDB`` and ``chain`` fields in the pipeline state. If the
    `duplicate_extra` parameter is `True`, then all fields in the current state
    are copied into the new objects. Otherwise, the objects will just contain
    the ``PDB`` and ``chain`` fields.

    :param bool duplicate_extra: If `True`, copy all fields from the
        current pipeline state into each of the new chains.

    .. warning::

        This component returns a *list* of templates, so it should only be used
        when wrapped in a :py:class:`phyre_engine.component.component.Map`
        component, or some other component that can handle a pipeline state of
        type `list`.
    """

    ADDS = []
    REMOVES = []
    REQUIRED = ["template_db", "PDB", "chain"]

    def __init__(self, duplicate_extra=True):
        self.duplicate_extra = duplicate_extra

    def run(self, data, config=None, pipeline=None):
        """Expand sequence cluster."""
        template_db, pdb_id, chain_id = self.get_vals(data)
        clus_ids = template_db.expand_seq_reps(pdb_id, chain_id)
        new_state = [data]
        for new_pdb_id, new_chain_id in clus_ids:
            new_member = copy.copy(data) if self.duplicate_extra else {}
            new_member["PDB"] = new_pdb_id
            new_member["chain_id"] = new_chain_id
            new_state.append(new_member)
        self.logger.info("%s_%s expanded to %d members",
                         pdb_id, chain_id, len(new_state))
        return new_state


class BuildProfiles(Component):
    """
    Build all profiles (A3M, HHM and CS219) for the template in the
    ``template_obj`` field of the pipeline state.

    Files will be saved in subdirectories named for the file type, middle two
    letters of the PDB code, and the template ID. That is, the profiles for
    chain ``A`` of PDB code ``1xyz`` are saved in
    :file:`{ext}/xy/1xyz_A.{ext}``, where :file:`{ext}` indicates the type of
    the file.

    Internally, this class will call :py:class:`phyre_engine.component.hhsuite.HHBlits`,
    :py:class:`phyre_engine.component.hhsuite.AddPsipred`,
    :py:class:`phyre_engine.component.hhsuite.AddDssp`
    (and :py:class:`phyre_engine.component.secstruc.DSSP` beforehand),
    :py:class:`phyre_engine.component.hhsuite.HHMake` and
    :py:class:`phyre_engine.component.hhsuite.CSTranslate`.

    :param str base_dir: Base directory for storing profile files.
    :param str blits_db: Location of the hhblits sequence database.
    :param str hh_bin_dir: Directory containing hhsuite binaries.
    :param str HHLIB: HHLIB environment variable.

    :param str dssp_bin_dir: Location of the DSSP binary.

    :param int cpu: Number of CPUs to use for hhblits.
    :param int iterations: Number of iterations to use for hhblits.
    :param int verbose: Verbosity of hhsuite components.

    :param bool overwrite: Overwrite existing files.
    """
    REQUIRED = ["template_obj", "structure"]
    ADDS = ["sequence", "a3m", "hhm", "cs219"]
    REMOVES = []

    EXTENSIONS = ("a3m", "hhm", "cs219", "fasta")

    def __init__(self, profile_dir,
                 blits_db, hh_bin_dir=None, HHLIB=None,
                 dssp_bin_dir=None,
                 hh_options=None,
                 overwrite=True):
        self.profile_dir = profile_dir
        self.blits_db = blits_db
        self.hh_bin_dir = hh_bin_dir
        self.HHLIB = HHLIB
        self.hh_options = hh_options

        self.dssp_bin_dir = dssp_bin_dir

        self.overwrite = overwrite

    @classmethod
    def config(cls, params, config):
        return config.extract({
            "foldlib": ["profile_dir", "overwrite"],
            "hhsuite": [
                ("bin_dir", "hh_bin_dir"),
                "HHLIB",
                ("options", "hh_options"),
                ("database", "blits_db"),
            ],
            "dssp": [("bin_dir", "dssp_bin_dir")],
        }).merge_params(params)

    def write_fasta(self, filename, template):
        """Write a FASTA file containing the template sequence."""
        with filename.open("w") as fas_out:
            fas_out.write(">{}_{}\n{}\n".format(
                template.pdb_id, template.chain_id, template.canonical_seq))

    def output_files(self, pdb_id, chain_id):
        """Return a dictionary of output files, indexed by type."""
        pdb_id = pdb_id.lower()
        files = {}
        for ext in self.EXTENSIONS:
            basename = "{}_{}.{}".format(pdb_id, chain_id, ext)
            file_path = Path(self.profile_dir, ext, pdb_id[1:3], basename)
            files[ext] = file_path
        return files

    def cstranslate_opts(self, a3m_file, cs219_file):
        """
        Return the options dictionary for cstranslate.

        These are the options recommended in the hhsuite user manual.
        """
        data_dir = Path(self.HHLIB, "data")
        return {
            "context-data": data_dir / "context_data.lib",
            "alphabet": data_dir / "cs219.lib",
            "pc-admix": 0.3,
            "pc-ali": 4,
            "informat": "a3m",
            "outfile": cs219_file,
        }

    def components(self, output_files, pipe_config):
        """
        List of components to be called.

        In order, the following components will be called:

        :py:class:`phyre_engine.component.hhsuite.HHBlits`
            Generate the ``a3m`` multiple sequence alignment.

        :py:class:`phyre_engine.component.hhsuite.AddPsipred`
            Calculate predicted secondary structure of the template and add it
            to the ``a3m`` file.

        :py:class:`phyre_engine.component.secstruc.DSSP`
            Calculate secondary structure of the template.

        :py:class:`phyre_engine.component.hhsuite.AddDssp`
            Add calculated secondary structure to the ``a3m`` file.

        :py:class:`phyre_engine.component.hhsuite.HHMake`
            Build an ``hhm`` file from the ``a3m`` file.

        :py:class:`phyre_engine.component.hhsuite.CSTranslate`
            Build a 219-letter alphabet representation of the sequence
            profile for fast prefiltering by hhblits.

        :returns: List of components.
        """
        # Module alias for RSI reasons
        hhsuite = phyre_engine.component.hhsuite

        # Easy function for defining a TryCatch component.
        def trycatch(components):
            child_pipe = phyre_engine.Pipeline(
                config=pipe_config,
                components=components)
            return phyre_engine.component.component.TryCatch(
                pipeline=child_pipe,
                pass_through=True)

        # Define individual components
        hhblits_opts = self.hh_options if self.hh_options is not None else {}
        hhblits_opts["oa3m"] = output_files["a3m"]
        hhblits_opts["input"] = output_files["fasta"]

        hhblits = hhsuite.HHBlits(
            database=self.blits_db,
            bin_dir=self.hh_bin_dir,
            HHLIB=self.HHLIB,
            options=hhblits_opts,
            cache_dir=output_files["a3m"].parent)
        add_psipred = hhsuite.AddPsipred(HHLIB=self.HHLIB)
        dssp = phyre_engine.component.secstruc.DSSP(bin_dir=self.dssp_bin_dir)
        add_dssp = hhsuite.AddDssp()
        hhmake = hhsuite.HHMake(
            bin_dir=self.hh_bin_dir,
            HHLIB=self.HHLIB,
            options={
                "output": output_files["hhm"],
            },
            cache_dir=output_files["hhm"].parent)
        cstranslate = hhsuite.CSTranslate(
            flags=["binary"],
            bin_dir=self.hh_bin_dir,
            HHLIB=self.HHLIB,
            options=self.cstranslate_opts(
                output_files["a3m"], output_files["cs219"]),
            cache_dir=output_files["cs219"].parent)
        return [
            hhblits,
            trycatch([add_psipred]),
            trycatch([dssp, add_dssp]),
            hhmake,
            cstranslate
        ]

    def run(self, data, config=None, pipeline=None):
        """Build fold library profile files for template."""
        template, _ = self.get_vals(data)

        output_files = self.output_files(template.pdb_id, template.chain_id)
        # Make parent directories of output files and delete existing files.
        for out_file in output_files.values():
            out_file.parent.mkdir(parents=True, exist_ok=True)
            if out_file.exists() and self.overwrite:
                out_file.unlink()

        if not output_files["fasta"].exists():
            self.write_fasta(output_files["fasta"], template)

        for component in self.components(output_files, config):
            data["sequence"] = template.canonical_seq
            data = component.run(data, config, pipeline)

        # Remove this to keep state self contained and avoid bulking it out too
        # much.
        if "secondary_structure" in data:
            del data["secondary_structure"]

        return data
