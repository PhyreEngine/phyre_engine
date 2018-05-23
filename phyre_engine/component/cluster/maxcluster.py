"""
Compare structures in various ways using `MaxCluster
<http://www.sbg.bio.ic.ac.uk/maxcluster/index.html>`_.
"""
from pathlib import Path
import re
import subprocess
import tempfile

from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool

MAXCLUSTER = ExternalTool({
    "experiment": "e",
    "prediction": "p",
    "list": "l",
    "loglevel": "L",
    "log": "log",
    "distance_cutoff": "d",
    "tm_normalisation": "N",
    "maxsub_iterations": "i",
    "independent": "in",
    "cluster_method": "C",
    "cluster_threshold_init": "T",
    "cluster_threshold_max": "Tm",
    "cluster_threshold_adjust": "a",
    "cluster_size_init": "is",
    "cluster_size_min": "ms",
    "jury_score_threshold": "s",
    "jury_pair_threshold": "P",
}, long_prefix="-")

class MaxclusterComponent(Component):
    """
    Base class for components calling MaxCluster.

    :param list[str] flags: List of flags (i.e. boolean options) to pass to
        MaxCluster.

    :param dict options: List of options (i.e. flags with values) to pass to
        MaxCluster.

    :param str executable: Location of the maxcluster executable.
    """

    def __init__(self, flags=None, options=None, bin_dir=None, cache=None):
        self.bin_dir = bin_dir
        self.flags = [] if flags is None else flags
        self.options = {} if options is None else options
        self.cache = cache

    @staticmethod
    def write_model_list(model_list, templates):
        """
        Write model names to a list for use with the ``-l`` option.
        """
        models = {}
        for template in templates:
            print(template["model"], file=model_list)
        model_list.flush()

    def run_maxcluster(self, options=None, flags=None):
        """
        Run maxcluster, merging extra options and flags with the instance
        variables `options` and `flags`.

        :returns: Result of the maxcluster process.
        :rtype: :py:class:`subprocess.CompletedProcess`
        """

        if options is None:
            options = {}
        if flags is None:
            flags = []

        # Run maxcluster
        options.update(self.options)
        flags.extend(self.flags)

        command_line = MAXCLUSTER(
            (self.bin_dir, "maxcluster"),
            flags=flags, options=options)
        process = subprocess.run(
            command_line,
            stdout=subprocess.PIPE, check=True, universal_newlines=True)
        return process

    def read_cache(self):
        """
        Returns the value of the cached output, or `None` if the cache has not
        been initialised or is disabled.
        """
        if self.cache is not None:
            if Path(self.cache).exists():
                return Path(self.cache).open("r").read()
        return None

    def write_cache(self, text):
        """Dump `text` to the cache file if `cache` is not `None`."""
        if self.cache is not None:
            with open(self.cache, "w") as cache_out:
                cache_out.write(text)


class Jury(MaxclusterComponent):
    """
    Use `MaxCluster <http://www.sbg.bio.ic.ac.uk/maxcluster/index.html>`_ to
    calculated the number of aligned residue pairs between each structure in
    a list of models (i.e. use the 3D-JURY algorithm).

    This operation scales quadratically with the number of models to be
    aligned, and so can be quite expensive. Because of this, the default
    behaviour is to save the output of ``maxcluster`` to a cache file, and
    re-use the cache when possible. To disable this behaviour, set `cache =
    None`.

    :param str cache: Save maxcluster output in this file. Set to `None`
        to disable caching.

    .. seealso::

        :py:class:`.MaxclusterComponent`
            For constructor parameters.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("cache", "maxcluster.pairwise")
        super().__init__(*args, **kwargs)

    @staticmethod
    def parse_maxcluster_output(stdout):
        """
        Parse maxcluster output.

        :return: Dictionary indexed by model name giving the number of aligned
            residue pairs between that model and all other models.
        """
        DELIMITER = "INFO  : ======================================\n"
        sections = stdout.split(DELIMITER)
        for i, section in enumerate(sections):
            if section.startswith("INFO  : 3Djury"):
                jury_section = sections[i + 1]
                break

        results = {}
        for line in jury_section.split("\n"):
            if line.startswith("INFO  : Rank") or not line.strip():
                continue
            pairs, file = line.split()[-2:]
            results[file] = int(pairs)
        return results

    def run(self, data, config=None, pipeline=None):
        """Calculate 3D-JURY results using maxcluster."""
        # Write a list containing all models. Also create a dictionary mapping
        # models to templates so we can look up templates from the maxcluster
        # output.
        maxcluster_out = self.read_cache()
        models = {t["model"]: t for t in data["templates"]}
        if maxcluster_out is None:
            with tempfile.NamedTemporaryFile("w") as model_list:
                self.write_model_list(model_list, data["templates"])
                proc = self.run_maxcluster(options={
                    "list": model_list.name,
                    # Not necessary, but allows us to share cache files with
                    # the Pairwise component.
                    "loglevel": 3,
                })
                maxcluster_out = proc.stdout
                self.write_cache(maxcluster_out)

        # Parse results, getting a map of model names to results.
        results = self.parse_maxcluster_output(maxcluster_out)
        for model, pairs in results.items():
            models[model]["jury_pairs"] = pairs
        return data


class Pairwise(MaxclusterComponent):
    """
    Use `MaxCluster <http://www.sbg.bio.ic.ac.uk/maxcluster/index.html>`_ to
    calculate pairwise similarity between all models in the ``templates``
    list.

    This component will add the field ``pairwise_similarity`` to the pipeline
    state. This field will be a dictionary indexed by tuples of model names.
    Each pair of structures will be inserted twice, so that pairs may be looked
    up in either order.

    Each element in the ``pairwise_similarity`` dictionary will either be
    `None`, indicating that the two structures could not be superimposed, or a
    dictionary with the following elements:

    ``RMSD``
        Root Mean Square Deviation of the structures, in Å.

    ``pairs``
        Number of aligned residue pairs in the MaxSub.

    `length``
        Number of residue pairs in the alignment.

    ``maxsub``
        MaxSub score.

    ``gRMSD``
        Global RMSD using the MaxSub superposition.

    ``TM``
        TM (Template Modelling) score between the two models.

    .. seealso::

        :py:class:`.MaxclusterComponent`
            For constructor parameters.

        :py:class:`.Jury`
            For details of caching.
    """
    REQUIRED = ["templates"]
    ADDS = ["pairwise_similarity"]
    REMOVES = []

    FIELD_TYPES = {
        "RMSD": float,
        "pairs": int,
        "maxsub": float,
        "length": int,
        "gRMSD": float,
        "TM": float,
    }

    RECORD_REGEX = (
        r"Iter \d+: Pairs=\s*(?P<pairs>\d+), "
        r"RMSD=\s*(?P<RMSD>\d*\.\d+), "
        r"MAXSUB=(?P<maxsub>\d*\.\d+)\. "
        r"Len=\s*(?P<length>\d+)\. "
        r"gRMSD=\s*(?P<gRMSD>\d*\.\d+), "
        r"TM=\s*(?P<TM>\d*\.\d+)$"
    )

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("cache", "maxcluster.pairwise")
        super().__init__(*args, **kwargs)

    @classmethod
    def parse_maxcluster_output(cls, text):
        """
        Parse output of maxcluster. This parses the detailed information given
        by log level 3, because we want to read the TM-score.
        """
        alignments = {}
        for line in text.split("\n"):
            if " vs. " in line and not line.startswith("INFO "):
                pdb1, pdb2 = line.split(" vs. ")
                alignments[(pdb1, pdb2)] = None
                alignments[(pdb2, pdb1)] = None
            elif line.startswith("Iter "):
                match = re.search(cls.RECORD_REGEX, line)
                if match:
                    fields = {k: cls.FIELD_TYPES[k](v)
                              for k, v in match.groupdict().items()}
                    alignments[(pdb1, pdb2)] = fields
                    alignments[(pdb2, pdb1)] = fields
        return alignments

    def run(self, data, config=None, pipeline=None):
        """Generate pairwise alignments using maxcluster."""
        maxcluster_out = self.read_cache()
        if maxcluster_out is None:
            with tempfile.NamedTemporaryFile("w") as model_list:
                self.write_model_list(model_list, data["templates"])
                proc = self.run_maxcluster(options={
                    "list": model_list.name,
                    "loglevel": 3,
                })
                maxcluster_out = proc.stdout
                self.write_cache(maxcluster_out)

        results = self.parse_maxcluster_output(maxcluster_out)
        data["pairwise_similarity"] = results
        return data
