"""
Compare structures in various ways using MaxCluster
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
