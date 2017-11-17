import subprocess
import tempfile

from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool


class MaxCluster(Component):
    """
    Use `MaxCluster <http://www.sbg.bio.ic.ac.uk/maxcluster/index.html>`_ to
    calculated the number of aligned residue pairs between each structure in
    a list of models.

    :param str executable: Location of the maxcluster executable.
    """
    REQUIRED = ["templates"]
    ADDS = []
    REMOVES = []

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

    def __init__(self, flags=None, options=None, executable="maxcluster"):
        self.executable = executable
        self.flags = [] if flags is None else flags
        self.options = {} if options is None else options

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
        # Write a list containing all models. Also create a dictionary mapping
        # models to templates so we can look up templates from the maxcluster
        # output.
        models = {}
        with tempfile.NamedTemporaryFile("w") as model_list:
            for template in data["templates"]:
                print(template["model"], file=model_list)
                models[template["model"]] = template
            model_list.flush()

            # Run maxcluster
            options = {"list": model_list.name}
            options.update(self.options)

            flags = []
            flags.extend(self.flags)

            command_line = self.MAXCLUSTER(
                (None, self.executable),
                flags=flags, options=options)
            process = subprocess.run(
                command_line,
                stdout=subprocess.PIPE, check=True, universal_newlines=True)

        # Parse results, getting a map of model names to results.
        results = self.parse_maxcluster_output(process.stdout)
        for model, pairs in results.items():
            models[model]["jury_pairs"] = pairs
        return data
