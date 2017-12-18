import json
import subprocess
import tempfile

import jmespath

from phyre_engine.component.component import Component
from phyre_engine.tools.external import ExternalTool
from phyre_engine.tools.jmespath import JMESExtensions


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


class EM4GMM(Component):
    """
    Use `em4gmm <https://github.com/juandavm/em4gmm>`_for clustering via the
    Expectation Maximisation (EM) algorithm using Gaussian Mixture Models
    (GMMs).

    This component is used to cluster arbitrary lists in the pipeline state
    using GMMs. The list to be clustered is selected via the JMESPath query
    `select_expr`, and the dimensions are chosen via the JMESPath query
    `dimensions_expr`.

    For each sample, the fields ``class`` and ``lprob`` are added,
    corresponding to the fields in the log file from ``gmmclass``.

    :param str select_expr: JMESPath expresion selecting the list of samples
        to cluster.

    :param str dimensions_expr: JMESPath expression evaluated relative to each
        item in the list returned by `select_expr`, giving the values of each
        dimension.

    The remaining parameters for the component mirror the command line options
    of `gmmtrain` and `gmmclass`:


    :param str mixture_model: Name of the file used to save the trained mixture
        model (``-m`` option of ``gmmtrain``).

    :param str model_details: Log file name containing details of the model
        (``-r`` option of ``gmmtrain``).

    :param str sample_details: File name used to save details of sample
        classifications (``-r`` option of ``gmmclass``).

    :param int num_components: Optional number of components of the mixture (
        (``-n`` option of ``gmmtrain``).

    :param float merge: Optional merge threshold based on similarity (``-u``
        option of ``gmmtrain``).

    :param float stop: Optional stop criterion based on likelihood (``-s``
        option of ``gmmtrain``).

    :param int iterations: Optional maximum number of EM iterations (``-i``
        option of ``gmmtrain``).

    :param int threads: Optional maximum number of threads used (``-t`` option
        of ``gmmtrain``and ``gmmclass``). Default is 1.

    :param str world_model: Optional world model used for smoothing (``-w``
        option of ``gmmclass``).

    :param str bin_dir: Directory containing the executables ``gmmtrain`` and
        ``gmmclass``. By default, the executables are looked up on the system
        path.
    """
    GMMTRAIN = ExternalTool({
        "samples": "d",
        "mixture_model": "m",
        "model_details": "r",
        "num_components": "n",
        "merge": "u",
        "stop": "s",
        "iterations": "i",
        "threads": "t",
    })

    GMMCLASS = ExternalTool({
        "samples": "d",
        "mixture_model": "m",
        "sample_details": "r",
        "world_model": "w",
        "threads": "t",
    })

    ADDS = ["clusters"]
    REQUIRED = []
    REMOVES = []

    def __init__(self, select_expr, dimensions_expr,
                 mixture_model="model.gmm",
                 model_details="train.log.json",
                 sample_details="classify.log.json",
                 bin_dir=None, **kwargs):
        self.select_expr = select_expr
        self.dimensions_expr = dimensions_expr
        self.bin_dir = bin_dir

        self.gmmtrain_opts = {
            "mixture_model": mixture_model,
            "model_details": model_details,
        }
        for opt in self.GMMTRAIN.flag_map:
            if opt in kwargs:
                self.gmmtrain_opts[opt] = kwargs[opt]

        self.gmmclass_opts = {
            "mixture_model": mixture_model,
            "sample_details": sample_details,
        }
        for opt in self.GMMCLASS.flag_map:
            if opt in kwargs:
                self.gmmclass[opt] = kwargs[opt]

    def run(self, data, config=None, pipeline=None):
        """Run em4gmm for automatic clustering."""

        # Select sample
        jmes_opts = jmespath.Options(custom_functions=JMESExtensions(data))
        sample_list = jmespath.search(self.select_expr, data, jmes_opts)

        # Extract data points
        data_points = [
            jmespath.search(self.dimensions_expr, sample, jmes_opts)
            for sample in sample_list
        ]

        num_samples = len(data_points)
        num_dims = len(data_points[0])

        # Write samples to file
        with tempfile.NamedTemporaryFile("w") as sample_file:
            print("{} {}".format(num_dims, num_samples), file=sample_file)
            for sample in data_points:
                print(" ".join([str(i) for i in sample]), file=sample_file)
            sample_file.flush()

            # Run trainer
            gmmtrain_opts = {"samples": sample_file.name}
            gmmtrain_opts.update(self.gmmtrain_opts)
            gmmtrain = self.GMMTRAIN(
                (self.bin_dir, "gmmtrain"),
                options=gmmtrain_opts)
            self.logger.debug("Running %s", gmmtrain)
            subprocess.run(gmmtrain, check=True)

            # Run classifier
            gmmclass_opts = {"samples": sample_file.name}
            gmmclass_opts.update(self.gmmclass_opts)
            gmmclass = self.GMMCLASS(
                (self.bin_dir, "gmmclass"),
                options=gmmclass_opts)
            self.logger.debug("Running %s", gmmclass)
            subprocess.run(gmmclass, check=True)

        # Parse cluster definitions from trainer log file
        with open(self.gmmtrain_opts["model_details"], "r") as model_in:
            model = json.load(model_in)
            data["clusters"] = model

        # Parse sample data, adding to the samples
        with open(self.gmmclass_opts["sample_details"], "r") as samples_in:
            sample_details = json.load(samples_in)["samples_results"]

            for details in sample_details:
                i = details["sample"]
                sample_list[i]["class"] = details["class"]
                sample_list[i]["lprob"] = details["lprob"]
        return data
