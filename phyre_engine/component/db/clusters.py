"""
This module contains components that make use of sequence clusters when building
a fold library. Sequence clusters are stored in the ``clusters`` key of the
pipeline.

.. seealso::

    :py:class:`.ClusterParser`

        For details on how clusters are stored in the ``clusters`` list.

"""
from phyre_engine.component.component import Component
import urllib.request

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
    ADDS = ['cluster_file']
    REMOVES = []

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

    def run(self, data, config=None, pipeline=None):
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
    ADDS = ['clusters']
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
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
    ADDS = ['templates']
    REMOVES = []

    def run(self, data, config=None, pipeline=None):
        """Extract representatives."""

        clusters = self.get_vals(data)
        representatives = [clus[0] for clus in clusters]
        if "templates" not in data:
            data["templates"] = []

        for rep in representatives:
            (pdb, chain) = rep.split("_")
            data["templates"].append({"PDB":pdb, "chain":chain})
        return data
