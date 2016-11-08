from phyre_engine.component import Component
import urllib.request


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
