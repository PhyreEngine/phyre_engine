"""Test components in the :py:mod:`phyre_engine.component.cluster` module."""
from pathlib import Path
import tempfile
import unittest

import phyre_engine.component.cluster as cluster
import phyre_engine.test
from phyre_engine.test.data.minimal_template import MINIMAL_PDB


class TestMaxCluster(unittest.TestCase):
    """Test :py:class:`phyre_engine.component.cluster.MaxCluster`."""

    @classmethod
    def setUpClass(cls):
        """Write two minimal PDBs for alignment purposes."""
        cls.tmpdir = Path(tempfile.mkdtemp("-test", "phyreengine-"))
        cls.PIPELINE = {"templates": []}
        for i in range(2):
            model = cls.tmpdir / "{}.pdb".format(i)
            with model.open("w") as pdb_out:
                pdb_out.write(MINIMAL_PDB)
                cls.PIPELINE["templates"].append({"model": str(model)})

    @classmethod
    def tearDownClass(cls):
        """Remove temporary directory."""
        shutil.rmtree(str(cls.tmpdir))

    SAMPLE_OUTPUT = textwrap.dedent("""\
    INFO  : Reading PDB list file '/dev/fd/63'
    INFO  : Successfully read 10 / 10 PDBs from list file '/dev/fd/63'
    INFO  : Successfully read 10 Chain structures
    INFO  : Processed 40 of 45 MAXSUBs
    INFO  : CPU time = 0.06 seconds
    INFO  : ======================================
    INFO  : 3Djury (Threshold: >    20 pairs @ > 0.200)
    INFO  : ======================================
    INFO  : Rank     Model    Pairs       File
    INFO  :     1 :        5      916       04-1qob_A.pdb
    INFO  :     2 :        6      916       04-1qob_B.pdb
    INFO  :     3 :        8      915       06-1qof_A.pdb
    INFO  :     4 :        9      915       06-1qof_B.pdb
    INFO  :     5 :        4      910       03-1j7a_A.pdb
    INFO  :     6 :        7      909       05-1j7c_A.pdb
    INFO  :     7 :        2      905       02-2cjn_A.pdb
    INFO  :     8 :        3      900       02-2cjo_A.pdb
    INFO  :     9 :       10      846       07-1fxi_A.pdb
    INFO  :    10 :        1      645       02-1roe_A.pdb
    INFO  : ======================================
    INFO  : Pairwise maximum linkage clustering
    INFO  : ======================================
    INFO  : Hierarchical Tree
    INFO  : ======================================
    INFO  : Node     Item 1   Item 2      Distance
    INFO  :     0 :        7        4        0.000  05-1j7c_A.pdb  03-1j7a_A.pdb
    INFO  :    -1 :        9        6        0.001  06-1qof_B.pdb  04-1qob_B.pdb
    INFO  :    -2 :        8        5        0.001  06-1qof_A.pdb  04-1qob_A.pdb
    INFO  :    -3 :       -1       -2        0.007
    INFO  :    -4 :       -3        0        0.014
    INFO  :    -5 :        3        2        0.016  02-2cjo_A.pdb  02-2cjn_A.pdb
    INFO  :    -6 :       -4       -5        0.125
    INFO  :    -7 :       10       -6        0.153  07-1fxi_A.pdb
    INFO  :    -8 :       -7        1        0.443                 02-1roe_A.pdb
    INFO  : ======================================
    INFO  : 1 Clusters @ Threshold  0.800 (0.8)
    INFO  : ======================================
    INFO  : Item     Cluster
    INFO  :     1 :        1                        02-1roe_A.pdb
    INFO  :     2 :        1                        02-2cjn_A.pdb
    INFO  :     3 :        1                        02-2cjo_A.pdb
    INFO  :     4 :        1                        03-1j7a_A.pdb
    INFO  :     5 :        1                        04-1qob_A.pdb
    INFO  :     6 :        1                        04-1qob_B.pdb
    INFO  :     7 :        1                        05-1j7c_A.pdb
    INFO  :     8 :        1                        06-1qof_A.pdb
    INFO  :     9 :        1                        06-1qof_B.pdb
    INFO  :    10 :        1                        07-1fxi_A.pdb
    INFO  : ======================================
    INFO  : Centroids
    INFO  : ======================================
    INFO  : Cluster  Centroid  Size        Spread
    INFO  :     1 :        9       10        0.093  06-1qof_B.pdb
    INFO  : ======================================
    """)

    def test_parse(self):
        """Parse sample output."""
        self.assertEqual(
            cluster.MaxCluster.parse_maxcluster_output(self.SAMPLE_OUTPUT), {
                "04-1qob_A.pdb": 916,
                "04-1qob_B.pdb": 916,
                "06-1qof_A.pdb": 915,
                "06-1qof_B.pdb": 915,
                "03-1j7a_A.pdb": 910,
                "05-1j7c_A.pdb": 909,
                "02-2cjn_A.pdb": 905,
                "02-2cjo_A.pdb": 900,
                "07-1fxi_A.pdb": 846,
                "02-1roe_A.pdb": 645,
            })

    @phyre_engine.test.requireFields(["bin_dir"], ["tools", "maxcluster"])
    def test_run(self):
        """Run MaxCluster to calculate number of aligned residue pairs."""
        bin_dir = phyre_engine.test.config["tools"]["maxcluster"]["bin_dir"]
        maxcluster = cluster.MaxCluster(
            options={"jury_pair_threshold": 0},
            bin_dir=bin_dir)
        results = maxcluster.run(copy.deepcopy(self.PIPELINE))
        self.assertEqual(results["templates"][0]["jury_pairs"], 12)
        self.assertEqual(results["templates"][1]["jury_pairs"], 12)


@phyre_engine.test.requireFields(["bin_dir"], ["tools", "em4gmm"])
class TestEM4GMM(unittest.TestCase):
    """Test the EM4GMM component."""

    SAMPLES = [
        {"a": -5458.943, "b": -3746.670},
        {"a": 8619.275, "b": -3768.136},
        {"a": 8662.408, "b": -4027.530},
        {"a": -5575.400, "b": -3801.035},
        {"a": 8588.518, "b": -3912.427},
        {"a": 8612.947, "b": -3799.723},
        {"a": 8711.633, "b": -3954.823},
        {"a": 8689.383, "b": -3880.307},
        {"a": -5464.863, "b": -3586.280},
        {"a": -389.876, "b": 8330.014},
    ]

    # Generated on the samples above
    SAMPLE_CLASS = [
        {"lprob": [ -11576.3397417498, -16.7851236393 ], "class": 1 },
        {"lprob": [ -12945.1427112632, -15.9364833638 ], "class": 1 },
        {"lprob": [ -13439.9561536653, -17.0935529759 ], "class": 1 },
        {"lprob": [ -11703.7496472569, -16.6048579782 ], "class": 1 },
        {"lprob": [ -13195.4888249515, -16.0241773898 ], "class": 1 },
        {"lprob": [ -12999.8502011259, -15.8377932709 ], "class": 1 },
        {"lprob": [ -13326.9026297243, -16.3261393795 ], "class": 1 },
        {"lprob": [ -13180.2042333290, -15.8952609991 ], "class": 1 },
        {"lprob": [ -11288.9425377661, -18.5322030763 ], "class": 1 },
        {"lprob": [ -13.5152050419, -4914.6062961773 ], "class": 0 }
    ]

    MODEL_DETAILS = {
        "dimension": 2,
        "classes": 2,
        "minimum_dcov": [ 20855.4468057661, 6661.6023938132 ],
        "model": [
             {"class": 0, "lprior": -2.3025850930,
              "means": [-389.8760000000, 8330.0140000000],
              "dcov": [20855.4468057661, 6661.6023938132]},
             {"class": 1, "lprior": -0.1053605157,
              "means": [3931.6620000000, -3830.7701111111],
              "dcov": [44477868.2776024342, 15093.8554352075]}
        ]
    }

    def test_cluster(self):
        """Cluster some synthetic data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            mixture_model = Path(tmpdir, "model.gmm")
            model_details = Path(tmpdir, "train.details.json")
            sample_details = Path(tmpdir, "classify.details.json")
            em4gmm = cluster.EM4GMM(
                "data", "[a, b]",
                mixture_model, model_details, sample_details,
                num_components=2,
                bin_dir=phyre_engine.test.config["tools"]["em4gmm"]["bin_dir"])
            results = em4gmm.run({"data": self.SAMPLES})

        self.assertEqual(results["clusters"], self.MODEL_DETAILS)
        for sample, result in zip(results["data"], self.SAMPLE_CLASS):
            self.assertEqual(sample["class"], result["class"])
            self.assertEqual(sample["lprob"], result["lprob"])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
