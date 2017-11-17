"""Test components in the :py:mod:`phyre_engine.component.cluster` module."""
import copy
from pathlib import Path
import shutil
import tempfile
import textwrap
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

    @phyre_engine.test.requireFields(["maxcluster"], ["tools"])
    def test_run(self):
        """Run MaxCluster to calculate number of aligned residue pairs."""
        executable = phyre_engine.test.config["tools"]["maxcluster"]
        maxcluster = cluster.MaxCluster(
            options={"jury_pair_threshold": 0},
            executable=executable)
        results = maxcluster.run(copy.deepcopy(self.PIPELINE))
        self.assertEqual(results["templates"][0]["jury_pairs"], 12)
        self.assertEqual(results["templates"][1]["jury_pairs"], 12)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
