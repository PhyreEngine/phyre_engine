"""Test components in the :py:mod:`phyre_engine.component.cluster` module."""
from pathlib import Path
import tempfile
import unittest

import phyre_engine.component.cluster as cluster
import phyre_engine.test


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
