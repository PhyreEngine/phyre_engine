"""Test components of :py:mod:`phyre_engine.component.pbs.qsub` module."""
import unittest
from unittest.mock import patch, MagicMock
import phyre_engine.component.pbs.qsub as qsub
import textwrap

class TestFunctions(unittest.TestCase):
    """Test functions defined at the module pbs.qsub module level."""

    # Sample XML from 'qstat -x' (trimmed)
    _XML = textwrap.dedent(
        r"""<?xml version="1.0" encoding="UTF-8"?>
        <Data>
          <Job>
             <Job_Id>25289[].bamboo.bc.ic.ac.uk</Job_Id>
             <Job_Name>GenomeBuild</Job_Name>
             <job_state>R</job_state>
          </Job>
          <Job>
             <Job_Id>25290.bamboo.bc.ic.ac.uk</Job_Id>
             <Job_Name>STDIN</Job_Name>
             <job_state>R</job_state>
             <queue>medium</queue>
          </Job>
        </Data>
        """)

    # Corresponding job dictionary
    _JOBS = {
        "25289[].bamboo.bc.ic.ac.uk": "R",
        "25290.bamboo.bc.ic.ac.uk": "R"
    }

    @patch("subprocess.run", MagicMock(return_value=MagicMock(stdout=_XML)))
    def test_running_jobs(self):
        """Parse sample XML correctly."""
        jobs = qsub.running_jobs()
        self.assertDictEqual(jobs, self._JOBS)
