"""Test components of :py:mod:`phyre_engine.component.pbs.qsub` module."""
import unittest
from unittest.mock import patch, MagicMock, sentinel, call
import phyre_engine.component.pbs.qsub as qsub
import copy
import itertools
import types
import textwrap

class TestTrickle(unittest.TestCase):
    """Test the Trickle class."""

    _SAMPLE_STATE = {
        "foo": "bar",
        "elements": [{"val": i} for i in range(0, 10)]
    }

    def setUp(self):
        """Create a copy of the sample state we can modify freely."""
        self.state = copy.deepcopy(self._SAMPLE_STATE)

    @staticmethod
    def job_id_gen(upper=None, state_path="/state/path"):
        """
        Return a generator of job IDs increasing incrementally. If "upper" is
        supplied, the generator terminates at that upper bound. Otherwise, it
        will continue to infinity.
        """
        base_gen = itertools.count() if upper is None else range(0, upper)
        return (qsub.RemoteJob(str(i), state_path) for i in base_gen)


    def patched_trickle(self, max_jobs, num_elems, state_path="/state/path"):
        """
        Patch all methods with side effects.

        The job IDs generated by this patched object will by integers starting
        from zero and proceeding to infinity (stored as strings).

        The number of enqueued jobs will be the number of times "enqueue" was
        called.

        We set the queue interval to "sentinel.queue_interval" and poll interval
        to sentinel.poll_interval, so we can check the arguments of the calls to
        the "time.sleep" mock to see why it was called.
        """
        def _num_running_jobs(trickle_obj, _our_jobs):
            return trickle_obj.enqueue.call_count

        trickle = qsub.Trickle("elements", max_jobs, num_elems,
                               queue_interval=sentinel.queue_interval,
                               poll_interval=sentinel.poll_interval,
                               pipeline={}, name="test")
        job_id_gen = self.job_id_gen(state_path=state_path)
        trickle.enqueue = MagicMock(side_effect=job_id_gen)

        trickle.num_running_jobs = types.MethodType(_num_running_jobs, trickle)

        trickle._save_pipeline = MagicMock()
        trickle._save_state = MagicMock()
        return trickle

    @patch("time.sleep")
    def test_enqueues_all_jobs(self, sleep_mock):
        """
        All elements are enqueued when enqueuing fewer elements than max_jobs.
        """
        # Max_jobs = 100, enqueuing 1 at a time (i.e. 10 enqueued)
        trickle = self.patched_trickle(100, 1)

        num_elems = len(self.state["elements"])
        results = trickle.run(self.state)
        # Check that all jobs were enqueued
        self.assertListEqual(
            [j.id for j in results["qsub_jobs"]],
            [j.id for j in self.job_id_gen(num_elems)])

        # Check that sleep was called the 9 times, between each queue
        self.assertEqual(
            sleep_mock.call_args_list,
            [call(sentinel.queue_interval)] * (num_elems - 1))

    @patch("time.sleep")
    def test_waits_when_full(self, sleep_mock):
        """Stops trickling jobs when the queue is full."""
        # Max_jobs = 1, enqueuing 5 at a time (i.e. 2 enqueue)
        trickle = self.patched_trickle(1, 5)

        # When first called, tell Trickle we have zero jobs enqueued. After
        # that, tell it there is one enqueued twice, then drop back to zero.
        # This should cause the first element to be enqueued (when
        # num_running_jobs() = 0), then pause for two poll intervals
        trickle.num_running_jobs = MagicMock(side_effect=[0, 1, 1, 0, 1])

        trickle.run(self.state)

        # Should have had to poll twice (when num_running_jobs() == 1)
        poll_sleeps = len([i for i in sleep_mock.call_args_list
                           if i == call(sentinel.poll_interval)])
        self.assertEqual(poll_sleeps, 2)


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
