from . import jobscript
from . import exception
from phyre_engine.component.component import Component
import abc
import os
import math
import copy
import pickle
from pathlib import Path
import subprocess
import tempfile
import phyre_engine
import collections
import xml.etree.ElementTree as ET
import sys
import logging
from subprocess import CalledProcessError
from phyre_engine.component.pbs.jobscript import StartScript
import time

log = lambda: logging.getLogger(__name__)

r"""
Components for running sub-pipelines on PBS nodes.

The components in this module are responsible for splitting a pipeline on a
certain field, queueing a sub-pipeline on each subset of data data on separate
nodes, and gathering the results of each pipeline so serial processing can
resume.
"""

#: Sub-pipeline and data slice used by each worker.
SubPipeline = collections.namedtuple("SubPipeline", ["pipeline"])

#: Used to save data when a worker is complete
CompletedState = collections.namedtuple(
    "CompletedState", ["data"])

# File name templates
WORKER_STATE = "state-{}.pickle"
RESUME_PIPELINE = "pipeline.pickle"

class RemoteJobBase(metaclass=abc.ABCMeta):
    """
    Abstract base class representing remote jobs started by this module.
    """

    @property
    @abc.abstractmethod
    def id(self):
        """
        :return: ID of this job, as assigned by ``qsub``.
        """
        pass

    @property
    @abc.abstractmethod
    def state(self):
        """
        Retrieve pipeline state from remote job(s).

        :raise phyre_engine.component.pbs.exception.WorkerFailure: If the output
            of a worker was unexpected.
        """
        pass

class RemoteJob(RemoteJobBase):
    """
    Represents a remote job run using ``qsub``.

    The state of the job is read from a single pickle file.

    :param str job_id: Job ID, as assigned by ``qsub``.
    :param str state_file: File containing the pickled pipeline state.
    """

    def __init__(self, job_id, state_file):
        self.job_id = job_id
        self.state_file = Path(state_file)

    @property
    def id(self):
        return self.job_id

    @property
    def state(self):
        with self.state_file.open("rb") as state_in:
            data = pickle.load(state_in)
            if not isinstance(data, CompletedState):
                raise exception.UncompletedState(self.state_file)
            return data.data

class RemoteArrayJob(RemoteJobBase):
    """
    Represents a remote array job run using ``qsub`` with the ``-t`` option.

    The state of the job is read from multiple pickle files.

    :param str job_id: Job ID, as assigned by ``qsub``.
    :param list state_files: List of files containing the pickled pipeline
        state.
    :param str join_var: Join states on this variable.
    """

    def __init__(self, job_id, state_files, join_var):
        self.job_id = job_id
        self.state_files = [Path(s) for s in state_files]
        self.join_var = join_var

    @property
    def id(self):
        return self.job_id

    @property
    def state(self):
        # Read the first state, and then append the contents of each sliced
        # array from the remaining states.

        full_state = None
        for state_file in self.state_files:
            with state_file.open("rb") as state_in:
                partial_state = pickle.load(state_in)
                if not isinstance(partial_state, CompletedState):
                    raise exception.UncompletedState(state_file)

                if full_state is None:
                    full_state = partial_state.data
                else:
                    full_state[self.join_var].extend(
                        partial_state.data[self.join_var])
        return full_state

    def __repr__(self):
        return "<{name}(job_id={job_id}, join_var={join_var}>".format(
            name=type(self).__name__,
            job_id=repr(self.job_id),
            join_var=repr(self.join_var))

class BaseQsub(Component):
    def __init__(
            self, storage_dir, worker_pipeline,
            name, worker_config=None, qsub_args=None):
        self.storage_dir = Path(storage_dir)

        # Load pipeline if it isn't already a Pipeline object
        if not isinstance(worker_pipeline, phyre_engine.Pipeline):
            worker_pipeline = phyre_engine.Pipeline.load(worker_pipeline)
        self.worker_pipeline = worker_pipeline

        self.name = name
        self.worker_config = worker_config if worker_config is not None else {}
        self.qsub_args = qsub_args if qsub_args is not None else []

    def _save_pipeline(self, pickle_path, data, master_config):
        """
        Save the current pipeline to a pickle named after the node that we are
        going to be contacting.
        """

        # Overwrite any sections of our config with
        config = copy.copy(master_config)
        config.update(self.worker_config)

        with pickle_path.open("wb") as state_out:
            pipe = copy.copy(self.worker_pipeline)
            pipe.config = config
            pipe.start = data
            pickle.dump(SubPipeline(pipe), state_out)

    def _qsub(self, script, stdout_dir, stderr_dir=None, extra=None):

        # By default, place stderr files in same directory as stdout
        if stderr_dir is None:
            stderr_dir = stdout_dir

        if extra is None:
            extra = []

        start_qsub_cmd = [
            "qsub",
            "-o", str(stdout_dir),
            "-e", str(stderr_dir),
            "-d", os.getcwd(),
            "-N", self.name,
        ] + list(self.qsub_args) + list(extra)

        popen_args = {
            "stdout": subprocess.PIPE,
            "stdin": subprocess.PIPE,
            "universal_newlines": True
        }

        with subprocess.Popen(start_qsub_cmd, **popen_args) as process:
            stdout, _ = process.communicate(str(script))
        if process.returncode != 0:
            raise CalledProcessError(
                process.returncode, start_qsub_cmd, stdout, None)

        # Return job ID
        return stdout.strip()

class Slice(BaseQsub):
    r"""
    Run a pipeline in parallel on separate nodes.
    """
    REQUIRED = []
    ADDS = ["qsub_jobs"]
    REMOVES = []

    CONFIG_SECTION = "qsub"

    # We write file names formatted like this
    _STATE_WRITE_NAME = "state-{slice_index}.pickle"

    # The jobscript reads filenames formatted like this. In this case, the
    # braces are interpreted by the shell in which the jobscript is running, not
    # by python's str.format method.
    _STATE_READ_NAME = "state-${PBS_ARRAYID}.pickle"

    def __init__(self, split_var, join_var, max_jobs,
                 *qsub_args, **qsub_kwargs):

        super().__init__(*qsub_args, **qsub_kwargs)
        self.split_var = split_var
        self.join_var = join_var
        self.max_jobs = max_jobs

    def _slice_state(self, storage_dir, data, config):
        """
        Slice pipeline and save into a series of pickles containing the subset
        of data on which each remote job will operate. Each pickle contains a
        named tuple containing the elements "pipeline" and "data".

        :return: List of paths pointing to saved states.
        """

        data_len = len(data[self.split_var])
        # Number of array elements. If we have fewer tasks to run than
        # self.max_jobs, just run them all. Otherwise, we will run self.max_jobs
        # and process some tasks in serial.
        num_jobs = min(self.max_jobs, data_len)
        elems_per_job = math.ceil(data_len / num_jobs)

        saved_states = []
        for i in range(0, num_jobs):
            pickle_name = self._STATE_WRITE_NAME.format(slice_index=i)
            pickle_path = storage_dir / pickle_name

            start = i * elems_per_job
            end = (i + 1) * elems_per_job

            data_slice = data[self.split_var][start:end]
            state = copy.copy(data)
            state[self.split_var] = data_slice

            super()._save_pipeline(pickle_path, state, config)
            saved_states.append(pickle_path)
        return saved_states

    def run(self, data, config, pipeline):
        """
        Start an array job (``-t``) with qsub.
        """

        # Get a unique directory underneath the storage directory.
        job_dir = Path(tempfile.mkdtemp("slice", "qsub", str(self.storage_dir)))

        stdout_dir = job_dir / "stdout"
        stderr_dir = job_dir / "stderr"
        stdout_dir.mkdir()
        stderr_dir.mkdir()

        pickles = self._slice_state(job_dir, data, config)

        script = jobscript.StartScript(str(job_dir / self._STATE_READ_NAME))
        qsub_args = ["-t", "0-{}".format(len(pickles) - 1)]

        job_id = self._qsub(script, stdout_dir, stderr_dir, qsub_args)

        if "qsub_jobs" not in data:
            data["qsub_jobs"] = []

        data["qsub_jobs"].append(RemoteArrayJob(job_id, pickles, self.join_var))
        return data


class ContactNodes(BaseQsub):
    """
    Run a pipeline on several nodes in parallel.

    This component operates similarly to the :py:class:`.Qsub` component, but
    rather than slicing pipeline data and running multiple copies of the same
    pipeline on arbitrary nodes assigned by the queueing system, it runs
    multiple copies of the same pipeline with the *same* data separately on
    every node it can find.

    :param list nodes: List of node hostnames. By default, this component will
        contact all nodes found by the ``pbsnodes`` command.
    """

    ADDS = ["qsub_jobs"]
    REQUIRED = []
    REMOVES = []

    CONFIG_SECTION = "qsub"

    # Name of the pickle that we write and that is read on a worker node
    _PICKLE_NAME = "state-{node}.pickle"

    def __init__(self, nodes=None, *qsub_args, **qsub_kwargs):
        super().__init__(*qsub_args, **qsub_kwargs)
        if nodes is None:
            nodes = self._online_nodes()
        self.nodes = nodes

    def run(self, data, config=None, pipeline=None):
        # Get a unique directory in which to store state, output, etc
        job_dir = Path(tempfile.mkdtemp("contact", "qsub",
                                        str(self.storage_dir)))
        stdout_dir = job_dir / "stdout"
        stderr_dir = job_dir / "stderr"
        stdout_dir.mkdir()
        stderr_dir.mkdir()

        if "qsub_jobs" not in data:
            data["qsub_jobs"] = []

        for node in self.nodes:
            pickle_path = job_dir / self._PICKLE_NAME.format(node=node)
            self._save_pipeline(pickle_path, data, config)
            script = StartScript(pickle_path)

            qsub_args = ["-lnodes={}".format(node)]
            job_id = self._qsub(script, stdout_dir, stderr_dir, qsub_args)
            data["qsub_jobs"].append(RemoteJob(job_id, pickle_path))

        return data

    @staticmethod
    def _online_nodes():
        """
        Find online nodes using pbsnodes. Get XML output from pbsnodes and find
        all those nodes that are not either down or offline.
        """
        pbs_result = subprocess.run(
            ["pbsnodes", "-x"],
            check=True, stdout=subprocess.PIPE)
        xml_tree = ET.fromstring(pbs_result.stdout.decode(sys.stdout.encoding))
        nodes = []
        for node in xml_tree.findall("./Node"):
            name = node.find("name").text
            states = set(node.find("state").text.split(","))
            if "offline" not in states and "down" not in states:
                nodes.append(name)
        return nodes

class Wait(Component):
    """
    Wait for all previously-started qsub jobs to end.

    This component will wait until all jobs in the ``qsub_jobs`` list (in the
    pipeline state) are complete.    Jobs are considered complete when they are
    in the "C" state.

    .. node::

        This component is fairly stupid: it simply polls ``qstat`` until each
        job is complete.

    :param int poll_interval: Interval in seconds between calling ``qstat``.
    """

    REQUIRED = ["qsub_jobs"]
    ADDS = []
    REMOVES = []

    def __init__(self, poll_interval=30):
        self.poll_interval = poll_interval

    def run(self, data, config=None, pipeline=None):
        self._poll(data["qsub_jobs"])
        return data

    def _poll(self, jobs):
        # Get the job IDs that we are looking for
        looking_for = set([j.id for j in jobs])

        while True:
            # Read all jobs from qstat
            qs_result = subprocess.run(
                ["qstat", "-x"],
                check=True, stdout=subprocess.PIPE)
            stdout = qs_result.stdout.decode(sys.stdout.encoding)
            xml_tree = ET.fromstring(stdout)

            # We will stop polling if only have state "C" in here:
            job_states = set()
            for job in xml_tree.findall("./Job"):
                job_id = job.find("Job_Id").text
                job_state = job.find("job_state").text

                if job_id in looking_for:
                    job_states.add(job_state)

            if job_states == set(["C"]):
                break
            time.sleep(self.poll_interval)

class LoadState(Component):
    """
    Replace the state of the current pipeline with the (combined) state of any
    jobs previously launched via the :py:mod:`phyre_engine.component.pbs.qsub`
    module.

    The pipeline state will be replaced with the state of any remote jobs
    present in the ``qsub_jobs`` element of the pipeline state. If multiple jobs
    are present in the ``qsub_jobs`` list, then the parameter ``join_var`` must
    be supplied: it is assumed that all pipeline data except for ``join_var``
    are constant, and the field ``join_var`` is combined from each job. Note
    that the order in which jobs are combined is arbitrary. An array job counts
    as a single job.

    If the ``convert_to_list`` parameter is set to ``True``, then the value
    pointed to by ``join_var`` will be coerced into a list if it is not already.

    :param str join_var: If we are joining multiple jobs (*not* just a single
        array job), this must be defined. Pipeline state is merged on this key.

    :param bool convert_to_list: Convert ``join_var`` to list. This allows the
        jobs to be joined on scalar variables.
    """

    REQUIRED = ["qsub_jobs"]
    REMOVES = ["qsub_jobs"]
    ADDS = []

    def __init__(self, join_var=None, convert_to_list=False):
        self.join_var = join_var
        self.convert_to_list = convert_to_list

    def run(self, data, config=None, pipeline=None):
        if len(data["qsub_jobs"]) > 1 and self.join_var is None:
            raise exception.TooManyjobs(data["qsub_jobs"])

        # Read states even if replace_state is false so that exceptions are
        # thrown if any job is incomplete.
        full_state = None
        for job in data["qsub_jobs"]:
            if full_state is None:
                full_state = job.state
                if (self.convert_to_list
                        and not isinstance(full_state[self.join_var], list)):
                    full_state[self.join_var] = [full_state[self.join_var]]
            else:
                job_slice = job.state[self.join_var]
                if self.convert_to_list and not isinstance(job_slice, list):
                    job_slice = [job_slice]
                full_state[self.join_var].extend(job_slice)

        del data["qsub_jobs"]
        return full_state

class Detach(BaseQsub):
    """
    Run the specified pipeline on a remote node.

    This component will execute its ``worker_pipeline`` on a remote node as a
    job started using ``qsub``.

    If the ``depends`` parameter is any non-false value, then this job will be
    started with a dependency on any jobs in the ``qsub_jobs`` list in the
    pipeline state. The ``depends`` parameter should be a tuple of length 2.
    Each element is a verb passed to ``qsub`` as part of the ``-W depends=``
    parameter. The first element is used for non-array jobs, and the second is
    used for array jobs.

    For example, the tuple ``("afterok", "afterokarray")`` would result in a
    command line for ``qsub`` like
    ``-W depends=afterok:1234:789,afterokarray=101112[]``.

    :param tuple depends: Verbs to use to wait for previous jobs to complete.
        Set this to ``False`` to prevent any waiting.

    .. seealso::

        :py:class:`phyre_engine.component.pbs.qsub.BaseQsub`:
            Extra command line parameters.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    CONFIG_SECTION = "qsub"

    _PICKLE_NAME = "detached.pickle"

    def __init__(self, depends=("afterok", "afterokarray"),
                 *qsub_args, **qsub_kwargs):
        super().__init__(*qsub_args, **qsub_kwargs)
        self.depends = depends


    def _depend_args(self, data):
        """
        Build the command line argument specifying job dependencies.

        This has two sections separated by a comma, one for single jobs and one
        for array jobs. Each job ID in each section is separated by a colon. It
        is valid to have zero of either kind of job.
        """

        if self.depends and "qsub_jobs" in data and data["qsub_jobs"]:
            # Separate single and array jobs based on whether they contain the
            # string "[]" in their ID.
            single_jobs = []
            array_jobs = []
            for job in data["qsub_jobs"]:
                if "[]" in job.id:
                    array_jobs.append(job.id)
                else:
                    single_jobs.append(job.id)

            single_job_str = ""
            array_job_str = ""
            if single_jobs:
                single_job_str = ":".join([self.depends[0]] + single_jobs)
            if array_jobs:
                array_job_str = ":".join([self.depends[1]] + array_jobs)
            job_strs = [s for s in [single_job_str, array_job_str] if s]
            depend_str = "depend=" + ",".join(job_strs)
            return ["-W", depend_str]
        else:
            # No extra arguments if we are not using dependencies or there are
            # no jobs on which this job depends.
            return []

    def run(self, data, config=None, pipeline=None):
        # Store pipeline state in a unique directory
        job_dir = Path(
            tempfile.mkdtemp("detach", "qsub", str(self.storage_dir)))

        pickle_path = job_dir / self._PICKLE_NAME
        stdout_dir = job_dir / "stdout"
        stderr_dir = job_dir / "stderr"
        stdout_dir.mkdir()
        stderr_dir.mkdir()

        self._save_pipeline(pickle_path, data, config)
        script = jobscript.StartScript(pickle_path)
        qsub_args = self._depend_args(data)
        job_id = self._qsub(script, stdout_dir, stderr_dir, qsub_args)
        job = RemoteJob(job_id, pickle_path)

        # Remove existing job_ids if they had to complete before this job could
        # start.
        if self.depends or "qsub_jobs" not in data:
            data["qsub_jobs"] = []
        data["qsub_jobs"].append(job)
        return data
