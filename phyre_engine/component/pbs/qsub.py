from . import jobscript
from . import exception
from phyre_engine.component.component import Component, PipelineComponent
import abc
import os
import math
import copy
import pickle
from pathlib import Path
import subprocess
import tempfile
import collections
import xml.etree.ElementTree as ET
import sys
from subprocess import CalledProcessError
from phyre_engine.component.pbs.jobscript import StartScript
import time
import phyre_engine.tools.yaml as yaml


r"""
Components for running sub-pipelines on PBS nodes.

The components in this module are responsible for splitting a pipeline on a
certain field, queueing a sub-pipeline on each subset of data data on separate
nodes, and gathering the results of each pipeline so serial processing can
resume.
"""

# File name templates
WORKER_STATE = "state-{}.pickle"
RESUME_PIPELINE = "pipeline.pickle"

class RemoteJobBase(metaclass=abc.ABCMeta):
    """
    Abstract base class representing remote jobs started by this module.

    See :py:class:`.RemoteJob` for more information.
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

    The state of the job is read from a single pickle file. If the state does
    not contain the key ``qsub_complete``, which is set by
    :py:mod:`phyre_engine.component.pbs.run` when a remote pipeline is executed,
    a :py:exc:`.exception.UncompletedState` exception is raised.

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
            if "qsub_complete" not in data:
                raise exception.UncompletedState(self.state_file)
            del data["qsub_complete"]
            return data

    def __repr__(self):
        return "<{name}(job_id={job_id}, state_file={state_file}>".format(
            name=type(self).__name__,
            job_id=self.job_id, state_file=self.state_file)

class RemoteArrayJob(RemoteJobBase):
    """
    Represents a remote array job run using ``qsub`` with the ``-t`` option.

    The state of the job is read from multiple pickle files.

    See :py:class:`.RemoteJob` for more information.

    :param str job_id: Job ID, as assigned by ``qsub``.
    :param list state_files: List of files containing the pickled pipeline
        state.
    :param str join_var: Join states on this variable.
    """

    def __init__(self, job_id, state_files, join_var):
        self.job_id = job_id
        self.state_files = [Path(s) for s in state_files]
        self.join_var = join_var

    def _state_generator(self):
        for state_file in self.state_files:
            with state_file.open("rb") as state_in:
                state = pickle.load(state_in)
                if "qsub_complete" not in data:
                    raise exception.UncompletedState(self.state_file)
                del data["qsub_complete"]
                yield state

    @property
    def id(self):
        return self.job_id

    @property
    def state(self):
        return _load_state(self._state_generator(), self.join_var)

    def __repr__(self):
        return "<{name}(job_id={job_id}, join_var={join_var}>".format(
            name=type(self).__name__,
            job_id=repr(self.job_id),
            join_var=repr(self.join_var))

class BaseQsub(PipelineComponent):
    """
    Base class for components making use of the :manpage:`qsub(1)` command.

    This class contains tools for serialising a pipeline.

    :param str storage_dir: Storage directory in which to save pipeline state
        and standard output/error files.
    :param dict pipeline: Pipeline specification as expected by
        :py:meth:`phyre_engine.pipeline.Pipeline.load`.
    :param str name: Name of the job.
    :param list[str] qsub_args: Extra arguments to pass to qsub.
    """

    def __init__(self, name, *pipeline_args, storage_dir=".",
                 qsub_args=None, **pipeline_kwargs):
        super().__init__(*pipeline_args, **pipeline_kwargs)
        self.storage_dir = Path(storage_dir)
        self.name = name
        self.qsub_args = qsub_args if qsub_args is not None else []

    def _save_pipeline(self, pipe_path, base_config):
        """
        Save pipeline specification to a YAML file. If `base_config` is
        supplied, it is used as the "base" configuration: any configuration
        set in "pipeline.config" overwrites the base configuration.
        """
        pipeline_defn = self.pipeline_definition(base_config)
        with pipe_path.open("w") as pipe_out:
            yaml.dump(pipeline_defn, pipe_out)

    def _save_state(self, state_path, state):
        """Pickle pipeline state and save it to state_path."""
        with state_path.open("wb") as state_out:
            pickle.dump(state, state_out)

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

    :param str split_var: Name of a list to distribute across nodes.
    :param str join_var: Name of the list on which to merge jobs.
    :param int max_jobs: Maximum number of parallel jobs to run. Data will be
        grouped into a maximum of this number of approximately equal chunks.

    .. seealso::

        :py:class:`.BaseQsub`
            For arguments relating to qsub.

    .. note::

        If there are no items present in `split_var` a :py:exc:`ValueError` will
        be raised when this component is run.
    """
    ADDS = ["qsub_jobs"]
    REMOVES = []

    @property
    def REQUIRED(self):
        return [self.split_var]

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

            super()._save_state(pickle_path, state)
            saved_states.append(pickle_path)
        return saved_states

    def run(self, data, config=None, pipeline=None):
        """
        Start an array job (``-t``) with qsub.
        """
        if len(data[self.split_var]) == 0:
            raise ValueError(
                "Cannot split on list '{}' containing no items".format(
                    self.split_var))

        # Get a unique directory underneath the storage directory.
        job_dir = Path(tempfile.mkdtemp("slice", "qsub", str(self.storage_dir)))

        stdout_dir = job_dir / "stdout"
        stderr_dir = job_dir / "stderr"
        stdout_dir.mkdir()
        stderr_dir.mkdir()

        pipeline_path = job_dir / "pipeline.yml"
        self._save_pipeline(pipeline_path, config)
        pickles = self._slice_state(job_dir, data, config)

        script = jobscript.StartScript(
            str(pipeline_path),
            str(job_dir / self._STATE_READ_NAME))
        qsub_args = ["-t", "0-{}".format(len(pickles) - 1)]

        job_id = self._qsub(script, stdout_dir, stderr_dir, qsub_args)

        if "qsub_jobs" not in data:
            data["qsub_jobs"] = []

        data["qsub_jobs"].append(RemoteArrayJob(job_id, pickles, self.join_var))
        return data

class Trickle(BaseQsub):
    """
    Slowly trickle jobs into a queueing system using :command:`qsub`.

    This component slowly enqueues jobs to a queuing system. It will leave
    `queue_interval` seconds between jobs, helping to avoid overloading the
    scheduler. Once `max_jobs` jobs have been enqueued, this component will
    begin polling :command:`qstat`. If the number of jobs in the system falls
    below `max_jobs`, more jobs will be enqueued.

    This component operates similarly to :py:class:`.Slice`, in that several
    items are sliced out of an array and submitted as one job. The number of
    jobs in each slice is controlled by the parameter `num_elements`. For the
    sake of efficiency this should be chosen so that the first job submitted
    will still be running when the first batch of jobs finish trickling in.

    .. note::

        This component will wait until all jobs are *enqueued*. You should
        follow it with a :py:class:`.Wait` component if you wish to wait until
        all jobs are *finished*.

    .. warning::

        This component does not submit array jobs, so an explicit `join_var`
        parameter must be passed to any :py:class:`.LoadState` components
        following this.

    :param str split_var: The name of the list in the pipeline state that will
        be split into jobs.

    :param int max_jobs: Maximum number of jobs that will be allowed to run at
        once.

    :param int num_elements: Number of items to process in each job.

    :param float queue_interval: Time in seconds to wait between submitting
        each job.

     :param float poll_interval: Time in seconds between polling qstat while
         waiting.

    .. seealso::

        :py:class:`.BaseQsub`
            For the remaining constructor parameters.
    """
    ADDS = ["qsub_jobs"]
    REMOVES = []

    @property
    def REQUIRED(self):
        return [self.split_var]

    CONFIG_SECTION = "qsub"

    def __init__(self, split_var, max_jobs, num_elements, queue_interval=1.0,
                 poll_interval=30, *qsub_args, **qsub_kwargs):
        super().__init__(*qsub_args, **qsub_kwargs)
        self.split_var = split_var
        self.max_jobs = max_jobs
        self.num_elements = num_elements
        self.queue_interval = queue_interval
        self.poll_interval = poll_interval

    def enqueue(self, pipeline_path, pipeline_state, start_index):
        """
        Enque elements starting at `start_index` with the pipelines state
        `pipeline_state`.

        :return: Remote job details
        :rtype: :py:class:`.RemoteJob`
        """
        end_index = start_index + self.num_elements
        state_slice = pipeline_state[self.split_var][start_index:end_index]

        # Create a copy of the state with the correct slcie and save it.
        sub_state = pipeline_state.copy()
        sub_state[self.split_var] = state_slice
        state_path = self.storage_dir / "state.{}-{}.pickle".format(
            start_index, end_index)
        self._save_state(state_path, sub_state)

        # Generate the job script and run qsub
        job_script = jobscript.StartScript(pipeline_path, state_path)
        job_id = self._qsub(job_script, self.storage_dir)

        return RemoteJob(job_id, state_path)

    def num_running_jobs(self, our_jobs):
        """
        Returns the number of our jobs that are currently enqueued or running.
        """
        current_jobs = running_jobs()
        num_jobs = 0
        for job in our_jobs:
            if job.id in current_jobs and current_jobs[job.id] in {"Q", "R"}:
                num_jobs += 1
        return num_jobs

    def run(self, data, config=None, pipeline=None):
        """Trickle jobs into qsub."""

        qsub_jobs = []

        # The pipeline is static across all jobs, so we can write it now
        pipe_path = self.storage_dir / "pipeline.yml"
        self._save_pipeline(pipe_path, config)

        # Pointer to the next chunk to be enqueued
        slice_index = 0
        while slice_index < len(data[self.split_var]):
            num_jobs = self.num_running_jobs(qsub_jobs)
            if num_jobs < self.max_jobs:
                for _ in range(num_jobs, self.max_jobs):
                    job = self.enqueue(pipe_path, data, slice_index)
                    qsub_jobs.append(job)
                    slice_index += self.num_elements
                    if slice_index >= len(data[self.split_var]):
                        break
                    time.sleep(self.queue_interval)
            else:
                time.sleep(self.poll_interval)

        if "qsub_jobs" not in data:
            data["qsub_jobs"] = []
        data["qsub_jobs"].extend(qsub_jobs)
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

        pipeline_path = job_dir / "pipeline.yml"
        self._save_pipeline(pipeline_path, config)

        if "qsub_jobs" not in data:
            data["qsub_jobs"] = []

        for node in self.nodes:
            pickle_path = job_dir / self._PICKLE_NAME.format(node=node)
            self._save_state(pickle_path, data)
            script = StartScript(str(pipeline_path), str(pickle_path))

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

    .. note::

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
        qsub_jobs = self.get_vals(data)
        self.logger.debug("Waiting for %s to complete", qsub_jobs)
        self._poll(qsub_jobs)
        return data


    def _poll(self, our_jobs):
        while True:
            job_states = set()
            current_jobs = running_jobs()

            for looking_for in our_jobs:
                if looking_for.id in current_jobs:
                    job_states.add(current_jobs[looking_for.id])

            # If all our jobs are in the "C" state or they have disappeared
            # entirely (which happens when the scheduler crashes or the polling
            # interval is too long) then we assume that everything is complete.
            if not job_states or job_states == {"C"}:
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

    :param str join_var: If we are joining multiple jobs (*not* just a single
        array job), this must be defined. Pipeline state is merged on this key.

    :param bool update: If `True`, the pipeline state is updated rather than
        replaced. That is, only fields that were returned from the remote
        job are actually replaced.
    """

    REQUIRED = ["qsub_jobs"]
    REMOVES = ["qsub_jobs"]
    ADDS = []

    def __init__(self, join_var=None, update=True):
        self.join_var = join_var
        self.update = update

    def run(self, data, config=None, pipeline=None):
        if len(data["qsub_jobs"]) > 1 and self.join_var is None:
            raise exception.TooManyjobs(data["qsub_jobs"])

        # Read states even if replace_state is false so that exceptions are
        # thrown if any job is incomplete.
        full_state = _load_state(
            (j.state for j in data["qsub_jobs"]),
            self.join_var)

        del data["qsub_jobs"]
        if self.update:
            data.update(full_state)
            return data
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
        pipeline_path = job_dir / "pipeline.yml"
        stdout_dir = job_dir / "stdout"
        stderr_dir = job_dir / "stderr"
        stdout_dir.mkdir()
        stderr_dir.mkdir()

        self._save_pipeline(pipeline_path, config)
        self._save_state(pickle_path, data)
        script = jobscript.StartScript(str(pipeline_path), str(pickle_path))
        qsub_args = self._depend_args(data)
        job_id = self._qsub(script, stdout_dir, stderr_dir, qsub_args)
        job = RemoteJob(job_id, pickle_path)

        # Remove existing job_ids if they had to complete before this job could
        # start.
        if self.depends or "qsub_jobs" not in data:
            data["qsub_jobs"] = []
        data["qsub_jobs"].append(job)
        return data

def running_jobs():
    """
    Find the state of all jobs in the queue sysetm.

    This function will run :command:`qstat` to find all jobs currently running
    in the queue system.

    :return: Dictionary of job states (e.g. ``R`` for "running", ``Q`` for
        "queued") indexed by job ID.
    :rtype: dict
    """
    qs_results = subprocess.run(
        ["qstat", "-x"],
        check=True, stdout=subprocess.PIPE,
        universal_newlines=True)

    xml_tree = ET.fromstring(qs_results.stdout)
    jobs = {}
    for job in xml_tree.findall("./Job"):
        job_id = job.find("Job_Id").text
        job_state = job.find("job_state").text
        jobs[job_id] = job_state
    return jobs


def _load_state(states, join_var):
    """
    Load state from all the states in the `states` list. Fields are joined
    according to `join_var`. If `join_var` is a list of fields, then all listed
    variables are joined.

    It should be perfectly possible to pass a generator as `states`.
    """

    full_state = None
    # Convert join_var to list
    join_var = [join_var] if isinstance(join_var, str) else join_var

    for partial_state in states:
        if full_state is None:
            full_state = partial_state
        else:
            for var in join_var:
                full_state[var].extend(partial_state[var])
    return full_state
