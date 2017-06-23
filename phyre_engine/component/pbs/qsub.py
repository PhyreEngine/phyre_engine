from . import jobscript
from phyre_engine.component.component import Component
import os
import math
import copy
import pickle
from pathlib import Path
import subprocess
import tempfile
import phyre_engine
import collections
from phyre_engine.pipeline import ExpectedExit
import xml.etree.ElementTree as ET
import sys

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

#: Save the current pipeline and index to resume after workers are complete.
PipeState = collections.namedtuple(
    "PipeState", ["pipeline", "pipeline_index", "config"])

# File name templates
WORKER_STATE = "state-{}.pickle"
RESUME_PIPELINE = "pipeline.pickle"

class Qsub(Component):
    r"""
    Run a pipeline in parallel on separate nodes.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(
            self, slice_in, slice_out, max_jobs, storage_dir, sub_pipeline,
            name="phyre_engine", log_config=None, qsub_args=None):

        self.slice_in = slice_in
        self.slice_out = slice_out
        self.max_jobs = max_jobs
        self.storage_dir = Path(storage_dir)
        self.name = name
        self.log_config = log_config
        self.qsub_args = qsub_args

        if not isinstance(sub_pipeline, phyre_engine.Pipeline):
            sub_pipeline = phyre_engine.Pipeline.load(sub_pipeline)
        self.sub_pipeline = sub_pipeline

    def _slice_state(self, data, config):
        """
        Slice pipeline and save into a series of pickles containing the subset
        of data on which each remote job will operate. Each pickle contains a
        named tuple containing the elements "pipeline" and "data".
        """

        data_len = len(data[self.slice_in])
        # Number of array elements. If we have fewer tasks to run than
        # self.max_jobs, just run them all. Otherwise, we will run self.max_jobs
        # and process some tasks in serial.
        num_jobs = min(self.max_jobs, data_len)
        elems_per_job = math.ceil(data_len / num_jobs)

        # Overwrite the logging section of the pipeline config if it is set
        config = copy.copy(config)
        if self.log_config is not None:
            config["logging"] = self.log_config

        for i in range(0, num_jobs):
            file_name = WORKER_STATE.format(i)

            start = i * elems_per_job
            end = (i + 1) * elems_per_job

            data_slice = data[self.slice_in][start:end]
            state = copy.copy(data)
            state[self.slice_in] = data_slice

            with Path(self.storage_dir, file_name).open("wb") as state_out:
                sub_pipe = copy.copy(self.sub_pipeline)
                sub_pipe.config = config
                sub_pipe.start = state

                pickle.dump(SubPipeline(sub_pipe), state_out)
        return num_jobs

    def _save_pipeline(self, config, pipeline):
        """
        Save the pipeline and current component index so we can resume later.
        """
        pipeline_index = [i for i, c
                          in enumerate(pipeline.components)
                          if c is self][0]

        state = PipeState(
            pipeline=pipeline,
            pipeline_index=pipeline_index,
            config=config)
        with (self.storage_dir / RESUME_PIPELINE).open("wb") as out_fh:
            pickle.dump(state, out_fh)

    def run(self, data, config, pipeline):
        """
        Start an array job (``-t``) with qsub, and start a job to restart the
        pipeline that depends on the array job.
        """

        self._save_pipeline(config, pipeline)
        num_jobs = self._slice_state(data, config)
        main_job_id = self._qsub_start(num_jobs)
        resume_job_id = self._qsub_resume(num_jobs, main_job_id)
        raise ExpectedExit((
            "Stopping master process and running worker job {}. Resuming as "
            "PBS job {}.").format(main_job_id, resume_job_id))

    def _qsub_start(self, num_jobs):
        start_qsub_cmd = [
            "qsub",
            "-t", "0-{}".format(num_jobs - 1),
            "-o", str(self.storage_dir),
            "-e", str(self.storage_dir),
            "-d", os.getcwd(),
            "-N", self.name,
        ]
        start_qsub_cmd.extend(self.qsub_args)

        with tempfile.NamedTemporaryFile("w") as jobfile:
            # Write jobscript to temp file and append it to the command line.
            # The pickle file contains the shell variable $PBS_ARRAYID.
            pickle_file =  '"{storage}/state-${{PBS_ARRAYID}}.pickle"'.format(
                storage=str(self.storage_dir))

            script = jobscript.StartScript(self.storage_dir, pickle_file)
            jobfile.write(str(script))
            jobfile.flush()
            start_qsub_cmd.append(jobfile.name)

            process = subprocess.Popen(
                start_qsub_cmd,
                stdout=subprocess.PIPE,
                universal_newlines=True)

            stdout, _ = process.communicate()
            job_id = stdout.rstrip("\n")
            return job_id

    def _qsub_resume(self, num_jobs, job_id):
        resume_qsub_cmd = [
            "qsub",
            "-d", os.getcwd(),
            "-N", "{}-resume".format(self.name),
            "-W", "depend=afteranyarray:{0}".format(job_id),
        ]
        resume_qsub_cmd.extend(self.qsub_args)

        with tempfile.NamedTemporaryFile("w") as jobfile:
            # Write jobscript to temp file and append it to the command line
            script = jobscript.ResumeScript(
                num_jobs, self.slice_out, self.storage_dir)
            jobfile.write(str(script))
            jobfile.flush()
            resume_qsub_cmd.append(jobfile.name)

            process = subprocess.Popen(
                resume_qsub_cmd,
                stdout=subprocess.PIPE,
                universal_newlines=True)

            stdout, _ = process.communicate()
            job_id = stdout.rstrip("\n")
            return job_id

class ContactNodes(Component):
    ADDS = []
    REQUIRED = []
    REMOVES = []

    def __init__(
            self, storage_dir, sub_pipeline, nodes=None,
            name="phyre_engine", log_config=None, qsub_args=None):

        self.storage_dir = Path(storage_dir)
        self.name = name
        self.log_config = log_config
        self.qsub_args = qsub_args

        if nodes is None:
            nodes = self._online_nodes()
        self.nodes = nodes

        if not isinstance(sub_pipeline, phyre_engine.Pipeline):
            sub_pipeline = phyre_engine.Pipeline.load(sub_pipeline)
        self.sub_pipeline = sub_pipeline

    def run(self, data, config=None, pipeline=None):
        for node in self.nodes:
            serialised_state = self._save_pipeline(node, data, config)
            self._qsub_start(serialised_state, node)
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


    def _save_pipeline(self, node, data, config):
        """
        Save the current pipeline to a pickle named after the node that we are
        going to be contacting.
        """

        file_name = "state-{}".format(node)

        # Overwrite the logging section of the pipeline config if it is set
        config = copy.copy(config)
        if self.log_config is not None:
            config["logging"] = self.log_config

        pickle_path = Path(self.storage_dir, file_name)
        with pickle_path.open("wb") as state_out:
            sub_pipe = copy.copy(self.sub_pipeline)
            sub_pipe.config = config
            sub_pipe.start = data

            pickle.dump(SubPipeline(sub_pipe), state_out)
        return pickle_path

    def _qsub_start(self, pickled_state, node):
        start_qsub_cmd = [
            "qsub",
            "-o", str(self.storage_dir),
            "-e", str(self.storage_dir),
            "-d", os.getcwd(),
            "-N", self.name,
            "-lnodes={}".format(node)
        ]
        start_qsub_cmd.extend(self.qsub_args)

        with tempfile.NamedTemporaryFile("w") as jobfile:
            # Write jobscript to temp file and append it to the command line
            script = jobscript.StartScript(self.storage_dir, pickled_state)
            jobfile.write(str(script))
            jobfile.flush()
            start_qsub_cmd.append(jobfile.name)

            process = subprocess.Popen(
                start_qsub_cmd,
                stdout=subprocess.PIPE,
                universal_newlines=True)

            stdout, _ = process.communicate()
            job_id = stdout.rstrip("\n")
            return job_id
