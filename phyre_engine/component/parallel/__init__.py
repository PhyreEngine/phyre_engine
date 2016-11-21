from abc import ABC, abstractmethod
import math
import pathlib
import pexpect
import pickle
import os
import subprocess
import tempfile
import textwrap
import sys

from phyre_engine.component import Component

class ParallelComponent(Component):
    REQUIRED = []
    ADDS = []
    REMOVES = []

    """Take advantage of PBS to run a function in parallel."""

    TASK_RUNNER = textwrap.dedent(
        """\
        #!/usr/bin/env python
        import sys
        sys.path.extend({path_dirs})
        import phyre_engine.component.parallel
        from {submodule} import {subclass}
        phyre_engine.component.parallel.run_tasks()
        """)

    def __init__(self, component, max_jobs, storage_dir,
            slice_var_in, slice_var_out=None,
            path_dirs=None, submodule=None, subclass=None):
        """Initialise this parallel component.

        Args:
            ``component``: Component to run in parallel.

            ``max_jobs``: Maximum number of jobs to run in parallel.

            ``slice_var_in``: Array variable in data blob on which jobs should
                be split.

            ``slice_var_out``: Array variable on which to join output. If not
                set then ``slice_var_in`` will be used.

            ``storage_dir``: Node-accessible directory in which to store data.

            ``path_dirs``: Array of paths to add to the python module search
                path.

            ``submodule``: Module containing the subclass calling this
                constructor. Defaults to ``self.__module__``, so if the
                subclass is a module on the python path, there is no need to
                pass this parameter. However, if this is called from an
                executable file (such as a unit test), the module may be set to
                ``__main__``, in which case you should pass this parameter.

            ``subclass``: Similar to ``submodule``, except for the class.
                Defaults to ``self.__class__.__name__``.
        """
        self.component = component
        self.max_jobs = max_jobs
        self.slice_var_in = slice_var_in
        self.slice_var_out = slice_var_out if slice_var_out else slice_var_in
        self.storage_dir = storage_dir

        self.path_dirs = path_dirs if path_dirs else []
        self.submodule = submodule if submodule else component.__module__
        self.subclass  = subclass if subclass  else component.__class__.__name__

    def run(self, data):
        """Submit this job to the queue system and wait until all tasks are
        complete.

        This method will submit an array job (``qsub -t``).
        """

        # Note the fragile backslashes at the end of the lines. The "with"
        # statement does not allow items to be split across lines using
        # parentheses, so just be careful of those backslashes.
        #
        # "out_dir" is the output directory in which data is collected.
        # "pickle_file" is the file containing the current pipeline state.
        with \
            tempfile.TemporaryDirectory(dir=self.storage_dir) as out_dir, \
            tempfile.NamedTemporaryFile("wb", dir=self.storage_dir) as pickle_file:

            pickle_data = {
                "self": self,
                "data": data,
                "out_dir": out_dir,
            }
            pickle.dump(pickle_data, pickle_file)
            pickle_file.flush()

            # Number of array elements. If we have fewer tasks to run than
            # self.max_jobs, just run them all. Otherwise, we will run
            # self.max_jobs and process some tasks in serial.
            array_elems = min(self.max_jobs, len(data[self.slice_var_in]))

            with tempfile.NamedTemporaryFile("w") as jobscript:
                jobscript.write(self.TASK_RUNNER.format(
                    interpeter = sys.executable,
                    path_dirs = self.path_dirs,
                    submodule = self.submodule,
                    subclass  = self.subclass
                    ))
                jobscript.flush()

                qsub_cmd = [
                    "qsub",
                    "-t", "0-{}".format(array_elems - 1),
                    "-v", "pickle={}".format(pickle_file.name),
                    "-d", os.getcwd(),
                    jobscript.name
                ]
                process = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE,
                        universal_newlines=True)
                stdout, stderr = process.communicate()
                job_id = stdout.rstrip("\n")

                # Now start an interactive job waiting for the array jobs to
                # finish.  We need to use (p)expect here because qsub checks if
                # it is running on a tty.
                monitor_args = [
                    "-W", "depend=afteranyarray:{0}".format(job_id),
                    "-I"
                ]
                waiter = pexpect.spawn("qsub", monitor_args)
                waiter.expect("^qsub: waiting for job \S+ to start\r\n")
                waiter.expect("^qsub: job \S+ ready\r\n", timeout=None)
                waiter.sendline("exit")

            # At this point, all the node jobs should have run. We now just need to
            # collect the data from each one.
            collected = []
            for i in range(0, array_elems):
                file = pathlib.Path(out_dir, str(i))
                if not file.exists():
                    raise NoOutputError(i, file)
                with file.open("rb") as f:
                    child_data = pickle.load(f)
                    collected.extend(child_data)
            data[self.slice_var_out] = collected
            return data

class NoOutputError(Exception):
    def __init__(self, child, file):
        err_msg = "No output from child {} (file: {})"
        super().__init__(err_msg.format(child, file))

def run_tasks():
    pickle_data = pickle.load(open(os.environ["pickle"], "rb"))

    obj     = pickle_data["self"]
    data    = pickle_data["data"]
    out_dir = pickle_data["out_dir"]

    pbs_arrayid = int(os.environ["PBS_ARRAYID"])
    elems_per_job = math.ceil(len(data[obj.slice_var_in]) / obj.max_jobs)
    i = int(pbs_arrayid) * elems_per_job

    data[obj.slice_var_in] = data[obj.slice_var_in][i:i + elems_per_job]

    results = obj.component.run(data)

    out_file = pathlib.Path(out_dir, str(pbs_arrayid))
    with out_file.open("wb") as out:
        pickle.dump(results[obj.slice_var_out], out)
