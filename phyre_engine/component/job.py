"""
This module contains tools for treating the execution of a PhyreEngine pipeline
as a single run of a job.
"""
import os
import tempfile
from phyre_engine.component.component import Component

class Directory(Component):
    """
    Create a unique directory under a given base directory.

    The path of the job directory is stored in the ``job_directory`` pipeline
    variable.

    :param str basedir: Base directory under which to store jobs.
    :param str prefix: Prefix of each directory.
    :param str suffix: Suffix of each directory.
    :param bool chdir: If ``True``,  move to the newly-created directory.
    :param int perms: Use these permissions for the newly-created file.
        By default, permissions are ``0700`` i.e. readable only by the owner.
    """
    REQUIRED = []
    REMOVES = []
    ADDS = ["job_directory"]

    def __init__(
            self, basedir, prefix=None, suffix=None, chdir=False,
            perms=0o700):
        self.basedir = basedir
        self.prefix = prefix
        self.suffix = suffix
        self.chdir = chdir
        self.perms = perms

    def run(self, data, config=None, pipeline=None):
        """Create a unique job directory."""
        job_dir = tempfile.mkdtemp(self.suffix, self.prefix, self.basedir)
        os.chmod(job_dir, self.perms)
        if self.chdir:
            os.chdir(job_dir)
        data["job_directory"] = job_dir
        return data
