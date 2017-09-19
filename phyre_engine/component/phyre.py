"""
This module contains components that are specific to the Phyre protein structure
prediction server. They may be useful in other contexts, but less care is paid
to ensuring that components in this module are reusable and generic.
"""
from datetime import datetime
import json
import os
import pathlib
import time

from phyre_engine.component.component import Component

class JobInfoComponent(Component):
    """
    Base class for components that append keys to a JSON file.

    :param str file: System information is appended to this file.
    """
    REQUIRED = []
    ADDS = []
    REMOVES = []

    def __init__(self, file="sysinfo.json"):
        self.file = pathlib.Path(file)

    def write(self, fields):
        """Append fields to a JSON file."""
        file_contents = {}

        # We want to append to the file, so read it first.
        if self.file.exists():
            with self.file.open("r") as json_in:
                file_contents = json.load(json_in)

        file_contents.update(fields)
        with self.file.open("w") as json_out:
            json.dump(file_contents, json_out, indent=4)


class StartingJobInfo(JobInfoComponent):
    """
    Write PBS job ID to the ``pbs_id`` field of a JSON file.
    """
    REQUIRED = ["qsub_jobs"]

    def run(self, data, config=None, pipeline=None):
        self.write({"pbs_id": data["qsub_jobs"][0].id})
        return data

class StartedJobInfo(JobInfoComponent):
    """
    Record various pieces of information about the system to a file. The
    following fields will be appended to the file.

    :hostname: The host name of the current system, as reported by
        ``os.uname()``

    :start_unix_time: The current time, in seconds since the UNIX epoch.

    :start_human_time: Human-readable time in ISO 8601 format.
    """

    def run(self, data, config=None, pipeline=None):
        """Reading system information."""
        fields = {
            "hostname": os.uname().nodename,
            "start_unix_time": time.time()
        }
        fields["start_human_time"] = datetime.fromtimestamp(
            fields["start_unix_time"]).isoformat()
        self.write(fields)
        return data

class EndedJobInfo(JobInfoComponent):
    """
    Write ``end_unix_time`` and ``end_human_time`` to a file.
    """

    def run(self, data, config=None, pipeline=None):
        fields = {"end_unix_time": time.time()}
        fields["end_human_time"] = datetime.fromtimestamp(
            fields["end_unix_time"]).isoformat()
        self.write(fields)
        return data
