"""
Tools for checking a git repository for changes.
"""
from phyre_engine.component.component import Component
import os
import contextlib
import subprocess
from pathlib import Path
import sys

#: Value of the ``source`` attribute of each element of ``update_required``.
GIT_SOURCE = "git"

class Check(Component):
    """
    Check for updates to a git repository.

    This module will fetch (but not merge) any new updates to a git repository.
    If updates are available, the tool will be added to the ``update_required``
    list in the pipeline state. Each element added to the ``update_required``
    list is a dictionary containing the following keys:

    new_version:
        ID of the newest commit.

    repo_dir:
        Directory containing the git repository (the ``src_dir`` attribute of
        the constructor).

    name:
        Name of the software (the ``name`` attribute of the constructor).

    source:
        Method used to retrieve software. In this case, "git".

    Any errors encountered when running git will result in a
    :py:class:`subprocess.CalledProcessError` being raised.

    :param str name: Name of the software.
    :param str repo_dir: Source directory containing git repository.
    :param str git: Git executable.
    :param bool force: Always act as though an update is available.

    .. note::

        This class assumes that the git repository has previously been checked
        out: that is, that we can do a ``git fetch`` rather than having to do a
        ``git clone``.

    .. note::

        This module directly calls the ``git`` executable: ``git`` must be
        installed.
    """
    ADDS = ["update_required"]
    REMOVES = []
    REQUIRED = []

    _GIT_FETCH = ["fetch", "origin"]
    _GIT_DIFF = ["diff", "--quiet", "HEAD", "origin/master"]
    _GIT_GET_HEAD = ["rev-list", "--max-count=1", "--abbrev-commit",
                     "HEAD", "origin/master"]

    _NOT_FOUND_ERR = (
        "Source directory {src_dir} does not exist. You must clone and "
        "configure the repository before it can be updated by this component.")

    def __init__(self, name, src_dir, git="git", force=False):
        self.name = name
        self.src_dir = Path(src_dir)
        self.git = git
        self.force = force

    def run(self, data, config=None, pipeline=None):
        """
        :raises FileNotFoundError: If the configured ``src_dir`` does not exist.
        """
        if not Path(self.src_dir).exists():
            raise FileNotFoundError(self._NOT_FOUND_ERR.format(self.src_dir))

        # Create the update_required list even if no update is found, because it
        # is probably required by components later in the pipeline.
        if "update_required" not in data:
            data["update_required"] = []

        new_commit = self.check()
        if self.force or new_commit is not None:
            data["update_required"].append({
                "new_version": new_commit,
                "name": self.name,
                "repo_dir": self.src_dir,
                "source": GIT_SOURCE
            })

        return data

    def check(self):
        """
        Returns the hash of the most recent commit if updates are available.
        """
        # Check using git whether an update is available
        with chdir(self.src_dir):
            subprocess.run([self.git] + self._GIT_FETCH, check=True)
            diff_result = subprocess.run(
                [self.git] + self._GIT_DIFF,
                check=False)
            if diff_result.returncode == 1:
                # There was a difference: return the commit ID
                rev_list_result = subprocess.run(
                    [self.git] + self._GIT_GET_HEAD,
                    check=True, stdout=subprocess.PIPE)
                commit = rev_list_result.stdout.decode(sys.stdout.encoding)
                return commit.strip()
            else:
                return None

class Get(Component):
    """
    Update git repositories.

    :param bool stash: If true, stash changes before merge and then reapply
        them.
    :param str git: Git executable.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["update_required"]

    _GIT_STASH_PUSH = ["stash"]
    _GIT_STASH_POP = ["stash", "pop"]
    _GIT_MERGE = ["merge", "origin/HEAD"]

    def __init__(self, stash=True, git="git"):
        self.stash = stash
        self.git = git

    def run(self, data, config=None, pipeline=None):
        for tool in data["update_required"]:
            if tool["source"] != GIT_SOURCE:
                continue

            with chdir(tool["repo_dir"]):
                if self.stash:
                    subprocess.run(
                        [self.git] + self._GIT_STASH_PUSH,
                        check=True)
                subprocess.run([self.git] + self._GIT_MERGE, check=True)
                if self.stash:
                    subprocess.run([self.git] + self._GIT_STASH_POP, check=True)
        return data

@contextlib.contextmanager
def chdir(path):
    """
    Context manager for restoring the working directory after changing
    directory.
    """
    orig_dir = os.getcwd()
    try:
        os.chdir(str(path))
        yield
    finally:
        os.chdir(orig_dir)
