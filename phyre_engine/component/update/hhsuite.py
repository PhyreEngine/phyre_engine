"""
Tools for updating an hhsuite version.
"""
from phyre_engine.component.component import Component
import os
import contextlib
import subprocess
from pathlib import Path
import tempfile

class GitUpdate(Component):
    """
    Check for updates to the hh-suite git repository.

    This module will fetch any new updates to hh-suite. If any updates have
    occurred, it will stash any local changes in order to preserve any
    configuration information (such as HHPaths.pm or any changes to
    CMakeLists.txt), merge the remote branch, and reapply the stashed changes.

    Any errors encountered when running git will result in a
    :py:class:`subprocess.CalledProcessError` being raised.

    :param str src_dir: Source directory containing hh-suite.

    .. warning::

        This component will likely be fairly fragile: if non-trivial changes
        have been made to any files in the source tree, merging with the
        upstream copy or popping the stashed changes may well result in merge
        conflicts. It is strongly recommended that your production setup not
        rely on the source directory being sane. That is, breaking the source
        directory should not affect the installed binaries of hh-suite.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    _GIT_FETCH = ["git", "fetch", "origin"]
    _GIT_DIFF = ["git", "diff", "--quiet", "HEAD", "origin/master"]
    _GIT_SUBMODULE_INIT = ["git", "submodule", "init"]
    _GIT_SUBMODULE_UPDATE = ["git", "submodule", "update"]
    _GIT_STASH_PUSH = ["git", "stash"]
    _GIT_MERGE = ["git", "merge", "origin/master"]
    _GIT_STASH_POP = ["git", "stash", "pop"]

    _NOT_FOUND_ERR = (
        "Source directory {src_dir} does not exist. You must clone and "
        "configure hh-suite before it can be updated by this component.")

    def __init__(self, src_dir):
        self.src_dir = src_dir

    def run(self, data, config=None, pipeline=None):
        """
        :raises FileNotFoundError: If the configured ``src_dir`` does not exist.
        """
        if not Path(self.src_dir).exists():
            raise FileNotFoundError(self._NOT_FOUND_ERR.format(self.src_dir))
        self._update()
        return data

    def _update(self):
        # Check using git whether an update is available
        with chdir(self.src_dir):
            subprocess.run(self._GIT_FETCH, check=True)
            diff_result = subprocess.run(self._GIT_DIFF, check=False)
            if diff_result.returncode == 1:
                # There was a diff, so stash our local changes
                # (scripts/HHPaths.pm and any other configuration files), merge
                # the master, and pop our stash.
                subprocess.run(self._GIT_STASH_PUSH, check=True)
                subprocess.run(self._GIT_MERGE, check=True)
                subprocess.run(self._GIT_STASH_POP, check=True)
                subprocess.run(self._GIT_SUBMODULE_INIT, check=True)
                subprocess.run(self._GIT_SUBMODULE_UPDATE, check=True)


class Compile:
    """
    Compile hh-suite.

    HH-suite will be compiled in a temporary directory by calling ``cmake``.
    Arbitrary options may be added to the cmake command line by adding them to
    the ``cmake_opts`` parameter *without* a ``-D`` prefix (i.e. ``FOO: "bar"``
    will add the parameter ``-DFOO=bar`` to the command line). Compilation will
    be done by ``make``. Hh-suite will be installed to ``src_dir`` by setting
    the ``CMAKE_INSTALL_PREFIX`` flag and running ``make install``.

    :param str src_dir: Directory containing hh-suite source.
    :param str dest_dir: Installation prefix.
    :param dict cmake_opts: Dictionary containing options to pass to ``cmake``,
        without the ``-D`` prefix.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, src_dir, dest_dir, cmake_opts=None):
        self.src_dir = src_dir
        self.dest_dir = dest_dir
        self.cmake_opts = cmake_opts if cmake_opts is not None else {}

    def run(self, data, config=None, pipeline=None):
        with tempfile.TemporaryDirectory() as tmpdir:
            with chdir(tmpdir):
                subprocess.run(self._cmake_cmd_line(), check=True)
                subprocess.run(["make"], check=True)
                subprocess.run(["make", "install"], check=True)
        return data

    def _cmake_cmd_line(self):
        cmd_line = [
            "cmake",
            "-DCMAKE_INSTALL_PREFIX:PATH={}".format(self.dest_dir)]
        for var, value in self.cmake_opts.items():
            cmd_line.append("-D{}={}".format(var, value))
        cmd_line.append(self.src_dir)
        return cmd_line


@contextlib.contextmanager
def chdir(path):
    """
    Context manager for restoring the working directory after changing
    directory.
    """
    orig_dir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(orig_dir)
