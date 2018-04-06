"""PhyreEngine installer."""
from setuptools import setup, find_packages
import configparser
import os
import subprocess

def _read_setup_cfg_version():
    """Returns the version number, or raises a KeyError if it is not set."""
    setup_cfg = configparser.ConfigParser()
    setup_cfg.read("setup.cfg")
    try:
        return setup_cfg["metadata"]["version"]
    except KeyError:
        return None


def _write_setup_cfg(version):
    """Read ``setup.cfg``, set version number and overwrite."""
    setup_cfg = configparser.ConfigParser()
    setup_cfg.read("setup.cfg")
    setup_cfg.setdefault("metadata", {})["version"] = version
    with open("setup.cfg", "w") as setup_out:
        setup_cfg.write(setup_out)


def _git_version_string():
    """
    Call ``git describe --tags`` to get the version number.

    Any leading "v" is stripped, and all hyphens replaced with underscores.
    """
    try:
        git_result = subprocess.run(
            ["git", "describe", "--tags"],
            stdout=subprocess.PIPE, check=True, universal_newlines=True)
        git_version = git_result.stdout.strip()

        if git_version[0] == "v":
            git_version = git_version[1:]
        git_version = git_version.replace("-", "_")
        return git_version
    except subprocess.CalledProcessError:
        return None


def _conda_version_string():
    """Return the ``PKG_VERSION`` environment var or raise a KeyError."""
    try:
        return os.environ["PKG_VERSION"]
    except KeyError:
        return None

def version():
    """
    Deduce version number.

    This works in three stages:

    1. The ``PKG_VERSION`` environment variable is used if it is set. This is
       the version number set by Conda.

    2. If a call to ``git describe`` succeeds, the output of that is used. A
       leading "v" is stripped, and any dashes are replaced with underscores.

    3. Finally, attempt to read a version string from ``setup.cfg``.

    If either of the first two methods succeed, the results are also written to
    ``setup.cfg``.
    """
    version_getters = (
            _conda_version_string,
            _git_version_string,
            _read_setup_cfg_version)

    ver = None
    for getter in version_getters:
        ver = getter()
        if ver is not None:
            break
    _write_setup_cfg(ver)
    return ver

setup(
    name="PhyreEngine",
    version=version(),
    packages=find_packages(),
    include_package_data=True,
    author="Stefans Mezulis",
    author_email="stefans.mezuli08@imperial.ac.uk",
    description="PhyreEngine is a tool for building bioinformatics pipelines.",
    url="http://www.sbg.bio.ic.ac.uk/phyreengine",
    install_requires=[
        "biopython>=1.71",
        "numpy",
        "scipy",
        "PyYAML",
        "jmespath>=0.9.3",
        "appdirs",
    ],
    tests_require=[
        "xmlrunner",
    ],
    entry_points={
        "console_scripts": [
            "phyre_engine = phyre_engine.run:main"
    ]},
    test_suite="phyre_engine.test"
)
