"""PhyreEngine installer."""
from setuptools import setup, find_packages
import configparser
import subprocess

def version():
    """
    Use ``git describe`` to get a version number and write it to ``setup.cfg``.

    If ``git describe`` cannot be called do nothing and let ``setup.cfg`` give
    the version.
    """
    setup_cfg = configparser.ConfigParser()
    setup_cfg.read("setup.cfg")
    if "metadata" not in setup_cfg:
        setup_cfg["metadata"] = {}

    try:
        git_result = subprocess.run(
            ["git", "describe", "--tags"],
            stdout=subprocess.PIPE, check=True)

        git_version = git_result.stdout.decode("UTF-8").strip()

        # Strip leading "v" from git tag
        if git_version[0] == "v":
            git_version = git_version[1:]

        # Write to setup.cfg
        setup_cfg["metadata"]["version"] = git_version
        with open("setup.cfg", "w") as setup_out:
            setup_cfg.write(setup_out)
    except Exception as e:
        # Not a problem, just use the value in setup.cfg
        pass

    return setup_cfg["metadata"]["version"]

setup(
    name="PhyreEngine",
    version=version(),
    packages=find_packages(exclude=[
        "phyre_engine.test",
        "phyre_engine.test.*"
    ]),
    author="Stefans Mezulis",
    author_email="stefans.mezuli08@imperial.ac.uk",
    description="PhyreEngine is a tool for building bioinformatics pipelines.",
    url="http://www.sbg.bio.ic.ac.uk/phyreengine",
    install_requires=[
        "biopython>=1.70",
        "numpy",
        "scipy",
        "PyYAML",
        "jmespath>=0.9.3"
    ],
    entry_points={
        "console_scripts": [
            "phyre_engine = phyre_engine.run:main"
    ]},
    test_suite="phyre_engine.test"
)
