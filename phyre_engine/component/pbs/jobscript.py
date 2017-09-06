"""
Module containing classes that, when converted to a string, give a jobscript for
running pipelines or sub-pipelines via qsub.
"""

import textwrap

class JobScript:
    r"""
    Base class for jobscripts that, when converted to a string, give a
    jobscript for running the parallel jobs started by :py:class:`.Qsub`.

    :param str python: Path to the python interpreter to use.
        Default: "python".
    :param list python_path: Extra directories to append to the
        ``PYTHONPATH`` environment variable.
    """

    HEADER = textwrap.dedent("""\
    #!/bin/bash
    # Job script run by phyre_engine.component.pbs.qsub
    """)

    PATH_SETUP = 'export PYTHONPATH="{python_path}:$PYTHONPATH"'

    def __init__(self, python="python", python_path=None):
        self.python = python
        self.python_path = python_path

    def __str__(self):
        script = self.HEADER
        if self.python_path is not None:
            script += self.PATH_SETUP.format(python_path=self.python_path)
        return script

class StartScript(JobScript):
    MODULE = "phyre_engine.component.pbs.run"
    RUN_COMMAND = "'{python}' -m{module} {pipeline} {state}"

    def __init__(self, pipeline, state, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pipeline = pipeline
        self.state = state

    def __str__(self):
        script = super().__str__()
        script += self.RUN_COMMAND.format(
            python=self.python,
            module=self.MODULE,
            pipeline=self.pipeline,
            state=self.state)
        return script
