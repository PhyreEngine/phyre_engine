"""
Tests for components in the :py:mod:`phyre_engine.component.pymol` module.

These tests will only be run if the ``pymol`` entry is present in the ``tools``
dictionary in the test configuration. The following fields are understood:

``executable``
    Path to the pymol executable to use.

``free_port``
    Network port number known to be free.
"""
import tempfile
import unittest
import xmlrpc.client

import phyre_engine.component.pymol as pymol
import phyre_engine.test
from phyre_engine.test import config as CONFIG

# Pymol command to create a file.
TOUCH = r"""o = open("{}", "w"); o.write("1\n"); o.flush()"""

@phyre_engine.test.requireFields(["executable", "free_port"],
                                 ["tools", "pymol"])
class PymolTestBase(unittest.TestCase):
    """Base class for Pymol tests."""

    @classmethod
    def setUpClass(cls):
        cls.EXECUTABLE = CONFIG["tools"]["pymol"]["executable"]
        cls.PORT = CONFIG["tools"]["pymol"]["free_port"]

    def init(self, command):
        """Return an Init component."""
        return pymol.Init(command, self.EXECUTABLE, self.PORT)

    def quit(self, server):
        """Close connection to pymol server."""
        try:
            server.quit()
        except ConnectionRefusedError:
            # Expected, "quit" causes connection to abort.
            pass


class TestInit(PymolTestBase):
    """Test the Init component."""

    def test_good_port(self):
        """Start Pymol on a known-good port."""
        with tempfile.NamedTemporaryFile("r") as tmp:
            pymol_init = self.init(TOUCH.format(tmp.name))
            pipeline_state = pymol_init.run({})
            self.assertEqual(tmp.read(), "1\n")
            self.assertIn("port", pipeline_state["pymol"])
            self.assertIn("pid", pipeline_state["pymol"])
            self.quit(pymol_init.server)


class TestRun(PymolTestBase):
    """Test the Run component."""

    def test_run_command(self):
        """Run command using Run component."""
        pymol_init = self.init("")
        pipeline_state = pymol_init.run({})
        with tempfile.NamedTemporaryFile("r") as tmp:
            pymol_run = pymol.Run(TOUCH.format(tmp.name))
            pipeline_state = pymol_run.run(pipeline_state)
            self.assertEqual(tmp.read(), "1\n")
            self.quit(pymol_init.server)


class TestQuit(PymolTestBase):
    """Test the Quit component."""

    def test_quit(self):
        """Quit a pymol server using the Quit component."""
        pymol_init = self.init("")
        pymol_quit = pymol.Quit()

        pipeline_state = pymol_init.run({})
        pipeline_state = pymol_quit.run(pipeline_state)
        self.assertNotIn("pymol", pipeline_state)

        with self.assertRaises(ConnectionRefusedError):
            pymol_init.server.do("print 123")
