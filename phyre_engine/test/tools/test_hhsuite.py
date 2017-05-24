import unittest
import subprocess

import phyre_engine.tools.hhsuite as hhsuite

# Example tool for testing
class ExampleTool(hhsuite.HHSuiteTool):
    flags = {
        "test1": "t",
        "test2": "test2",
        "test3": "-test3",
        "test4": "--test4",
    }
    def __init__(self, program="test_program", *args, **flags):
        super().__init__(program, *args, **flags)

class TestHHSuite(unittest.TestCase):
    """Test common methods for hhsuite tools."""

    def verify_cmd_param(self, cmd, param, value=None):
        self.assertIn(param, cmd, "{} is in {}".format(param, cmd))

        param_index = cmd.index(param)
        if value is not None:
            self.assertEqual(
                cmd[param_index + 1], value,
                "Value of {}".format(param))

    def test_cmd_line(self):
        """Check that command lines are processed correctly."""

        test_tool = ExampleTool(
                "/path/to/test_program",
                "positional_arg_1", "positional_arg_2",
                test1="val_1", test2="val_2",
                test3=True, test4=123,
                test5=456, b=789,
                **{"-----test-long":"101112"})

        cmd = test_tool.command_line

        program = cmd.pop(0)
        self.assertEqual(
            program, "/path/to/test_program",
            "Can specify program with path")
        self.verify_cmd_param(cmd, "-t", "val_1")
        self.verify_cmd_param(cmd, "--test2", "val_2")
        self.verify_cmd_param(cmd, "-test3", None)
        self.verify_cmd_param(cmd, "--test4", "123")
        self.verify_cmd_param(cmd, "--test5", "456")
        self.verify_cmd_param(cmd, "-b", "789")
        self.verify_cmd_param(cmd, "-----test-long", "101112")
        self.assertEqual(
            cmd[-2], "positional_arg_1",
            "Positional arg appended")
        self.assertEqual(
            cmd[-1], "positional_arg_2",
            "Positional arg appended")

    # These sanity checks are pretty much just so we can check that the program
    # name is correct and the class can be instantiated properly. Seems trivial,
    # but it's bitten me in the past. Very easy to copy/paste a new tool and
    # forget to change the program name.

    def test_hhblits_cmd_line(self):
        """HHBlits sanity check"""
        hhblits = hhsuite.HHBlits(database="test_db")
        self.assertListEqual(
            hhblits.command_line,
            ["hhblits", "-d", "test_db"],
            "Sanity check on hhblits command line")

    def test_hhsearch_cmd_line(self):
        """HHSearch sanity check"""
        hhsearch = hhsuite.HHSearch(output="test_out")
        self.assertListEqual(
            hhsearch.command_line,
            ["hhsearch", "-o", "test_out"],
            "Sanity check on hhsearch command line")

    def test_hhmake_cmd_line(self):
        """HHMake sanity check"""
        hhmake = hhsuite.HHMake("input", output="test_out", verbose=3)

        program = hhmake.command_line.pop(0)
        self.assertEqual(program, "hhmake", "correct program")
        self.verify_cmd_param(hhmake.command_line, "-i", "input")
        self.verify_cmd_param(hhmake.command_line, "-o", "test_out")
        self.verify_cmd_param(hhmake.command_line, "-v", "3")

    def test_cstranslate_cmd_line(self):
        """cstranslate sanity check"""
        cstranslate = hhsuite.CSTranslate(
                infile="in_file",
                pc_admix=0,
                binary=True)
        program = cstranslate.command_line.pop(0)
        self.assertEqual(program, "cstranslate", "correct program")
        self.verify_cmd_param(cstranslate.command_line, "--infile", "in_file")
        self.verify_cmd_param(cstranslate.command_line, "--pc-admix", "0")
        self.verify_cmd_param(cstranslate.command_line, "--binary")

    def test_ffindex_build_cmd_line(self):
        """Sanity check on ffindex_build command line."""
        ffindex_build = hhsuite.FFIndexBuild("DATA", "INDEX", append=True)
        self.assertListEqual(
            ffindex_build.command_line,
            ["ffindex_build", "-a", "DATA", "INDEX"],
            "Sanity check on ffindex_build command line")
