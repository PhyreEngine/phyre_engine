#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.test.run -- Run PhyreEngine tests

This module is the entry point to running unit tests for PhyreEngine. This
module allows a configuration file to be specified which tells tests where
to find various system dependencies. Without setting a configuration file using
this script, long-running tests or tests with dependencies on external tools
will not be run.
'''

import logging
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
import unittest
import phyre_engine.test
import phyre_engine.tools.yaml as yaml


def arg_parser():
    """Set up argument parser."""
    parser = ArgumentParser(
        description=__import__('__main__').__doc__,
        formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument(
        "-v", "--verbose", dest="verbosity", action="count", default=1,
        help="set verbosity level [default: %(default)s]")
    parser.add_argument(
        "-d", "--test-dir", dest="test_dir", default=str(Path(__file__).parent),
        help="Directory in which to search for tests. [default: %(default)s]")
    parser.add_argument(
        "-p", "--pattern", dest="pattern", default="test*.py",
        help="Pattern matching test files [default: %(default)s].")
    parser.add_argument(
        "-c", "--config", dest="config",
        help="YAML configuration file.")
    return parser


def _main():
    # Capture warnings with the logger
    logging.captureWarnings(True)
    try:
        parser = arg_parser()
        args = parser.parse_args()

        if args.config is not None:
            with open(args.config, "r") as yml_in:
                phyre_engine.test.config = yaml.load(yml_in)

        loader = unittest.TestLoader()
        tests = loader.discover(args.test_dir, args.pattern)
        runner = unittest.TextTestRunner(verbosity=args.verbosity)
        runner.run(tests)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(_main())
