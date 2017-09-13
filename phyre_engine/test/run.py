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

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeLoader as SafeLoader
except ImportError:
    from yaml import SafeLoader

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

def construct_yaml_tuple(self, node):
    """
    Used to convert sequences from lists to tuples. Only applies to lists
    without any nested structures.
    """
    # TODO: This is shared with phyre_engine.run: refactor shared code
    seq = self.construct_sequence(node)
    if any(isinstance(e, (list, tuple, dict)) for e in seq):
        return seq
    return tuple(seq)

def _main():
    # Capture warnings with the logger
    logging.captureWarnings(True)
    try:
        parser = arg_parser()
        args = parser.parse_args()

        if args.config is not None:
            with open(args.config, "r") as yml_in:
                SafeLoader.add_constructor(
                    'tag:yaml.org,2002:seq',
                    construct_yaml_tuple)
                phyre_engine.test.config = yaml.load(yml_in, SafeLoader)

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
