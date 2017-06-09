#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.run -- Run a PhyreEngine pipeline

PhyreEngine is a framework for building and running common structural biology
pipelines. This script is responsible for reading, setting up and running a
PhyreEngine pipeline. Unless you're an expert, this script is probably all that
you need to run a pipeline.
'''

import logging.config
import sys

from argparse import ArgumentParser, Action
from phyre_engine.pipeline import Pipeline

try:
    # Use libyaml if it is available
    import yaml
    from yaml import CSafeLoader as SafeLoader, CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader, SafeDumper

class DumpAction(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        config = {}
        config.update({"logging": default_log_config()})
        config.update(dummy_pipeline())
        yaml.dump(config, sys.stdout, SafeDumper, default_flow_style=False)
        sys.exit(0)

def arg_parser():
    # Setup argument parser
    parser = ArgumentParser(description=__import__('__main__').__doc__)

    parser.add_argument(
        "-v", "--verbose", dest="verbose", action="count", default=1,
        help="set verbosity level [default: %(default)s]")
    parser.add_argument(
        "-e", "--example", dest="dump", action=DumpAction, nargs=0,
        help="Dump a sample pipeline and exit.")
    parser.add_argument(
        dest="pipeline", metavar="pipeline",
        help="YAML file describing the pipeline")
    return parser

def default_log_config():
    return {
        'version': 1,
        'formatters': {
            'simple': {
                'format': '{levelname} : {asctime} : {name} : {message}',
                'style': '{',
            },
        },
        'handlers': {
            'console': {
                'class': 'logging.StreamHandler',
                'level': 'WARNING',
                'formatter': 'simple',
                'stream': 'ext://sys.stderr',
            },
        },
        'root': {
            'level': 'WARNING',
            'handlers': ['console']
        }
    }

def dummy_pipeline():
    return {
        'pipeline': {
            'checkpoint': 'checkpoint_file.chk',
            'components': [
                'phyre_engine.component.FastaInput.FastaInput',
                'phyre_engine.component.SeqValidator.SeqValidator', {
                    'phyre_engine.component.hhsuite.HHBlits': [
                        '/path/to/data/uniclust30',
                        {'iterations': 3},
                    ]
                }
            ]
        }
    }

def init_logging(logging_dict):
    if logging_dict is not None:
        logging.config.dictConfig(logging_dict)
    else:
        logging.config.dictConfig(default_log_config())

def main():  # IGNORE:C0111
    '''Command line options.'''

    try:
        parser = arg_parser()
        args = parser.parse_args()
        if args.dump:
            return 0

        # Parse pipeline descriptor from YAML file
        with open(args.pipeline, "r") as yml_in:
            config = yaml.load(yml_in, SafeLoader)

        # Set up logging if a logging section was given in the pipeline
        init_logging(config.get("logging", None))

        # Load a pipeline
        pipeline = Pipeline.load(config["pipeline"])

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())
