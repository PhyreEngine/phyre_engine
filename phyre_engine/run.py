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
        yaml.dump(
            dummy_pipeline(),
            sys.stdout, SafeDumper, default_flow_style=False)
        sys.exit(0)

class StoreStartingValue(Action):
    """
    Parse a starting value for the pipeline state.

    For example, assume that this action is associated with the ``--start``
    option and the ``start`` ``dest``. The user passes the following options:

    .. code-block:: none
        --start foo:bar --start baz:qux1 --start baz:qux2

    This will result in the namespace returned by the argument parser having a
    ``start`` (set from the ``dest`` parameter) attribute that looks like the
    following:

    .. code-block::

       {"foo": "bar", "baz": ["qux1", "qux2"]}

    Values for this option are split on the first ``:``. If several of the same
    values are passed, a list is built.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        start_vals = getattr(namespace, self.dest, None)
        if start_vals is None:
            start_vals = {}
        key, value = values[0].split(":", maxsplit=1)

        if key not in start_vals:
            start_vals[key] = value
        else:
            # Convert to list if we saw a scalar
            if not isinstance(start_vals[key], list):
                start_vals[key] = [start_vals[key]]
            start_vals[key].append(value)
        setattr(namespace, self.dest, start_vals)


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
        "-s", "--start", dest="start", action=StoreStartingValue, nargs=1,
        help="Add a value to the initial pipeline state.")
    parser.add_argument(
        dest="pipeline", metavar="pipeline",
        help="YAML file describing the pipeline")
    return parser

def default_log_config():
    return {
        'version': 1,
        'disable_existing_loggers': False,
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
            'config': {'logging': default_log_config()},
            'components': [
                'phyre_engine.component.input.FastaInput',
                'phyre_engine.component.validate.SeqValidator', {
                    'phyre_engine.component.hhsuite.HHBlits': [
                        '/path/to/data/uniclust30',
                        {'iterations': 3},
                    ]
                }
            ]
        }
    }

def init_logging(pipeline):
    config = pipeline.get("config", {})
    LOGGING = "logging"
    DISABLE_LOGS = "disable_existing_loggers"

    if LOGGING not in config:
        config[LOGGING] = default_log_config()
    else:
        # We usually want to keep (and reload) loggers, so we set
        # disable_existing_loggers to False by default. If the user really
        # wants, it can be overridden.
        if DISABLE_LOGS not in config[LOGGING]:
            config[LOGGING][DISABLE_LOGS] = False

    pipeline["config"] = config
    logging.config.dictConfig(config[LOGGING])

def construct_yaml_tuple(self, node):
    # Used to convert sequences from lists to tuples. Only applies to lists
    # without any nested structures.
    seq = self.construct_sequence(node)
    if any(isinstance(e, (list, tuple, dict)) for e in seq):
        return seq
    return tuple(seq)

def main():  # IGNORE:C0111
    '''Command line options.'''

    try:
        parser = arg_parser()
        args = parser.parse_args()
        if args.dump:
            return 0

        # Parse pipeline descriptor from YAML file
        with open(args.pipeline, "r") as yml_in:
            SafeLoader.add_constructor(
                'tag:yaml.org,2002:seq',
                construct_yaml_tuple)
            config = yaml.load(yml_in, SafeLoader)

        # Update starting values if any were supplied on the command line
        if "start" not in config["pipeline"]:
            config["pipeline"]["start"] = {}
        if args.start is not None:
            config["pipeline"]["start"].update(args.start)

        # Set up logging if a logging section was given in the pipeline
        init_logging(config["pipeline"])

        # Load a pipeline
        pipeline = Pipeline.load(config["pipeline"])
        pipeline.run()

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())
