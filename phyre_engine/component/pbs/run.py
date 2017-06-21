#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.component.pbs.run -- Run a sub-pipeline on a worker node

This script is used to run a portion of pipeline on a worker node. It should
only be called automatically by the phyre_engine.component.pbs.qsub.Qsub class.
'''

import sys

from argparse import ArgumentParser
import pickle
import phyre_engine.run
from phyre_engine.component.pbs import qsub
import logging.config

def arg_parser():
    # Setup argument parser
    parser = ArgumentParser(description=__import__('__main__').__doc__)

    parser.add_argument(
        dest="state", metavar="state",
        help="Pickle containing the state of the pipeline")
    return parser


def main():  # IGNORE:C0111
    '''Command line options.'''
    try:
        parser = arg_parser()
        args = parser.parse_args()

        with open(args.state, "r+b") as state_fh:
            state = pickle.load(state_fh)

            # Initialise loggers
            if "logging" in state.pipeline.config:
                logging.config.dictConfig(state.pipeline.config["logging"])
            else:
                logging.config.dictConfig(phyre_engine.run.default_log_config())

            data = state.pipeline.run()
            state_fh.seek(0)
            state_fh.truncate()
            pickle.dump(qsub.CompletedState(data), state_fh)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())
