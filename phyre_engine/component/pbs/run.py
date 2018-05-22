#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.component.pbs.run -- Run a sub-pipeline on a worker node

This script is used to run a portion of pipeline on a worker node. It should
only be called automatically by classes in the phyre_engine.component.pbs.qsub
module.
'''

import sys

from argparse import ArgumentParser
import pickle
import phyre_engine.run
from phyre_engine.component.pbs import qsub
import phyre_engine.logutils
import phyre_engine.tools.yaml as yaml
import logging
import logging.config

ROOT_LOGGER = "phyre_engine.compnent.pbs.run"

def arg_parser():
    # Setup argument parser
    parser = ArgumentParser(description=__import__('__main__').__doc__)

    parser.add_argument(
        dest="pipeline", metavar="pipeline",
        help="JSON file containing the pipeline definition.")
    parser.add_argument(
        dest="state", metavar="state",
        help="Pickle containing the state of the pipeline")
    return parser


def main():  # IGNORE:C0111
    '''Command line options.'''
    try:
        logging.setLoggerClass(phyre_engine.logutils.PhyreEngineLogger)
        parser = arg_parser()
        args = parser.parse_args()
        with open(args.pipeline, "r") as pipeline_fh:
            pipeline_dict = yaml.load(pipeline_fh)
        pipeline = phyre_engine.pipeline.Pipeline.load(pipeline_dict)


        with open(args.state, "r+b") as state_fh:
            state = pickle.load(state_fh)

            # Initialise loggers
            if "logging" in pipeline.config:
                log_config = pipeline.config["logging"]
            else:
                log_config = phyre_engine.run.default_log_config()

            if "disable_existing_loggers" not in log_config:
                log_config["disable_existing_loggers"] = False
            logging.config.dictConfig(log_config)

            pipeline.start = state
            data = pipeline.run()

            data["qsub_complete"] = True

            state_fh.seek(0)
            state_fh.truncate()
            pickle.dump(data, state_fh)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as error:
        logging.getLogger(ROOT_LOGGER).error(
            "Uncaught exception encountered: exiting.",
            exc_info=error)
        raise error

if __name__ == "__main__":
    sys.exit(main())
