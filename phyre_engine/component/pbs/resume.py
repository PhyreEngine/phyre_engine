#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.component.pbs.resume -- Resume a pipeline after workers finish

This script is used to resume a pipeline after workers have finished processing
their sub-pipelines.
'''

import sys

from argparse import ArgumentParser
import pickle
from pathlib import Path
from phyre_engine.component.pbs import qsub
import logging.config

log = lambda: logging.getLogger(__name__)

def arg_parser():
    # Setup argument parser
    parser = ArgumentParser(description=__import__('__main__').__doc__)

    parser.add_argument(
        dest="storage_dir", metavar="dir", type=str,
        help="Directory containing pipeline state.")
    parser.add_argument(
        dest="num_jobs", metavar="num_jobs", type=int,
        help="Number of sub-pipeline states.")
    parser.add_argument(
        dest="join_var", metavar="join_var", type=str,
        help="Join sub-pipe data on this variable.")
    return parser


def main():  # IGNORE:C0111
    '''Command line options.'''
    try:
        parser = arg_parser()
        args = parser.parse_args()
        storage_dir = Path(args.storage_dir)

        # Get the pipeline to resume
        with (storage_dir / qsub.RESUME_PIPELINE).open("rb") as pipe_in:
            pipe_state = pickle.load(pipe_in)
            # Initialise loggers
            if "logging" in pipe_state.config:
                logging.config.dictConfig(pipe_state.config["logging"])


        # Generate the new state
        state = None
        for i in range(0, args.num_jobs):
            state_file = storage_dir / qsub.WORKER_STATE.format(i)
            try:
                with state_file.open("rb") as state_in:
                    sub_state = pickle.load(state_in)
                    if not isinstance(sub_state, qsub.CompletedState):
                        raise ValueError("No data in {}".format(state_file))

                    if state is None:
                        state = sub_state.data
                    else:
                        sub_slice = sub_state.data[args.join_var]
                        state[args.join_var].extend(sub_slice)
            except Exception as err:
                log().error(
                    "Error reading data from worker %d (%s)",
                    i, str(state_file),
                    exc_info=err)
                raise err

        pipe_state.pipeline.start = state
        pipe_state.pipeline.run(pipe_state.pipeline_index + 1)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())
