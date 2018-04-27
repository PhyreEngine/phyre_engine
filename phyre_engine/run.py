#!/usr/bin/env python
# encoding: utf-8
'''
phyre_engine.run -- Run a PhyreEngine pipeline

PhyreEngine is a framework for building and running common structural biology
pipelines. This script is responsible for reading, setting up and running a
PhyreEngine pipeline. Unless you're an expert, this script is probably all that
you need to run a pipeline.
'''

DETAILS = """\
The pipeline configuration is pre-processed with a template engine before the
pipeline is loaded. Templates use Python string formatting and are specified
with the "!template" tag. Items in the pipeline configuration are exposed to
the templates, along with the current environment in "ENV".

For example, user-specific tools may be configured like this:

    pipeline:
      config:
        PREFIX !template '/data/{{ENV[USER]}}/conda/env/phyreengine'
        dssp:
          bin_dir: !template '{{PREFIX}}/bin'
      disopred:
          bin_dir: !template '{{PREFIX}}/bin'
          data_dir: !template '{{PREFIX}}/share/disopred/data'
          dso_lib_dir: !template '{{PREFIX}}/share/disopred/dso_lib'

This would set "PREFIX" to "/data/$USER/conda/env/phyreengine", where "$USER"
is an environment variable. The remaining items then use "PREFIX" to configure
dssp and disopred.

You can also set a default pipeline configuration in the file
"{default_config}".
This is merged with the pipeline config given in the pipeline file, with the
pipeline file taking precedence. Templates may also be used in the default
configuration file.
"""

import logging
import logging.config
import os
import pathlib
import sys

import appdirs

import argparse
import phyre_engine.logutils
import phyre_engine.pipeline
from phyre_engine.tools.util import apply_dotted_key
import phyre_engine.tools.yaml as yaml

APP_SHORTAUTHOR = "imperial_college"
APP_SHORTNAME = "phyreengine"

class DumpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        yaml.dump(dummy_pipeline(), sys.stdout, default_flow_style=False)
        sys.exit(0)

class StoreStartingValue(argparse.Action):
    """
    Parse a starting value for the pipeline state.

    For example, assume that this action is associated with the ``--start``
    option and the ``start`` ``dest``. The user passes the following options:

    .. code-block:: none

        --start foo:bar --start baz:qux1 --start baz:qux2

    This will result in the namespace returned by the argument parser having a
    ``start`` (set from the ``dest`` parameter) attribute that looks like the
    following:

    .. code-block:: python

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
    parser = argparse.ArgumentParser(
        description=__import__('__main__').__doc__,
        epilog=DETAILS.format(
            default_config=default_config_file()),
        formatter_class=argparse.RawDescriptionHelpFormatter)

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
        "-c", "--config", dest="config", action=StoreStartingValue, nargs=1,
        default={}, help="Modify pipeline configuration.")
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
                'format': '{levelname} : {asctime} : {hostname} : {name} : {message}',
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


def init_logging(config):
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

    logging.setLoggerClass(phyre_engine.logutils.PhyreEngineLogger)
    logging.config.dictConfig(config[LOGGING])

def resolve_yml(document, getter, *extra_fields):
    """
    Process the :py:class:`phyre_engine.tools.yaml.ImplicitDocument` `document`
    until all templates are resolved.

    The `*extra_fields` arguments must be `dict`s, and will be merged together
    along the results of the callback `getter`, which is passed the document
    as an argument.
    """
    while isinstance(document, yaml.ImplicitDocument):
        fields = {}
        for extra in extra_fields:
            fields.update(extra)
        fields.update(getter(document.unresolved))
        document = document.resolve(fields, allow_unresolved=True)
    return document


def default_config_file():
    """Path of the default configuration file."""
    dirs = appdirs.AppDirs(APP_SHORTNAME, APP_SHORTAUTHOR)
    config_dir = pathlib.Path(dirs.user_config_dir)
    return config_dir / "config.yml"


def load_default_config():
    """
    Load default configuration file, if any. This is pre-processed in the
    same way as the main pipeline configuration.
    """
    config_file = default_config_file()
    if not config_file.exists():
        return {}
    else:
        with config_file.open("r") as conf_in:
            return yaml.load(conf_in, yaml.ImplicitLoader)

def pipeline_description(pipeline_file):
    """
    Load pipeline description file and pre-process the config.
    """

    pipe_config = phyre_engine.pipeline.PipelineConfig(
        resolve_yml(
            load_default_config(),
            lambda doc: doc,
            {"ENV": os.environ}))

    # Parse pipeline descriptor from YAML file
    with open(pipeline_file, "r") as yml_in:
        pipeline_desc = yaml.load(yml_in, yaml.ImplicitLoader)

    # Merge default configuration with the current unresolved config
    pipe_config = pipe_config.merge_params(
        phyre_engine.pipeline.PipelineConfig(
            pipeline_desc.unresolved["pipeline"].get("config", {})))
    pipeline_desc.unresolved["pipeline"]["config"] = pipe_config

    pipeline_desc = resolve_yml(
        pipeline_desc,
        lambda doc: doc["pipeline"]["config"],
        {"ENV": os.environ})

def main():  # IGNORE:C0111
    '''Command line options.'''

    # Before doing anything, tell the loggers to capture warnings.
    logging.captureWarnings(True)
    try:
        parser = arg_parser()
        args = parser.parse_args()
        if args.dump:
            return 0

        pipe_config = phyre_engine.pipeline.PipelineConfig(
            resolve_yml(
                load_default_config(),
                lambda doc: doc,
                {"ENV": os.environ}))

        # Parse pipeline descriptor from YAML file
        with open(args.pipeline, "r") as yml_in:
            pipeline_desc = yaml.load(yml_in, yaml.ImplicitLoader)

        # Merge default configuration with the current unresolved config
        pipe_config = pipe_config.merge_params(
            phyre_engine.pipeline.PipelineConfig(
                pipeline_desc.unresolved["pipeline"].get("config", {})))
        pipeline_desc.unresolved["pipeline"]["config"] = pipe_config

        pipeline_desc = resolve_yml(
            pipeline_desc,
            lambda doc: doc["pipeline"]["config"],
            {"ENV": os.environ})

        # Set default start values"
        pipeline_desc["pipeline"].setdefault("start", {})
        if args.start is not None:
            pipeline_desc["pipeline"]["start"].update(args.start)

        # Set up logging if a logging section was given in the pipeline
        init_logging(pipeline_desc["pipeline"]["config"])

        # Add extra configuration keys to the pipeline config
        for dotted_key, value in args.config.items():
            apply_dotted_key(pipeline_desc["pipeline"]["config"],
                             dotted_key, value)

        # Load a pipeline
        pipeline = phyre_engine.pipeline.Pipeline.load(
            pipeline_desc["pipeline"])
        pipeline.run()

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as error:
        logging.error(
            "Uncaught exception encountered: exiting.",
            exc_info=error)
        raise error

if __name__ == "__main__":
    sys.exit(main())
