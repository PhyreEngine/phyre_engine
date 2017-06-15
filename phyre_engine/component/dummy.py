from phyre_engine.component.component import Component
import pprint
import logging

log = lambda: logging.getLogger(__name__)

class Dummy(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, prefix):
        self.prefix = prefix

    def run(self, data, config=None, pipeline=None):
        log().warning("Dummy warning")
        for ln in pprint.pformat(data).split("\n"):
            print("{}: {}".format(self.prefix, ln))
        return data
