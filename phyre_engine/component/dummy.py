from phyre_engine.component.component import Component
import time
import pickle
import os

class Dummy(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, prefix, add=None):
        self.prefix = prefix
        self.add = add if add is not None else []

    def run(self, data, config=None, pipeline=None):
        self.logger.warning("Dummy warning: %s", self.prefix)
        data[self.prefix] = self.add
        return data

class Multiply(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, var, by):
        self.var = var
        self.by = by

    def run(self, data, config=None, pipeline=None):
        data[self.var] = [x * self.by for x in data[self.var]]
        return data

class Sleep(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, seconds):
        self.seconds = seconds

    def run(self, data, config=None, pipeline=None):
        time.sleep(self.seconds)
        return data

class Dump(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, file):
        self.file = file

    def run(self, data, config=None, pipeline=None):
        with open(self.file, "wb") as file_out:
            pickle.dump(data, file_out)
        return data

class Hostname(Component):
    ADDS = ["hostname"]
    REMOVES = []
    REQUIRED = []

    def run(self, data, config=None, pipeline=None):
        data["hostname"] = os.uname().nodename
        return data

class LogTest(Component):
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def run(self, data, config=None, pipeline=None):
        self.logger.debug("Debug message")
        self.logger.info("Info message")
        self.logger.warn("Warning message")
        self.logger.error("Error message")
        self.logger.critical("Critical error message")
        return data
