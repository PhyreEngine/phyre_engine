import contextlib
import io
import logging
import os


class PhyreEngineLogger(logging.Logger):
    """
    Custom logger used by PhyreEngine to log events.

    This logger adds the following attributes to each log record:

    ``hostname``
        Name of the current host.
    """
    def __init__(self, name):
        super().__init__(name)
        self.addFilter(HostnameFilter())


class HostnameFilter(logging.Filter):
    """
    This filter injects hostname information into each log record. The hostname
    can be accessed through the "hostname" attribute."
    """

    def filter(self, record):
        record.hostname = os.uname().nodename
        return True

@contextlib.contextmanager
def capture_log(logger, level=logging.DEBUG):
    log_capture_buf = io.StringIO()
    log_capture_handler = logging.StreamHandler(log_capture_buf)
    log_capture_handler.setLevel(level)

    log_capture_formatter = logging.Formatter('%(levelname)s - %(message)s')
    log_capture_handler.setFormatter(log_capture_formatter)

    try:
        logger.addHandler(log_capture_handler)
        yield log_capture_buf
    finally:
        logger.removeHandler(log_capture_handler)

def name(obj):
    """
    Standardised logger name for an object.

    This is just the ``__module__`` and ``__qualname__`` separated by a dot.
    """
    return ".".join((type(obj).__module__, type(obj).__qualname__))
