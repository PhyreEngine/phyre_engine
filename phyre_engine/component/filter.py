"""Components for filtering keys from the pipeline state."""
from phyre_engine.component.component import Component

class Whitelist(Component):
    """
    Keep only the specified keys in the pipeline state.

    :param list[str] *args: List of keys to keep.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []

    def __init__(self, *args):
        self.whitelist = set(args)

    def run(self, data, config=None, pipeline=None):
        """Apply a whitelist to the pipeline state."""
        new_state = {}
        for key, value in data.items():
            if key in self.whitelist:
                new_state[key] = value
        return new_state
