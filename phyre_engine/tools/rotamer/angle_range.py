class AngleRange:
    """
    Representation of an angular range. This only needs to exist because an
    angular range can fall across the zero, complicating checking using python's
    `range()` or manually checking a tuple. Any angle greater than 360° is
    calculated modulo 360.

    >>> ang_range = AngleRange((0, 15), (345, 360))
    >>> 0 in ang_range
    True
    >>> 15 in ang_range
    False
    >>> 45.5 in ang_range
    False
    >>> 350 in ang_range
    True
    >>> 360 in ang_range
    True
    >>> 361 in ang_range
    True
    """

    def __init__(self, *args):
        """
        Initialise a new angular range given a list of tuples defining start and
        end angles. The first angle of each tuple is considered to be included
        in the range, and the second angle excluded.

        :param *args: List of tuples defining ranges of angles.
        """
        self.ranges = args

    def __contains__(self, angle):
        """
        Check whether the given angle is contained within this angle range.

        :param angle: Angle to check.
        """
        if angle >= 360:
            angle = angle % 360

        for ang_range in self.ranges:
            if ang_range[0] <= angle < ang_range[1]:
                return True
        return False
