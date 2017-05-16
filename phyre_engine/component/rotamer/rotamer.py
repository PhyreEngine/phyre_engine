"""
Component for calculating which rotamer a residue adopts.
"""
from phyre_engine.component import Component
from phyre_engine.tools.rotamer import Rotamer

class CalculateRotamer(Component):
    """
    Find which rotamer is adopted by a given residue.

    This component requires the ``residues`` key to be present. Each element of
    the ``residues`` list is a dictionary containing the keys specified in
    :py:meth:`phyre_engine.component.rotamer.extract.AngleExtractor.run`. This
    component adds the ``rotamer`` key to each element.

    """

    REQUIRED = ["residues"]
    ADDS = []
    REMOVES = []

    def __init__(self, rotamers):
        """
        :param rotamers: Dictionary of rotamer definitions.

        .. seealso::
        Variable `ROTAMERS`
            In modules `tools.rotamer.data.dunbrack` and
            `tools.rotamer.data.molprobity`.
        """
        self.rotamers = rotamers

    def run(self, data):
        """
        Calculate which rotamer each residue in the ``residues`` key of `data`
        adopts.

        The following elements will be added to `data`:

        ``rotamer``
            A :py:class:`phyre_engine.tools.rotamer.Rotamer` object.
        """

        residues = self.get_vals(data)
        for residue in residues:
            rotamer = Rotamer.find(residue["sidechain"], self.rotamers)
            residue["rotamer"] = rotamer
        return data
