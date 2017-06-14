"""
Module containing components for parsing a list of residues.

The :py:mod:`phyre_engine.component.rotamer.extract` module contains tools for
extracting angles and side-chain information from  PDB files, but it may
sometimes be useful to parse a pre-generated data file. This module provides
components that read in a file containing residue information and parses it into
a list of residues.
"""
import csv
import re

from Bio.PDB.Residue import Residue

from phyre_engine.component import Component
from phyre_engine.tools.rotamer.rotamer import Sidechain
from phyre_engine.tools.rotamer.data.generic import NUM_CHI_ANGLES


class CsvParser(Component):
    """
    Parse a file in CSV (comma separated value) format.
    """

    REQUIRED = ["residue_csv"]
    ADDS = ["residues"]
    REMOVES = []

    def __init__(self, mapping, final_chi_range, **csv_params):
        """
        Initialise a new CsvPaser object.

        If `mapping` is supplied, it should be a dictionary of strings mapping
        the names of columns in the CSV file to the required keys of each
        element in the ``residues``  keys of the pipeline data.


        :param dict mapping: Mapping of fields in the CSV file to the required
            fields in each residue. The ``residue`` key defines the ID and res_type
            of a :py:mod:`Bio.PDB.Residue`. See the example for details.
        :param dict final_chi_range: Allowed ranges of final χ angles. See
            :ref:`description-of-rotamer-variables` for definition.
        :param **csv_params: Parameters to be parsed to :py:func:`csv.reader`.

        Consider a CSV file that looks like the following:

        .. code-block:: none
            AA RES_ID PHI PSI CHI1 CHI2 CHI3 CHI4
            SER 1 -60 60 30.12 0 0 0
            VAL 2 -13 30 10.69 0 0 0

        You would parse the following parameters to the constructor:

        >>> parser = CsvParser(mapping={
                    "resiue": {"id": (None, "RES_ID", None), "resname": "AA"},
                    "phi": "PHI",
                    "psi": "PSI",
                    "chi": ("CHI1", "CHI2", "CHI3", "CHI4")
                }, delimiter=' ')

        The `None` values in the residue ID indicate that there are no columns
        giving the residue hetero-flag or insertion code. These will be replaced
        with sensible default values. Some extra care is taken when parsing the
        residue ID: if the residue ID is not purely numeric, the non-numeric
        parts will be treated as the insertion code.

        The ``residue`` key of `mapping` is mandatory, but phi, psi and chi are
        optional and default to the following:

        .. code-block:: python
            mapping = {
                "phi": "phi",
                "psi": "psi",
                "chi": ("chi1", "chi2", "chi3", "chi4")
            }

        Any extra keys in mapping are also added to the residue: e.g.
        ``{"source": "PDB_FILE"}``.

        :raise KeyError: If the ``residue`` key is not present in `mapping` or
            is invalid.

        .. seealso::

        Function :py:func:`csv.reader`
            For the allowed parameters in `csv_params`.

        Function :py:func:`phyre_engine.component.rotamer.extract.AngleExtractor.run`
            For the keys required for each residue element.

        Section :ref:`description-of-rotamer-variables`
            For the definition of the `final_chi_range` parameter.
        """
        self.mapping = mapping
        self.final_chi_range = final_chi_range
        self.csv_params = csv_params

        # Raise exception if the "residue" key is not present.
        _ = mapping["residue"]
        _ = mapping["residue"]["id"]
        _ = mapping["residue"]["resname"]

        # Default keys for mapping
        if "phi" not in self.mapping:
            self.mapping["phi"] = "phi"
        if "psi" not in self.mapping:
            self.mapping["psi"] = "psi"
        if "chi" not in self.mapping:
            self.mapping["chi"] = ("chi1", "chi2", "chi3", "chi4")

    def run(self, data, config=None, pipeline=None):
        """
        Parse a CSV file and add the ``residues`` key to `data`.

        See :py:func:`phyre_engine.component.rotamer.extract.AngleExtractor.run`
        for the format of the ``residues`` list that will be added by this
        method.
        """
        residue_csv = self.get_vals(data)

        residues = []
        with open(residue_csv, "r", newline="") as csv_fh:
            reader = csv.DictReader(csv_fh, **self.csv_params)
            for values in reader:

                # Parse out residue ID
                res_id = [' ', None, ' ']
                for i, id_elem in enumerate(self.mapping["residue"]["id"]):
                    if id_elem is not None and values[id_elem] is not None:
                        res_id[i] = values[id_elem]
                try:
                    res_id[1] = int(res_id[1])
                except ValueError as err:
                    # Try and split numeric and non-numeric parts of the ID into
                    # residue number and insertion code.
                    numeric_match = re.search(r"(\d+)", res_id[1])
                    non_numeric_match = re.search(r"(\D+)", res_id[1])

                    # Bail if the matches fail (i.e. there are junk characters)
                    if numeric_match is None or non_numeric_match is None:
                        raise err

                    res_id[1] = int(numeric_match.group(0))
                    res_id[2] = non_numeric_match.group(0)

                res_id = tuple(res_id)

                # Get residue name
                res_name = values[self.mapping["residue"]["resname"]]

                # Build Bio.PDB.Residue object
                residue = Residue(res_id, res_name, 1)

                # Get phi, psi, chi angles
                num_chi = NUM_CHI_ANGLES[res_name]
                phi = float(values[self.mapping["phi"]])
                psi = float(values[self.mapping["psi"]])

                chi = [float(values[name])
                       for name in self.mapping["chi"][0:num_chi]]

                # Force final chi angle into correct range
                if (res_name in self.final_chi_range
                    and chi[-1] not in self.final_chi_range[res_name]):
                    chi[-1] = (chi[-1] + 180) % 360

                residue = {
                    "torsion": (phi, psi),
                    "sidechain": Sidechain(res_name, tuple(chi)),
                    "residue": residue
                }

                # Parse remaining mapped columns
                already_parsed = ("phi", "psi", "chi", "residue")
                for to_key, from_key in self.mapping.items():
                    if to_key not in already_parsed:
                        residue[to_key] = values[from_key]

                residues.append(residue)
        data["residues"] = residues
        return data
