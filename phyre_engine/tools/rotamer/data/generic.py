r"""
Generic constants for defining rotamers. This module exposes the following
variables:

:var AMINO_ACIDS: Tuple of three-letter amino acid names.
:vartype AMINO_ACIDS: Tuple of strings.

:var NUM_CHI_ANGLES: Dictionary containing the number of χ angles for each
    residue type. The keys are three-letter residue names.
:vartype NUM_CHI_ANGLES: Dictionary of amino acids to integers.

:var CHI_ATOMS: Dictionary containing tuples of atoms defining each χ angle. The
    keys to the dictionary are three-letter amino acid names, and the
    corresponding values are a list of tuples. The tuples each contain four
    elements, corresponding to the atom names defining that χ angle.
:vartype CHI_ATOMS: Dictionary of lists of tuples.

>>> CHI_ATOMS["VAL"][0]
("N", "CA" , "CB" , "CG1")

.. _description-of-rotamer-variables:

Rotamer definitions
-------------------

Different rotamer libraries use different rotamer definitions, and even
different definitions of side-chain angles. For example, Dunbrack's 1997 rotamer
library defines the final rotamer of ASN, ASP, GLN and GLU only on the range
-90--90°, while the MolProbity rotamer library considers ASN and GLN to be
defined on the full range of -180--180° but defines the ranges for ASP and GLU
on the ranges 0--180°.

This module supplies generic information regarding rotamer libraries: it is
unlikely that any implementation will need to alter the number of χ angles on
each amino acid, or the atoms over which each angle is defined. Rotamer
libraries should be defined in separate modules, and should contain variables
defining the angular ranges over which each rotamer is defined and the ranges
over which each non-standard final χ atom is defined. These variables can be
named arbitrarily, as they must be explicitly passed to any methods making use
of them, but the modules included with this distribution name them ``ROTAMERS``
and ``FINAL_CHI_RANGE``. In this discussion, I will use this convention to refer
to these variables.

The ``ROTAMERS`` variable
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ROTAMERS`` variable is a two-level dictionary indexed first by the
(three-letter) amino acid name, and then by the rotamer name. The value should
be a tuple of :py:class:`phyre_engine.tools.rotamer.AngleRange`\ s defining the
ranges of χ angles over which a given rotamer is defined. The rotamer names may
be any hashable value. For example, here is a rotamer for arginine:

.. code-block:: python

    ROTAMERS["ARG"]["mmm-85"] = (
        AngleRange((240, 360)),
        AngleRange((240, 360)),
        AngleRange((240, 360)),
        AngleRange((240, 360)),
    )

Each angle range in the tuple describes the allowed range for the corresponding
χ angles. In this case, each χ angle is allowed to range between 240 and 360°.
For amino acids with fewer χ angles, fewer elements need be included in the
tuple.

The ``FINAL_CHI_RANGE`` variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some amino acids have a final rotating unit on their side-chain that is
(roughly) symmetric. In the case of PHE and TYR, this is a ring; for ASP and
GLU, this is a carboxyl group. These amino acids are typically only defined in a
180° range. Some rotamer libraries take this further, and also limit HIS, ASN
and GLN to a 180° angular range.

This variable must be a dictionary indexed by (three-letter) amino acid names,
each value of which is an :py:class:`phyre_engine.tools.rotamer.AngleRange`.
These definitions may be used when calculating side-chain angles to flip the
final χ angle by 180°. For example:

.. code-block:: python

    FINAL_CHI_RANGE["PHE"] = AngleRange((0, 180))

"""

AMINO_ACIDS = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)

NUM_CHI_ANGLES = {
    "GLY": 0, "ALA": 0, "VAL": 1, "LEU": 2, "ILE": 2, "PRO": 1, "PHE": 2,
    "TYR": 2, "TRP": 2, "SER": 1, "THR": 1, "CYS": 1, "MET": 3, "MSE": 3,
    "LYS": 4, "HIS": 2, "ARG": 4, "ASP": 2, "ASN": 2, "GLN": 3, "GLU": 3,
}


CHI_ATOMS = {aa: [None] * n for aa, n in NUM_CHI_ANGLES.items() if n > 0}
CHI_ATOMS["VAL"][0] = ("N", "CA" , "CB" , "CG1")
CHI_ATOMS["LEU"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["LEU"][1] = ("CA", "CB" , "CG" , "CD1")
CHI_ATOMS["ILE"][0] = ("N", "CA" , "CB" , "CG1")
CHI_ATOMS["ILE"][1] = ("CA", "CB" , "CG1", "CD1")
CHI_ATOMS["PRO"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["PHE"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["PHE"][1] = ("CA", "CB" , "CG" , "CD1")
CHI_ATOMS["TYR"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["TYR"][1] = ("CA", "CB" , "CG" , "CD1")
CHI_ATOMS["TRP"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["TRP"][1] = ("CA", "CB" , "CG" , "CD1")
CHI_ATOMS["SER"][0] = ("N", "CA" , "CB" , "OG")
CHI_ATOMS["THR"][0] = ("N", "CA" , "CB" , "OG1")
CHI_ATOMS["CYS"][0] = ("N", "CA" , "CB" , "SG")
CHI_ATOMS["MET"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["MET"][1] = ("CA", "CB" , "CG" , "SD")
CHI_ATOMS["MET"][2] = ("CB", "CG" , "SD" , "CE")
CHI_ATOMS["MSE"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["MSE"][1] = ("CA", "CB" , "CG" , "SE")
CHI_ATOMS["MSE"][2] = ("CB", "CG" , "SE"  , "CE")
CHI_ATOMS["LYS"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["LYS"][1] = ("CA", "CB" , "CG" , "CD")
CHI_ATOMS["LYS"][2] = ("CB", "CG" , "CD" , "CE")
CHI_ATOMS["LYS"][3] = ("CG", "CD" , "CE" , "NZ")
CHI_ATOMS["HIS"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["HIS"][1] = ("CA", "CB" , "CG" , "ND1")
CHI_ATOMS["ARG"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["ARG"][1] = ("CA", "CB" , "CG" , "CD")
CHI_ATOMS["ARG"][2] = ("CB", "CG" , "CD" , "NE")
CHI_ATOMS["ARG"][3] = ("CG", "CD" , "NE" , "CZ")
CHI_ATOMS["ASP"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["ASP"][1] = ("CA", "CB" , "CG" , "OD1")
CHI_ATOMS["ASN"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["ASN"][1] = ("CA", "CB" , "CG" , "OD1")
CHI_ATOMS["GLN"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["GLN"][1] = ("CA", "CB" , "CG" , "CD")
CHI_ATOMS["GLN"][2] = ("CB", "CG" , "CD" , "OE1")
CHI_ATOMS["GLU"][0] = ("N", "CA" , "CB" , "CG")
CHI_ATOMS["GLU"][1] = ("CA", "CB" , "CG" , "CD")
CHI_ATOMS["GLU"][2] = ("CB", "CG" , "CD" , "OE1")
