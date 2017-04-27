"""
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

:var SYMMETRIC_FINAL_CHI: Tuple of amino acids for which the final rotating unit
    is a ring or ``O-C=O``. In these cases, it doesn't make any sense to
    differentiate between angles in the 0--180° range and those in the 180--360°
    range, so the angles of these residues should be taken modulo 180°.

>>> CHI_ATOMS["VAL"][0]
("N", "CA" , "CB" , "CG1"
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

SYMMETRIC_FINAL_CHI = ("PHE", "TYR", "ASP", "GLU")

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
