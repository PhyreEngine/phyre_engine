"""
Rotamer definitions according to Dunbrack.

The rotameric amino acids (ARG, CYS, ILE, LEU, LYS, MET, PRO, SER, THR, and VAL)
are treated as usual. Each chi angle is split into three regions:

    g+: 0  ≤ χ < 120
    t:  120 ≤ χ < 240
    g-: 240 ≤ χ < 360

In the parlance of Dunbrack's rotamer library, these are rotamers 1, 2 and 3
respectively.

Eight of the standard amino acids have non-rotameric degrees of freedom in their
final chi angle: ASN, ASP, GLN, GLU, HIS, TRP, PHE and TYR. The final rotamers
of ASN, ASP, GLN and GLU are defined as (from 1--3):

    g+: 30  ≤ χ < 90
    t:  330 ≤ χ < 360 and 0 ≤ χ < 30
    g-: 270 ≤ χ < 330

The final rotamers of PHE, TYR and HIS are defined as (from 1--2):

    g:  30  ≤ χ < 150
    t:  330 ≤ χ < 360 and 0 ≤ χ < 30

The final rotamers of TRP are:

    g+: 180  ≤ χ < 300
    t:  300 ≤ χ < 360 and 0 ≤ χ < 60
    g-: 60 ≤ χ < 180

Finally, the rotamers of PRO are defined as defined as:

    g+ (endo): 0  ≤ χ < 90
    g- (exo): 270 ≤ χ < 360

Note that amino acids with rings or O-C=O terminating units are only determined
on an interval of 180°; that is, they are considered to be equivalent when
flipped. Following the original paper, we do the same for HIS, ASN and GLN.

This module exposes the following variables:

:var ROTAMERS: Dictionary indexed by three-letter amino acid codes.

:var FINAL_CHI_RANGE: Set of amino acids for which the final rotating unit
    is only determined within 0--180 degrees. This module uses the rotamer
    definitions of Dunbrack's 1997 paper, so ASN, ASP, GLN, GLU, PHE, TYR, HIS,
    PRO and TRP are defined only over 180 degrees.

.. seealso::

    Section :ref:`description-of-rotamer-variables`
        For a description of the :py:data:`ROTAMERS` and
        :py:data:`FINAL_CHI_RANGE` variables.
"""

from ..angle_range import AngleRange
from .generic import NUM_CHI_ANGLES

FINAL_CHI_RANGE = {}
for _aa in ("ASN", "ASP", "GLN", "GLU", "PRO"):
    FINAL_CHI_RANGE[_aa] = AngleRange((0, 90), (270, 360))
for _aa in ("PHE", "TYR", "HIS"):
    FINAL_CHI_RANGE[_aa] = AngleRange((0, 150), (330, 360))

# Mapping of Dunbrack's IDs (1-3) for rotameric chi angles to the corresponding
# ranges
ROTAMERIC_ANGLES = {
    1: AngleRange((0, 120)),
    2: AngleRange((120, 240)),
    3: AngleRange((240, 360))
}
# Some of the final rotamers have different definitions:
ASN_ASP_GLN_GLU_FINAL = {
    1: AngleRange((30, 90)),
    2: AngleRange((330, 360), (0, 30)),
    3: AngleRange((270, 330)),
}
PHE_TYR_HIS_FINAL = {
    1: AngleRange((30, 150)),
    2: AngleRange((330, 360), (0, 30))
}
TRP_FINAL = {
    1: AngleRange((180, 300)),
    2: AngleRange((300, 360), (0, 60)),
    3: AngleRange((60, 180))
}
PRO_FINAL = {
    1: AngleRange((0, 90)),
    2: AngleRange((270, 360))
}

ALTERED_LAST_ROTAMERS = {
    "ASN": ASN_ASP_GLN_GLU_FINAL,
    "ASP": ASN_ASP_GLN_GLU_FINAL,
    "GLN": ASN_ASP_GLN_GLU_FINAL,
    "GLU": ASN_ASP_GLN_GLU_FINAL,
    "PHE": PHE_TYR_HIS_FINAL,
    "TYR": PHE_TYR_HIS_FINAL,
    "HIS": PHE_TYR_HIS_FINAL,
    "TRP": TRP_FINAL,
    "PRO": PRO_FINAL
}

ROTAMERS = {aa: {} for aa, n in NUM_CHI_ANGLES.items() if n > 0}

# Recursive function to define all rotamers.
# The "key" parameter is a list indicating the rotamer of interest (e.g. (1,1))
# The "value" parameter is a list of AngleRanges defining each rotamer.
def init_rotamers(aa, key=None, value=None, index=0):
    # Initialise an empty list of the correct length
    if key is None:
        key = [None] * NUM_CHI_ANGLES[aa]
    if value is None:
        value = [None] * NUM_CHI_ANGLES[aa]

    # Default to using the standard rotameric angles
    chi_list = ROTAMERIC_ANGLES
    if index >= NUM_CHI_ANGLES[aa]:
        # We're done recursing
        # When we exceed the number of chi angles, assign the list of IDs and
        # angle ranges to the ROTAMERS dict.
        ROTAMERS[aa][tuple(key)] = tuple(value)
        return
    elif index == NUM_CHI_ANGLES[aa] - 1:
        # Only the last chi angle can be non-rotameric, so here we check whether
        # we need to use an altered rotamer definition
        if aa in ALTERED_LAST_ROTAMERS:
            chi_list = ALTERED_LAST_ROTAMERS[aa]

    # Set the current key and value and recurse
    for rotamer, ang_range in chi_list.items():
        key[index] = rotamer
        value[index] = ang_range
        init_rotamers(aa, key, value, index + 1)

for _aa in ROTAMERS.keys():
    init_rotamers(_aa)
