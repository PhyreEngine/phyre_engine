"""
Rotamer definitions according to MolProbity (and therefore from the Penultimate
Rotamer Library and (Son of) Penutlimate Rotamer Library.

This module exposes the `ROTAMERS` variable, which is a dictionary indexed by
three-letter amino acid codes. Each key points to another dictionary, indexed by
the name of the rotamer according to MolProbity nomenclature. The values of this
dictionary are a tuple of `AngleRange`s indicating the actual range of each chi
angle for this rotamer.

:var ROTAMERS: Dictionary indexed by three-letter amino acid codes.

:var SYMMETRIC_FINAL_CHI: Set of amino acids for which the final rotating unit
    is a ring or ``O-C=O``. In these cases, it doesn't make any sense to
    differentiate between angles in the 0--180° range and those in the 180--360°
    range, so the angles of these residues should be taken modulo 180°.

>>> ROTAMERS["ARG"]["ppp80"]
"""

from ..angle_range import AngleRange
from .generic import NUM_CHI_ANGLES

# These are considered to be uniquely determined only within 0-180 degrees.
SYMMETRIC_FINAL_CHI = set(("PHE", "TYR", "ASP", "GLU"))

ROTAMERS = {aa: {} for aa, n in NUM_CHI_ANGLES.items() if n > 0}
ROTAMERS["THR"]["p"] = (
    AngleRange((0 , 120)),
)
ROTAMERS["THR"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["THR"]["m"] = (
    AngleRange((240, 360)),
)

ROTAMERS["VAL"]["p"] = (
    AngleRange((0 , 120)),
)
ROTAMERS["VAL"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["VAL"]["m"] = (
    AngleRange((240, 360)),
)

ROTAMERS["SER"]["p"] = (
    AngleRange((0 , 120)),
)
ROTAMERS["SER"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["SER"]["m"] = (
    AngleRange((240, 360)),
)

ROTAMERS["CYS"]["p"] = (
    AngleRange((0 , 120)),
)
ROTAMERS["CYS"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["CYS"]["m"] = (
    AngleRange((240, 360)),
)

ROTAMERS["PRO"]["Cg_endo"] = (
    AngleRange((0, 180)),
)
ROTAMERS["PRO"]["Cg_exo"] = (
    AngleRange((180, 360)),
)

ROTAMERS["PHE"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["PHE"]["t80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180)),
)
ROTAMERS["PHE"]["m-85"] = (
    AngleRange((240, 360)),
    AngleRange((35, 150)),
)
ROTAMERS["PHE"]["m-30"] = (
    AngleRange((240, 360)),
    AngleRange((0, 35), (150, 180)),
)

ROTAMERS["TYR"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["TYR"]["t80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180)),
)
ROTAMERS["TYR"]["m-85"] = (
    AngleRange((240, 360)),
    AngleRange((35, 150)),
)
ROTAMERS["TYR"]["m-30"] = (
    AngleRange((240, 360)),
    AngleRange((0, 35), (150, 180)),
)

ROTAMERS["TRP"]["p-90"] = (
    AngleRange((0, 120)),
    AngleRange((180, 360)),
)
ROTAMERS["TRP"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["TRP"]["t-105"] = (
    AngleRange((120, 240)),
    AngleRange((180, 305)),
)
ROTAMERS["TRP"]["t90"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180), (305, 360)),
)
ROTAMERS["TRP"]["m-90"] = (
    AngleRange((240, 360)),
    AngleRange((180, 305)),
)
ROTAMERS["TRP"]["m0"] = (
    AngleRange((240, 360)),
    AngleRange((0, 45), (305, 360)),
)
ROTAMERS["TRP"]["m95"] = (
    AngleRange((240, 360)),
    AngleRange((45, 180)),
)

ROTAMERS["HIS"]["p-80"] = (
    AngleRange((0, 120)),
    AngleRange((180, 360)),
)
ROTAMERS["HIS"]["p80"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["HIS"]["t-160"] = (
    AngleRange((120, 240)),
    AngleRange((130, 225)),
)
ROTAMERS["HIS"]["t-80"] = (
    AngleRange((120, 240)),
    AngleRange((225, 360)),
)
ROTAMERS["HIS"]["t60"] = (
    AngleRange((120, 240)),
    AngleRange((0, 130)),
)
ROTAMERS["HIS"]["m-70"] = (
    AngleRange((240, 360)),
    AngleRange((0, 20), (225, 360)),
)
ROTAMERS["HIS"]["m170"] = (
    AngleRange((240, 360)),
    AngleRange((130, 225)),
)
ROTAMERS["HIS"]["m80"] = (
    AngleRange((240, 360)),
    AngleRange((20, 130)),
)

ROTAMERS["LEU"]["pp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LEU"]["pt?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LEU"]["tp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LEU"]["tt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LEU"]["tm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LEU"]["mp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LEU"]["mt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LEU"]["mm?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)

ROTAMERS["ILE"]["pp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["ILE"]["pt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["ILE"]["tp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ILE"]["tt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["ILE"]["tm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["ILE"]["mp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["ILE"]["mt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["ILE"]["mm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)

ROTAMERS["ASN"]["p-10"] = (
    AngleRange((0, 120)),
    AngleRange((180, 360)),
)
ROTAMERS["ASN"]["p30"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["ASN"]["t-20"] = (
    AngleRange((120, 240)),
    AngleRange((0, 10), (180, 360)),
)
ROTAMERS["ASN"]["t30"] = (
    AngleRange((120, 240)),
    AngleRange((10, 180)),
)
ROTAMERS["ASN"]["m-20"] = (
    AngleRange((240, 360)),
    AngleRange((300, 360), (0, 40)),
)
ROTAMERS["ASN"]["m-80"] = (
    AngleRange((240, 360)),
    AngleRange((200, 300)),
)
ROTAMERS["ASN"]["m120"] = (
    AngleRange((240, 360)),
    AngleRange((40, 200)),
)

ROTAMERS["ASP"]["p-10"] = (
    AngleRange((0, 120)),
    AngleRange((90, 180)),
)
ROTAMERS["ASP"]["p30"] = (
    AngleRange((0, 120)),
    AngleRange((0, 90)),
)
ROTAMERS["ASP"]["t0"] = (
    AngleRange((120, 240)),
    AngleRange((0, 45), (120, 180)),
)
ROTAMERS["ASP"]["t70"] = (
    AngleRange((120, 240)),
    AngleRange((45, 120)),
)
ROTAMERS["ASP"]["m-20"] = (
    AngleRange((240, 360)),
    AngleRange((0, 180)),
)

ROTAMERS["GLN"]["pt20"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["pm0"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["pp0?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["tp-100"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((150, 300)),
)
ROTAMERS["GLN"]["tp60"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 150), (300, 360)),
)
ROTAMERS["GLN"]["tt0"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["tm0?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["mp0"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["mt-30"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["GLN"]["mm-40"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 60), (210, 360)),
)
ROTAMERS["GLN"]["mm100"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((60, 210)),
)

ROTAMERS["GLU"]["pp20?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["pt-20"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["pm0"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["tp10"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["tt0"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["tm-20"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["mp0"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["mt-10"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 180)),
)
ROTAMERS["GLU"]["mm-40"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 180)),
)

ROTAMERS["MET"]["ppp?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["MET"]["ptp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["MET"]["ptt?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["MET"]["ptm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["MET"]["pmm?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["MET"]["tpp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120), (330, 360)),
)
ROTAMERS["MET"]["tpt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 330)),
)
ROTAMERS["MET"]["ttp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["MET"]["ttt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["MET"]["ttm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["MET"]["tmt?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["MET"]["tmm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["MET"]["mpp?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["MET"]["mpt?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["MET"]["mtp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["MET"]["mtt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["MET"]["mtm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["MET"]["mmt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((30, 240)),
)
ROTAMERS["MET"]["mmm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 30), (240, 360)),
)

ROTAMERS["LYS"]["pppp?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pppt?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pppm?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["pptp?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pptt?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pptm?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["ppmp?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["ppmt?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["ppmm?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["ptpp?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["ptpt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["ptpm?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["pttp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pttt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pttm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["ptmp?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["ptmt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["ptmm?"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["pmpp?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pmpt?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pmpm?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["pmtp?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pmtt?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pmtm?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["pmmp?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["pmmt?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["pmmm?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tppp?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tppt?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tppm?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tptp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tptt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tptm"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tpmp?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tpmt?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tpmm?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["ttpp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["ttpt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["ttpm?"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tttp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tttt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tttm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["ttmp?"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["ttmt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["ttmm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tmpp?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tmpt?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tmpm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tmtp?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tmtt?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tmtm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["tmmp?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["tmmt?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["tmmm?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mppp?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mppt?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mppm?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mptp?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mptt"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mptm?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mpmp?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mpmt?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mpmm?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mtpp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mtpt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mtpm?"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mttp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mttt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mttm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mtmp?"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mtmt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mtmm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mmpp?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mmpt?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mmpm?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mmtp"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mmtt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mmtm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["LYS"]["mmmp?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["LYS"]["mmmt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["LYS"]["mmmm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)

ROTAMERS["ARG"]["ppp_?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["ppt_?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["ppm_?"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["ptp85"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ptp180"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 360)),
)
ROTAMERS["ARG"]["ptt85"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ptt180"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["ptt-85"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["ptm85"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ptm180"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 360)),
)
ROTAMERS["ARG"]["pmp_?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["pmt_?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["pmm_?"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["tpp85"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["tpp180"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 360)),
)
ROTAMERS["ARG"]["tpt85"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["tpt180"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 360)),
)
ROTAMERS["ARG"]["tpm_?"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["ttp85"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ttp180"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["ttp-105"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["ttt85"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ttt180"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["ttt-85"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["ttm105"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["ttm180"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["ttm-85"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["tmp_?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["tmt_?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["tmm_?"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["mpp_?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["mpt_?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["mpm_?"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["mtp85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["mtp180"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["mtp-105"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["mtt85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["mtt180"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["mtt-85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["mtm105"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["mtm180"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["mtm-85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["mmp_?"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 360)),
)
ROTAMERS["ARG"]["mmt85"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
)
ROTAMERS["ARG"]["mmt180"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
)
ROTAMERS["ARG"]["mmt-85"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
)
ROTAMERS["ARG"]["mmm180"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 240)),
)
ROTAMERS["ARG"]["mmm-85"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
)
