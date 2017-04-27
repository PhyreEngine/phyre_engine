"""
Rotamer definitions according to MolProbity (and therefore from the Penultimate
Rotamer Library and (Son of) Penutlimate Rotamer Library.

This module exposes the `ROTAMERS` variable, which is a dictionary indexed by
three-letter amino acid codes. Each key points to another dictionary, indexed by
the name of the rotamer according to MolProbity nomenclature. The values of this
dictionary are a tuple of `AngleRange`s indicating the actual range of each chi
angle for this rotamer.

>>> ROTAMERS["ARG"]["ppp80"]
"""

from ..angle_range import AngleRange
from .generic import NUM_CHI_ANGLES

ROTAMERS = {aa: {} for aa, n in NUM_CHI_ANGLES.items() if n > 0}
ROTAMERS["ARG"]["ppp80"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["ppp-140"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((130, 360))
)
ROTAMERS["ARG"]["ppt90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 125)),
    AngleRange((120, 240)),
    AngleRange((0, 125))
)
ROTAMERS["ARG"]["ppt170"] = (
    AngleRange((0, 120)),
    AngleRange((0, 130)),
    AngleRange((120, 240)),
    AngleRange((125, 230))
)
ROTAMERS["ARG"]["ppt-90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 125)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ptp90"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["ptp-170"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((130, 230))
)
ROTAMERS["ARG"]["ptp-110"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ptt90"] = (
    AngleRange((0, 120)),
    AngleRange((125, 250)),
    AngleRange((120, 240)),
    AngleRange((0, 125))
)
ROTAMERS["ARG"]["ptt180"] = (
    AngleRange((0, 120)),
    AngleRange((130, 240)),
    AngleRange((120, 240)),
    AngleRange((125, 230))
)
ROTAMERS["ARG"]["ptt-90"] = (
    AngleRange((0, 120)),
    AngleRange((125, 250)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ptm-80"] = (
    AngleRange((0, 120)),
    AngleRange((105, 240)),
    AngleRange((240, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ptm160"] = (
    AngleRange((0, 120)),
    AngleRange((105, 250)),
    AngleRange((240, 360)),
    AngleRange((0, 230))
)
ROTAMERS["ARG"]["pmt-80"] = (
    AngleRange((0, 120)),
    AngleRange((250, 360)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["pmt170"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 250)),
    AngleRange((135, 230))
)
ROTAMERS["ARG"]["pmt100"] = (
    AngleRange((0, 120)),
    AngleRange((250, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 135))
)
ROTAMERS["ARG"]["pmm-80"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["pmm150"] = (
    AngleRange((0, 120)),
    AngleRange((250, 360)),
    AngleRange((250, 360)),
    AngleRange((0, 230))
)
ROTAMERS["ARG"]["tpp80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["tpp-160"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((0, 120)),
    AngleRange((130, 360))
)
ROTAMERS["ARG"]["tpt90"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["ARG"]["tpt170"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((120, 240)),
    AngleRange((120, 230))
)
ROTAMERS["ARG"]["tpt-90"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["tpm170"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((240, 360)),
    AngleRange((0, 230))
)
ROTAMERS["ARG"]["tpm-80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 110)),
    AngleRange((240, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ttp80"] = (
    AngleRange((120, 240)),
    AngleRange((110, 260)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["ttp-170"] = (
    AngleRange((120, 240)),
    AngleRange((110, 250)),
    AngleRange((0, 120)),
    AngleRange((130, 225))
)
ROTAMERS["ARG"]["ttp-110"] = (
    AngleRange((120, 240)),
    AngleRange((110, 250)),
    AngleRange((0, 120)),
    AngleRange((225, 360))
)
ROTAMERS["ARG"]["ttt90"] = (
    AngleRange((120, 240)),
    AngleRange((110, 235)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["ARG"]["ttt180"] = (
    AngleRange((120, 240)),
    AngleRange((110, 235)),
    AngleRange((120, 240)),
    AngleRange((120, 230))
)
ROTAMERS["ARG"]["ttt-90"] = (
    AngleRange((120, 240)),
    AngleRange((110, 235)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["ttm110"] = (
    AngleRange((120, 240)),
    AngleRange((110, 240)),
    AngleRange((235, 360)),
    AngleRange((0, 135))
)
ROTAMERS["ARG"]["ttm170"] = (
    AngleRange((120, 240)),
    AngleRange((110, 240)),
    AngleRange((235, 360)),
    AngleRange((135, 230))
)
ROTAMERS["ARG"]["ttm-80"] = (
    AngleRange((120, 240)),
    AngleRange((110, 240)),
    AngleRange((235, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["tmt90"] = (
    AngleRange((120, 240)),
    AngleRange((235, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["tmt170"] = (
    AngleRange((120, 240)),
    AngleRange((235, 360)),
    AngleRange((120, 240)),
    AngleRange((130, 240))
)
ROTAMERS["ARG"]["tmt-80"] = (
    AngleRange((120, 240)),
    AngleRange((235, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["ARG"]["tmm160"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((235, 360)),
    AngleRange((0, 230))
)
ROTAMERS["ARG"]["tmm-80"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((235, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["mmm-85"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((225, 360))
)
ROTAMERS["ARG"]["mmm160"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 225))
)
ROTAMERS["ARG"]["mmp-170"] = (
    AngleRange((240, 360)),
    AngleRange((250, 360)),
    AngleRange((0, 120)),
    AngleRange((130, 360))
)
ROTAMERS["ARG"]["mmp80"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["mmt-90"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["mmt180"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((130, 230))
)
ROTAMERS["ARG"]["mmt90"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["mpp-170"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((130, 360))
)
ROTAMERS["ARG"]["mpp80"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["ARG"]["mpt-90"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["mpt180"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 230))
)
ROTAMERS["ARG"]["mpt90"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["ARG"]["mtm-85"] = (
    AngleRange((240, 360)),
    AngleRange((110, 240)),
    AngleRange((240, 360)),
    AngleRange((225, 360))
)
ROTAMERS["ARG"]["mtm110"] = (
    AngleRange((240, 360)),
    AngleRange((110, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 135))
)
ROTAMERS["ARG"]["mtm180"] = (
    AngleRange((240, 360)),
    AngleRange((110, 240)),
    AngleRange((240, 360)),
    AngleRange((135, 225))
)
ROTAMERS["ARG"]["mtp-110"] = (
    AngleRange((240, 360)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["mtp180"] = (
    AngleRange((240, 360)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((135, 230))
)
ROTAMERS["ARG"]["mtp85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 135))
)
ROTAMERS["ARG"]["mtt-85"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["ARG"]["mtt180"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((130, 230))
)
ROTAMERS["ARG"]["mtt90"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 130))
)
ROTAMERS["VAL"]["p"] = (
    AngleRange((0, 120)),
)
ROTAMERS["VAL"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["VAL"]["m"] = (
    AngleRange((240, 360)),
)
ROTAMERS["PRO"]["Cg_endo"] = (
    AngleRange((0, 180)),
    AngleRange((0, 360)),
    AngleRange((0, 360))
)
ROTAMERS["PRO"]["Cg_exo"] = (
    AngleRange((180, 360)),
    AngleRange((0, 360)),
    AngleRange((0, 360))
)
ROTAMERS["PHE"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["PHE"]["t80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180))
)
ROTAMERS["PHE"]["m-80"] = (
    AngleRange((240, 360)),
    AngleRange((40, 140))
)
ROTAMERS["PHE"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((0, 40))
)
ROTAMERS["PHE"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((140, 180))
)
ROTAMERS["THR"]["p"] = (
    AngleRange((0, 120)),
)
ROTAMERS["THR"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["THR"]["m"] = (
    AngleRange((240, 360)),
)
ROTAMERS["TRP"]["p-90"] = (
    AngleRange((0, 120)),
    AngleRange((180, 360))
)
ROTAMERS["TRP"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["TRP"]["t-100"] = (
    AngleRange((120, 240)),
    AngleRange((180, 315))
)
ROTAMERS["TRP"]["t60"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180))
)
ROTAMERS["TRP"]["t60"] = (
    AngleRange((120, 240)),
    AngleRange((315, 360))
)
ROTAMERS["TRP"]["m-90"] = (
    AngleRange((240, 360)),
    AngleRange((180, 305))
)
ROTAMERS["TRP"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((305, 360))
)
ROTAMERS["TRP"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((0, 45))
)
ROTAMERS["TRP"]["m100"] = (
    AngleRange((240, 360)),
    AngleRange((45, 180))
)
ROTAMERS["LYS"]["pmmt"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["pmtt"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["pptt"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["ptmm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["ptmt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["ptpp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["ptpt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["pttm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["pttp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["pttt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["tmmm"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["tmmt"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["tmtm"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["tmtp"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["tmtt"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["tppp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["tppt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["tptm"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["tptp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["tptt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 260)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["ttmm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["ttmp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["ttmt"] = (
    AngleRange((120, 240)),
    AngleRange((110, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["ttpm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["ttpp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 250)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["ttpt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 260)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["tttm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["tttp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["tttt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mmmm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mmmt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mmpt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mmtm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mmtp"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["mmtt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mppt"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mptm"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mptp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["mptt"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mtmm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mtmp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["mtmt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mtpm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mtpp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["mtpt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LYS"]["mttm"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LYS"]["mttp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LYS"]["mttt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["MSE"]["pp-130"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 360))
)
ROTAMERS["MSE"]["ppp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["MSE"]["ptp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120))
)
ROTAMERS["MSE"]["ptt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["MSE"]["ptm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["MSE"]["pmt"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["MSE"]["pmm"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["MSE"]["tpp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["MSE"]["tpt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((130, 330))
)
ROTAMERS["MSE"]["ttp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["MSE"]["ttt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 230))
)
ROTAMERS["MSE"]["ttm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 235)),
    AngleRange((230, 360))
)
ROTAMERS["MSE"]["tmt"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((110, 230))
)
ROTAMERS["MSE"]["tmm"] = (
    AngleRange((120, 240)),
    AngleRange((235, 360)),
    AngleRange((230, 360))
)
ROTAMERS["MSE"]["mpp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["MSE"]["mpt"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 220))
)
ROTAMERS["MSE"]["mpm"] = (
    AngleRange((240, 360)),
    AngleRange((0, 115)),
    AngleRange((220, 360))
)
ROTAMERS["MSE"]["mtp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["MSE"]["mtt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["MSE"]["mtm"] = (
    AngleRange((240, 360)),
    AngleRange((115, 240)),
    AngleRange((240, 360))
)
ROTAMERS["MSE"]["mmp"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 130))
)
ROTAMERS["MSE"]["mmt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((130, 230))
)
ROTAMERS["MSE"]["mmm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((230, 360))
)
ROTAMERS["MET"]["pp-130"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((120, 360))
)
ROTAMERS["MET"]["ppp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["MET"]["ptp"] = (
    AngleRange((0, 120)),
    AngleRange((120, 250)),
    AngleRange((0, 120))
)
ROTAMERS["MET"]["ptt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["MET"]["ptm"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["MET"]["pmt"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["MET"]["pmm"] = (
    AngleRange((0, 120)),
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["MET"]["tpp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 130))
)
ROTAMERS["MET"]["tpt"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((130, 330))
)
ROTAMERS["MET"]["ttp"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["MET"]["ttt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240)),
    AngleRange((120, 230))
)
ROTAMERS["MET"]["ttm"] = (
    AngleRange((120, 240)),
    AngleRange((120, 235)),
    AngleRange((230, 360))
)
ROTAMERS["MET"]["tmt"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360)),
    AngleRange((110, 230))
)
ROTAMERS["MET"]["tmm"] = (
    AngleRange((120, 240)),
    AngleRange((235, 360)),
    AngleRange((230, 360))
)
ROTAMERS["MET"]["mpp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["MET"]["mpt"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((120, 220))
)
ROTAMERS["MET"]["mpm"] = (
    AngleRange((240, 360)),
    AngleRange((0, 115)),
    AngleRange((220, 360))
)
ROTAMERS["MET"]["mtp"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["MET"]["mtt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["MET"]["mtm"] = (
    AngleRange((240, 360)),
    AngleRange((115, 240)),
    AngleRange((240, 360))
)
ROTAMERS["MET"]["mmp"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 130))
)
ROTAMERS["MET"]["mmt"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((130, 230))
)
ROTAMERS["MET"]["mmm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360)),
    AngleRange((230, 360))
)
ROTAMERS["ASN"]["p0"] = (
    AngleRange((0, 120)),
    AngleRange((0, 360))
)
ROTAMERS["ASN"]["t0"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["ASN"]["t0"] = (
    AngleRange((120, 240)),
    AngleRange((180, 360))
)
ROTAMERS["ASN"]["t160"] = (
    AngleRange((120, 240)),
    AngleRange((120, 180))
)
ROTAMERS["ASN"]["m-40"] = (
    AngleRange((240, 360)),
    AngleRange((220, 360))
)
ROTAMERS["ASN"]["m-40"] = (
    AngleRange((240, 360)),
    AngleRange((0, 50))
)
ROTAMERS["ASN"]["m110"] = (
    AngleRange((240, 360)),
    AngleRange((50, 220))
)
ROTAMERS["LEU"]["pp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["LEU"]["pt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["LEU"]["tp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["LEU"]["tt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["LEU"]["tm"] = (
    AngleRange((120, 240)),
    AngleRange((240, 360))
)
ROTAMERS["LEU"]["mp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120))
)
ROTAMERS["LEU"]["mt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["LEU"]["mm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["ASP"]["p0"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["ASP"]["t0"] = (
    AngleRange((120, 240)),
    AngleRange((0, 45))
)
ROTAMERS["ASP"]["t0"] = (
    AngleRange((120, 240)),
    AngleRange((120, 180))
)
ROTAMERS["ASP"]["t70"] = (
    AngleRange((120, 240)),
    AngleRange((45, 120))
)
ROTAMERS["ASP"]["m-30"] = (
    AngleRange((240, 360)),
    AngleRange((0, 180))
)
ROTAMERS["GLN"]["pt0"] = (
    AngleRange((0, 120)),
    AngleRange((130, 230)),
    AngleRange((0, 360))
)
ROTAMERS["GLN"]["pm20"] = (
    AngleRange((0, 120)),
    AngleRange((230, 360)),
    AngleRange((0, 360))
)
ROTAMERS["GLN"]["pp30"] = (
    AngleRange((0, 120)),
    AngleRange((0, 130)),
    AngleRange((0, 360))
)
ROTAMERS["GLN"]["tp-100"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((150, 310))
)
ROTAMERS["GLN"]["tp40"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 150))
)
ROTAMERS["GLN"]["tp40"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((310, 360))
)
ROTAMERS["GLN"]["tt0"] = (
    AngleRange((120, 240)),
    AngleRange((120, 230)),
    AngleRange((0, 360))
)
ROTAMERS["GLN"]["tm-30"] = (
    AngleRange((120, 230)),
    AngleRange((230, 360)),
    AngleRange((0, 80))
)
ROTAMERS["GLN"]["tm-30"] = (
    AngleRange((120, 230)),
    AngleRange((230, 360)),
    AngleRange((210, 360))
)
ROTAMERS["GLN"]["mp10"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((280, 360))
)
ROTAMERS["GLN"]["mp10"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 150))
)
ROTAMERS["GLN"]["mt0"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 360))
)
ROTAMERS["GLN"]["mm-40"] = (
    AngleRange((230, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 50))
)
ROTAMERS["GLN"]["mm-40"] = (
    AngleRange((230, 360)),
    AngleRange((240, 360)),
    AngleRange((210, 360))
)
ROTAMERS["GLN"]["mm110"] = (
    AngleRange((230, 360)),
    AngleRange((240, 360)),
    AngleRange((50, 210))
)
ROTAMERS["GLN"]["tm130"] = (
    AngleRange((120, 230)),
    AngleRange((230, 360)),
    AngleRange((80, 210))
)
ROTAMERS["GLN"]["mp-120"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((150, 280))
)
ROTAMERS["GLU"]["pp20"] = (
    AngleRange((0, 120)),
    AngleRange((0, 130)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["pt0"] = (
    AngleRange((0, 120)),
    AngleRange((130, 230)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["pm20"] = (
    AngleRange((0, 120)),
    AngleRange((230, 360)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["tp30"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["tt0"] = (
    AngleRange((120, 240)),
    AngleRange((120, 230)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["tm-30"] = (
    AngleRange((120, 230)),
    AngleRange((230, 360)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["mp0"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["mt-10"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240)),
    AngleRange((0, 180))
)
ROTAMERS["GLU"]["mm-30"] = (
    AngleRange((230, 360)),
    AngleRange((240, 360)),
    AngleRange((0, 180))
)
ROTAMERS["HIS"]["p-80"] = (
    AngleRange((0, 120)),
    AngleRange((200, 360))
)
ROTAMERS["HIS"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 200))
)
ROTAMERS["HIS"]["t-170"] = (
    AngleRange((120, 240)),
    AngleRange((140, 230))
)
ROTAMERS["HIS"]["t-90"] = (
    AngleRange((120, 240)),
    AngleRange((230, 360))
)
ROTAMERS["HIS"]["t70"] = (
    AngleRange((120, 240)),
    AngleRange((0, 140))
)
ROTAMERS["HIS"]["m-70"] = (
    AngleRange((240, 360)),
    AngleRange((225, 360))
)
ROTAMERS["HIS"]["m-70"] = (
    AngleRange((240, 360)),
    AngleRange((0, 20))
)
ROTAMERS["HIS"]["m170"] = (
    AngleRange((240, 360)),
    AngleRange((130, 225))
)
ROTAMERS["HIS"]["m90"] = (
    AngleRange((240, 360)),
    AngleRange((20, 130))
)
ROTAMERS["ILE"]["pp"] = (
    AngleRange((0, 120)),
    AngleRange((0, 120))
)
ROTAMERS["ILE"]["pt"] = (
    AngleRange((0, 120)),
    AngleRange((120, 240))
)
ROTAMERS["ILE"]["tp"] = (
    AngleRange((120, 240)),
    AngleRange((0, 120))
)
ROTAMERS["ILE"]["tt"] = (
    AngleRange((120, 240)),
    AngleRange((120, 240))
)
ROTAMERS["ILE"]["mp"] = (
    AngleRange((240, 360)),
    AngleRange((0, 120))
)
ROTAMERS["ILE"]["mt"] = (
    AngleRange((240, 360)),
    AngleRange((120, 240))
)
ROTAMERS["ILE"]["mm"] = (
    AngleRange((240, 360)),
    AngleRange((240, 360))
)
ROTAMERS["TYR"]["p90"] = (
    AngleRange((0, 120)),
    AngleRange((0, 180))
)
ROTAMERS["TYR"]["t80"] = (
    AngleRange((120, 240)),
    AngleRange((0, 180))
)
ROTAMERS["TYR"]["m-80"] = (
    AngleRange((240, 360)),
    AngleRange((40, 140))
)
ROTAMERS["TYR"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((0, 40))
)
ROTAMERS["TYR"]["m-10"] = (
    AngleRange((240, 360)),
    AngleRange((140, 180))
)
ROTAMERS["SER"]["p"] = (
    AngleRange((0, 120)),
)
ROTAMERS["SER"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["SER"]["m"] = (
    AngleRange((240, 360)),
)
ROTAMERS["CYS"]["p"] = (
    AngleRange((0, 120)),
)
ROTAMERS["CYS"]["t"] = (
    AngleRange((120, 240)),
)
ROTAMERS["CYS"]["m"] = (
    AngleRange((240, 360)),
)
