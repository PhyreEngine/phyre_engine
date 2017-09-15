"""
This module contains a minimal mmCIF file and the corresponding metadata
required to convert the mmCIF file into a template.
"""

#: Minimal mmCIF data.
MINIMAL_MMCIF = """\
loop_
_atom_type.symbol
C
N
O
P
S
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   2    C CA    . ALA A 1 4   ? 12.501 39.048 28.539  1.00 30.68 ? 4   ALA A CA    1
ATOM   7    C CA    . TYR A 1 5   A 15.552 39.410 26.282  1.00 8.51  ? 4   TYR A CA    1
ATOM   19   C CA    . ILE A 1 6   ? 17.791 39.281 29.375  1.00 23.27 ? 6   ILE A CA    1
ATOM   27   C CA    . ALA A 1 7   ? 16.004 36.186 30.742  1.00 5.20  ? 7   ALA A CA    1
ATOM   32   C CA    . LYS A 1 8   ? 16.137 34.313 27.425  1.00 4.96  ? 8   LYS A CA    1
ATOM   41   C CA    . GLN A 1 9   ? 19.794 35.327 26.885  1.00 8.15  ? 9   GLN A CA    1
ATOM   50   C CB    . ARG A 1 10  ? 20.706 33.656 30.141  1.00 19.12 ? 10  ARG A CA    1
HETATM 5128 P P     . AMP D 3 .   ? 24.511 18.911 8.472   1.00 27.61 ? 332 AMP A P     1
"""

#: Minimal template in PDB format
MINIMAL_PDB = """\
ATOM     2   CA  ALA A   4      12.501  39.048  28.539  1.00 30.68
ATOM     7   CA  TYR A   5      15.552  39.410  26.282  1.00 8.51
ATOM     19  CA  ILE A   6      17.791  39.281  29.375  1.00 23.27
ATOM     27  CA  ALA A   7      16.004  36.186  30.742  1.00 5.20
ATOM     32  CA  LYS A   8      16.137  34.313  27.425  1.00 4.96
ATOM     41  CA  GLN A   9      19.794  35.327  26.885  1.00 8.15
ATOM     50  CB  ARG A   10     20.706  33.656  30.141  1.00 19.12
HETATM   5128P   AMP D   100    24.511  18.911  8.472   1.00 27.61
"""


#: Original residue IDs of mmCIF file.
ORIG_MAPPING = [
    (' ', 4, ' '),
    (' ', 4, 'A'),
    (' ', 6, ' '),
    (' ', 7, ' '),
    (' ', 8, ' '),
    (' ', 9, ' '),
    (' ', 10, ' '),
    ('H_AMP', 332, ' ')
]

#: Canonical sequence of mmCIF file.
CANONICAL_SEQ = "AYIAKQ"

#: Renumbered indices of residues in canonical sequence.
CANONICAL_SEQ_INDICES = [1, 2, 3, 4, 5, 6]
