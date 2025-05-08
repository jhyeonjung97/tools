# # basic
# python ~/bin/tools/tetra/lr.py --Y form --output form

# # by row
# python ~/bin/tools/tetra/lr.py --Y form --row 3d --output form3d
# python ~/bin/tools/tetra/lr.py --Y form --row 4d --output form4d
# python ~/bin/tools/tetra/lr.py --Y form --row 5d --output form5d
# python ~/bin/tools/tetra/lr.py --Y form --row 4d 5d --output form4d5d

# # cohesive
# python ~/bin/tools/tetra/lr.py --Y coh --output coh
# python ~/bin/tools/tetra/lr.py --Y coh --row 3d --output coh3d
# python ~/bin/tools/tetra/lr.py --Y coh --row 4d --output coh4d
# python ~/bin/tools/tetra/lr.py --Y coh --row 5d --output coh5d
# python ~/bin/tools/tetra/lr.py --Y coh --row 4d 5d --output coh4d5d

# # essential descriptors
# python ~/bin/tools/tetra/lr.py --Y form --output form

python ~/bin/tools/tetra/lr.py --Y form --output form --X 'OS', 'CN', 'numb', 'chg', 'chgn', 'mag', 'volume', 'l_bond', 'madelung', \
'ICOHPm', 'ICOHPmn', 'ICOHPn', 'ICOBIm', 'ICOBImn', 'ICOBIn', 'ICOOPm', 'ICOOPmn', 'ICOOPn', \
'ion-1', 'ion', 'ion+1', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1', \
'pauling', 'Natom', 'mass', 'density', 'Vatom', 'dipole', 'Rcoval', 'Rmetal', 'Rvdw', \
'Tboil', 'Tmelt', 'Hevap', 'Hfus', 'Hform', 'n_electrons', 'd_electrons', \
'base_cfse', 'ee_repulsion', 'jt_effect', 'field_strength', 'cfse', 'exchange_stabilization'

python ~/bin/tools/tetra/lr.py --Y form --output form --X 'OS', 'CN', 'group', 'outer_e', 'Hevap', 'base_cfse', 'chg', 'chgn', 'mag', 'volume', 'l_bond', \
'ICOHPm', 'ICOHPmn', 'ICOHPn', 'ICOBIm', 'ICOBImn', 'ICOBIn', 'ICOOPm', 'ICOOPmn', 'ICOOPn', \
'ion-1', 'ion', 'ion+1', 'ion-1n', 'ionn', 'ion+1n', 'ionN-1', 'ionN', 'ionN+1'