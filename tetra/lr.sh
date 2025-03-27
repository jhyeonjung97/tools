# basic
python ~/bin/tools/tetra/lr.py --Y form --output form

# by row
python ~/bin/tools/tetra/lr.py --Y form --row 3d --output form3d
python ~/bin/tools/tetra/lr.py --Y form --row 4d --output form4d
python ~/bin/tools/tetra/lr.py --Y form --row 5d --output form5d
python ~/bin/tools/tetra/lr.py --Y form --row 4d 5d --output form4d5d

# cohesive
python ~/bin/tools/tetra/lr.py --Y coh --output coh
python ~/bin/tools/tetra/lr.py --Y coh --row 3d --output coh3d
python ~/bin/tools/tetra/lr.py --Y coh --row 4d --output coh4d
python ~/bin/tools/tetra/lr.py --Y coh --row 5d --output coh5d
python ~/bin/tools/tetra/lr.py --Y coh --row 4d 5d --output coh4d5d

# essential descriptors
python ~/bin/tools/tetra/lr.py --Y form --output form