#!/bin/bash

cd ~/Desktop/
sh ./vbulk.sh

python ~/bin/tools/tetra/mendeleev2tsv.py
python ~/bin/tools/tetra/parse_bulk.py
python ~/bin/tools/tetra/parse_bulk_comer.py
python ~/bin/tools/tetra/lr-cfse.py
python ~/bin/tools/tetra/lr-cfse.py --Y form --X group outer_e Hevap base_cfse ICOHPo ionNo OS CN mag volume l_bond chg --output compact
python ~/bin/tools/tetra/lr-cfse.py --Y form --X group outer_e Hform base_cfse ICOHPo ionNo OS CN mag volume l_bond chg --output compact
python ~/bin/tools/tetra/plot_bulk_cfse.py