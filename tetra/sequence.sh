#!/bin/bash

cd ~/Desktop/
sh ~/bin/tools/tetra/vbulk.sh
sh ~/bin/tools/tetra/vslab.sh

python ~/bin/tools/tetra/mendeleev2tsv.py
python ~/bin/tools/tetra/parse_bulk.py
python ~/bin/tools/tetra/parse_bulk_comer.py

python ~/bin/tools/tetra/lr-cfse.py
python ~/bin/tools/tetra/plot_bulk_cfse.py

sh ~/bin/tools/tetra/random.sh
sh ~/bin/tools/tetra/bulk_pred.sh

python ~/bin/tools/tetra/parse_slab.py
python ~/bin/tools/tetra/parse_slab_comer.py