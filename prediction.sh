#!/bin/bash
#rsync -avz \
#  --exclude='x/' --include='*.csv' --include='*.png' --include='*.pdf' --include='*/' --exclude='*' \
#  jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction .

rsync -avz \
  --exclude='seed_from_host.py' --exclude='slab-lobster.py' --exclude='bulk-lobster.py' --exclude='run_vasp*.py' --exclude='opt_slab.py' --exclude='opt_bulk*.py' \
  --include='*.png' --include='*.pdf' --include='*.csv' --include='*.py' --include='*/' --exclude='*' \
  jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/7_prediction .
