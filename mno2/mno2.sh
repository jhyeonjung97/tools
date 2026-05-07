#!/bin/bash
/usr/bin/rsync -e ssh -avz \
--include='*.cif' --include='*.vasp' --include='DONE' --include='CONTCAR' --include='*.png' --include='*.csv' --include='*.json' --include='*.traj' --include='*/' --exclude='*' \
jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/4_MnO2 .

/usr/bin/rsync -e "ssh -J jiuy97@s3dflogin.slac.stanford.edu" -avz \
--include='*.cif' --include='*.vasp' --include='DONE' --include='CONTCAR' --include='*.png' --include='*.csv' --include='*.json' --include='*.traj' --include='*/' --exclude='*' \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/4_MnO2 .
