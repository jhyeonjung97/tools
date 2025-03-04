#!/bin/bash

/usr/bin/rsync -e ssh --size-only -avzP \
--include='final_with_calculator.json' \
--include='POTCAR' \
--include='MadelungEnergies.lobster' \
--include='GROSSPOP.lobster' \
--include='icohp.txt' \
--include='icobi.txt' \
--include='zpe.txt' \
--include='lattice.cif' \
--include='atoms_bader_charge.json' \
--include='OUTCAR' \
--include='*.tsv' \
--include='*.png' \
--include='*/' \
--exclude='*' \
jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/3_V_bulk .