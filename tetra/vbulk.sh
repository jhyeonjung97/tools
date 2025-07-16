/usr/bin/rsync -e ssh -avzP \
--include='DONE' \
--include='POSCAR' \
--include='CONTCAR' \
--include='submit.sh' \
--include='unmatched' \
--include='*.json' --include='*.lobster' \
--include='*.txt' --include='*.png' \
--include='*/' --exclude='*' \
jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/7_V_bulk .
