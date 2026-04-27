/usr/bin/rsync -e ssh -avzP \
--include='*.txt' --include='*.png' --include='*.json' --include='*/' --exclude='*' \
jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/3_RuO2 .

/usr/bin/rsync -e ssh -avzP \
--include='*.txt' --include='*.png' --include='*.json' --include='*/' --exclude='*' \
-e "ssh -J jiuy97@s3dflogin.slac.stanford.edu" \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction/0_Hubbard_U/ ./3_RuO2/7_prediction/0_Hubbard_U

/usr/bin/rsync -e ssh -avzP \
--include='*.txt' --include='*.png' --include='*.json' --include='*/' --exclude='*' \
-e "ssh -J jiuy97@s3dflogin.slac.stanford.edu" \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction/3_surface_Ru/preemptable/ ./3_RuO2/7_prediction/3_surface_Ru
