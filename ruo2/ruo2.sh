/usr/bin/rsync -e ssh -avzpl \
--include='DONE' --include='CONTCAR' --include='*.traj' --include='pydos' --include='*.csv' --include='*.png' --include='*.json' --include='*/' --exclude='*' \
jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/3_RuO2 .

/usr/bin/rsync -e ssh -avzpl \
--include='*.json' --include='*/' --exclude='*' \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction/0_Hubbard_U/ ./3_RuO2/7_prediction/0_Hubbard_U

/usr/bin/rsync -e ssh -avzpl \
--include='*.json' --include='*/' --exclude='*' \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction/3_surface_Ru/preemptable/ ./3_RuO2/7_prediction/3_surface_Ru

/usr/bin/rsync -e ssh -avzpl \
--include='*.txt' --include='*.lobster' --include='*.json' --include='*/' --exclude='*' \
jiuy97@suncat:/sdf/data/suncat/suncatlustre/hailey/7_prediction/3_surface_Ru/icohp/ ./3_RuO2/7_prediction/5_surface_Ru_icohp
