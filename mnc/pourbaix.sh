scp jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/*_energies.tsv .
scp jiuy97@perlmutter.nersc.gov:/pscratch/sd/j/jiuy97/6_MNC/gas/metals.tsv .

python ~/bin/shoulder/pourbaix-fe.py
python ~/bin/shoulder/pourbaix-co.py
python ~/bin/shoulder/pourbaix-mo.py