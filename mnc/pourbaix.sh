#!/bin/bash

dir='/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/6_Fe/most_stable/relaxed'
cd ${dir}
python ~/bin/tools/mnc/add-oo.py
python ~/bin/tools/mnc/add-oo-oo.py
python ~/bin/tools/mnc/add-oo-ohh.py
python ~/bin/tools/mnc/add-oo-oh.py
python ~/bin/tools/mnc/add-oo-o.py
python ~/bin/tools/mnc/add-oo-ooh.py
python ~/bin/tools/mnc/add-ohh.py
python ~/bin/tools/mnc/add-ohh-ohh.py
python ~/bin/tools/mnc/add-ohh-oo.py
python ~/bin/tools/mnc/add-ohh-oh.py
python ~/bin/tools/mnc/add-ohh-o.py
python ~/bin/tools/mnc/add-ohh-ooh.py

pb_path='/pscratch/sd/j/jiuy97/6_MNC/pourbaix/2_Co'   
for ads in oo oo-oh oo-o oo-ooh
do
    for spin in LS IS HS
    do
        dest_dir="${pb_path}/${ads}/${spin}1"
        # mkdir -p ${dest_dir}
        cd ${dest_dir}
        cp "${dir}/start-${ads}.traj" ./start.traj
        cp ~/bin/tools/mnc/submit.sh .
        sed -i -e "s/mnc-sol.py/mnc-sol-${spin}-nupdown-oo.py/g" submit.sh
        sed -i -e "s/jobname/Co-${ads}-${spin}/g" submit.sh
        # sbatch submit.sh
        python "/pscratch/sd/j/jiuy97/6_MNC/scripts/mnc-sol-${spin}-nupdown-oo.py"
    done
done