#!/bin/bash

dir='/pscratch/sd/j/jiuy97/6_MNC/0_clean/3d/6_Fe/most_stable/2_'
cd ${dir}
python ~/bin/tools/mnc/add-oo.py
python ~/bin/tools/mnc/add-oo-oh.py
python ~/bin/tools/mnc/add-oo-o.py
python ~/bin/tools/mnc/add-oo-ooh.py

pb_path='/pscratch/sd/j/jiuy97/6_MNC/pourbaix/1_Fe'   
for ads in oo oo-oh oo-o oo-ooh
do
    for spin in LS IS HS
    do
        dest_dir="${pb_path}/${ads}/${spin}1"
        mkdir -p ${dest_dir}
        cp "${dir}/start-${ads}.traj" "${dest_dir}/start.traj"
        cp "~/bin/tools/mnc/submit.sh" "${dest_dir}"
        sed -i -e "s/mnc-sol.py/mnc-sol-${spin}-nupdown-oo.py/g" "${dest_dir}/submit.sh"
        sed -i -e "s/jobname/Fe-${spin}-{ads}/g" "${dest_dir}/submit.sh"
    done
done